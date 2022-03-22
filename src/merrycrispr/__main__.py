#!/usr/bin/env python3

from enum import Enum
from functools import partial
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pyfaidx
import typer
# from pathos.multiprocessing import ProcessingPool as Pool
from p_tqdm import p_umap
from pkg_resources import resource_filename
from tqdm.auto import tqdm

from . import merrycrispr
from .find_spacers import find_spacers
from .library_assembly import assemble_library, assemble_paired_library
from .off_target_scoring import off_target_discovery, off_target_scoring
from .on_target_scoring import on_target_scoring
from .seqextractor import (display_gtf_features, display_gtf_geneids,
                           display_gtf_genes, extract,
                           extract_for_tss_adjacent)
from .species_getter import (available_species, build_bowtie_index,
                             get_resources)


class LibraryType(str, Enum):
    knockout = "knockout"
    repressor = "repressor"
    activator = "activator"
    excision = "excision"
    cas13 = "Cas13"


@merrycrispr.command("prep_sequences")
def prep_sequences(
    gtf: Path = typer.Argument(..., "--gtf", "-g", help="Input GTF/GFF file"),
    library_type: LibraryType = typer.Option(
        LibraryType.knockout, "--library_type", "-l", help="Target library type"
    ),
    fasta: Path = typer.Argument(
        ...,
        "--fasta",
        "-f",
        help="FASTA sequence file",
    ),
    output: Optional[Path] = typer.Option(
        None, "--output", "-o", help="Location to write output file"
    ),
    feature_type: Optional[str] = typer.Option(
        None, "--feature_type", "-t", help="Feature type"
    ),
    gene_name: Optional[str] = typer.Option(
        None,
        "--gene_name",
        "-n",
        help="Gene(s) to extract.  To extract multiple genes, enclose the entire list in quotation marks and leave a space between each gene.",
    ),
    bound: Optional[int] = typer.Option(
        0,
        "--bound",
        help="Retrieve a given number of bases on either side of the feature instead of the sequence corresponding to a feature",
    ),
    show_features: Optional[bool] = typer.Option(
        False,
        "--show_features",
        help="Scan an annotation file to identify the features present",
    ),
    show_genes: Optional[bool] = typer.Option(
        False,
        "--show_genes",
        help="Scan an annotation file to identify the genes present",
    ),
    show_geneids: Optional[bool] = typer.Option(
        False,
        "--show_geneids",
        help="Scan an annotation file to identify the geneIDs present",
    ),
) -> None:
    """Generate target sequences to search for spacers
    \f
    Generally, a utility for extracting sequences from a FASTA file for a given GTF annotation.
    \f
    Parameters
    ----------
    gtf : `str`
        GTF/GFF file that matches the input FASTA file. Preferably one from Ensembl/GENCODE.
        Gzipped GTF/GFF files are acceptable, though their use may impose a performance penality.
    library_type : `str`, optional
        Target library type.  If provided, other parameters will be set to generate an appropriate
        set of features.
    fasta : `str`
        FASTA file corresponding to the input GTF/GFF.  Preferably one from Ensembl/GENCODE. Gzipped
        FASTA files will work, but their use results in a very substantial performance loss.
    output : `str`, optional (default: `None`)
        File to write extracted sequences to.  If not provided, seqextractor will instead return a
        `list` of :class:`~pyfaidx.Sequence` objects.
    feature_type : (default: `None`)
        A list of target feature types (i.e. gene, exon, intron, etc...).  Multiple types are
        allowed.
    gene_name : `str`, optional (default: `None`)
        If provided, seqextractor will only extract sequences for that list of genes.
    bound : `int`, optional (default: None)
        Instead of extracting the sequences corresponding to the feature itself, extract an area
        equal to the length of `bound` on either side of the feature.
    show_features : `bool` (default: `False`)
        Instead of extracting anything, display a list of the feature types present in the GTF/GFF.
    show_genes : `bool` (default: `False`)
        Instead of extracting anything, display a list of the genes present in the GTF/GFF.
    show_geneids : `bool` (default: `False`)
        Instead of extracting anything, display a list of the geneids present in the GTF/GFF.

    Returns
    -------
    `None`

    """
    if isinstance(gene_name, str):
        gene_name = gene_name.split()
    if library_type:
        if library_type == "knockout":
            feature_type = "exon"
            extract(
                gtffile=gtf,
                fastafile=fasta,
                feature_type=feature_type,
                outfile=output,
                gene_name=gene_name,
            )
        elif library_type == "repressor" or library_type == "activator":
            if bound is None:
                bound = 100
            extract_for_tss_adjacent(
                gtffile=gtf,
                fastafile=fasta,
                outfile=output,
                gene_name=gene_name,
                boundary=bound,
            )
        elif library_type == "excision":
            if bound is None:
                bound = 100
            extract(
                gtffile=gtf,
                fastafile=fasta,
                feature_type=feature_type,
                outfile=output,
                gene_name=gene_name,
                boundary=bound,
            )
        elif library_type == "Cas13":
            extract(
                gtffile=gtf,
                fastafile=fasta,
                feature_type="CDS",
                outfile=output,
                gene_name=gene_name,
            )

    elif show_features:
        if gtf:
            display_gtf_features(gtf)
        else:
            print(
                "Please enter the name of the GFF which you wish to scan for features"
            )
    elif show_genes:
        if gtf:
            display_gtf_genes(gtf, feature_type)
        else:
            print(
                "Please enter the name of the GFF which you wish to scan for features"
            )
    elif show_geneids:
        if gtf:
            display_gtf_geneids(gtf, feature_type)
        else:
            print(
                "Please enter the name of the GFF which you wish to scan for features"
            )
    elif not output:
        print("Please enter an output file name")
    elif not fasta:
        print(
            "Please enter the name of the file containing matching sequences (in FASTA format)"
        )
    elif not gtf:
        print("Please enter the name of the file containing features (in GFF format)")
    else:
        extract(
            gtffile=gtf,
            fastafile=fasta,
            feature_type=feature_type,
            outfile=output,
            gene_name=gene_name,
            boundary=bound,
        )


DATA_PATH = resource_filename("merrycrispr", "data/")
NUCLEASES = pd.read_csv(
    f"{DATA_PATH}nuclease_list.csv",
    dtype={
        "nuclease": str,
        "pam": str,
        "spacer_regex": str,
        "start": np.int8,
        "end": np.int8,
        "strand": str,
        "on_target_rule": str,
        "off_target_rule": str,
    },
    skip_blank_lines=True,
)

AVAILABLE_NUCLEASES = ", ".join(NUCLEASES["nuclease"])


@merrycrispr.command()
def create_library(
    input_sequences: Path = typer.Argument(
        ...,
        "--input_sequences",
        "-i",
        help="Input FASTA file containing sequences to target.",
    ),
    output_library: Path = typer.Argument(
        ...,
        "--output_library",
        "-p",
        help="Name of file to write library to (in CSV format).",
    ),
    reference: Path = typer.Argument(
        ...,
        "--reference",
        "-r",
        help="Path to the directory containing the appropriate Bowtie index.",
    ),
    restriction_sites: str = typer.Option(None),
    largeindex: bool = typer.Option(False, "--largeindex", help=""),
    on_target_rule_set: Optional[str] = typer.Option(
        None,
        "--on_target_rule_set",
        help="Override nuclease on-target rule set defined in nuclease.csv. Current options include '1', '2', and 'Azimuth'.",
    ),
    on_target_score_threshold: int = typer.Option(
        0,
        "--on_target_score_threshold",
        "-on",
        help="Spacers with an on-target score below this will be ignored.",
    ),
    off_target_rule_set: Optional[str] = typer.Option(
        None,
        "--off_target_rule_set",
        help="Override nuclease off-target rule set defined in nuclease.csv. Currently 'Hsu' is the only option",
    ),
    off_target_score_threshold: int = typer.Option(
        0,
        "--off_target_score_threshold",
        "-off",
        help="Spacers with an off-target score below this will be ignored.",
    ),
    off_target_count_threshold: int = typer.Option(
        100,
        "--off_target_count_threshold",
        help="Spacers with more than this many off-targets will be ignored.",
    ),
    number_mismatches_to_consider: int = typer.Option(
        3,
        "--number_mismatches_to_consider",
        help=(
            "Number of mismatches to allow in potential off-targets.  A number "
            "between 0 and 3.  Without setting `off_target_count_threshold` "
            "appropriately, higher values passed to this option may result in "
            "extremely long run times."
        ),
    ),
    nuclease: str = typer.Option(
        "SpCas9",
        "--nuclease",
        "-n",
        help=f"Cas nuclease to design for. "
        f"Current options include {AVAILABLE_NUCLEASES}",
    ),
    spacers_per_feature: int = typer.Option(
        9,
        "--spacers_per_feature",
        help=(
            "Number of spacers to find for each feature. Use '0' to return "
            "all spacers."
        ),
    ),
    reject: bool = typer.Option(False, "--reject", help=""),
    paired: bool = typer.Option(
        False,
        "--paired",
        help="Should spacers be designed to work as pairs (e.g. for excision)?",
    ),
    number_upstream_spacers: int = typer.Option(
        0,
        "--number_upstream_spacers",
        help=(
            "If designing paired spacers, number of spacers to design that target "
            "upstream of the feature."
        ),
    ),
    number_downstream_spacers: int = typer.Option(
        0,
        "--number_downstream_spacers",
        help=(
            "If designing paired spacers, number of spacers to "
            "design that target downstream of the feature."
        ),
    ),
    cores: int = typer.Option(
        0,
        "--cores",
        "-c",
        help="Number of processors to use. By default, will use all available.",
    ),
    chunks: int = typer.Option(
        8, "--chunks", "-k", help="Number of parts to split internal dataframes into."
    ),
    verbose: bool = typer.Option(
        False, "--verbose", "-v", help="Display messages indicative of progress."
    ),
    write_early_exit: bool = typer.Option(
        False,
        "--write_early_exit",
        help="Write the pre-on-target scoring spacer dataframe to a file for debuging purposes",
    ),
) -> None:
    """Build a CRISPR library
    \f

    Parameters
    ----------
    :param input_sequences :
    :param output_library :
    :param reference :
    :param restriction_sites :
    :param largeindex :
    :param on_target_score_threshold :
    :param off_target_score_threshold :
    :param off_target_count_threshold :
    number_mismatches_to_consider
    :param nuclease :
    :param spacers_per_feature :
    :param reject :
    :param paired :
    :param rule_set :
    :param number_upstream_spacers :
    :param number_downstream_spacers :
    :param cores :
    :param chunks:
    verbose : bool

    Return
    ------
    :type reference: object
    """
    targets = pyfaidx.Fasta(input_sequences)

    global NUCLEASES
    nuc = NUCLEASES[NUCLEASES["nuclease"] == nuclease].to_dict(orient="records")[0]

    spacers_df = find_spacers(
        itemlist=targets,
        nuclease_info=nuc,
        restriction_sites=restriction_sites,
        chunks=chunks,
    )
    if write_early_exit:
        spacers_df.to_csv("/Users/milessmith/workspace/mc_human_files/early_exit.csv")
        return None
    initialnumber = spacers_df.shape[0]

    # thank the gods for the tutorial at
    # https://www.machinelearningplus.com/python/parallel-processing-python/
    # scoring_pool = Pool(cores)
    chunked_spacer_dfs = np.array_split(spacers_df, chunks * 10)

    scoring_partial = partial(
        on_target_scoring,
        rule_set=on_target_rule_set,
        on_target_score_threshold=on_target_score_threshold,
    )

    spacers_df = pd.concat(p_umap(scoring_partial, chunked_spacer_dfs))

    # scoring_pool.close()
    # scoring_pool.join()
    # scoring_pool.clear()

    if spacers_df.shape[0] == 0:
        print("Sorry, no spacers matching that criteria were found")
        exit()
    else:
        if verbose:
            print(
                f"Finished scoring spacers. {spacers_df.shape[0]} of {initialnumber} "
                f"spacers have an on-target score above the cutoff threshold of "
                f"{on_target_score_threshold}."
            )

    tqdm.pandas(desc="Adding tracking hashes", unit="spacers")
    spacers_df["hash"] = spacers_df.progress_apply(lambda x: hash(tuple(x)), axis=1)
    if verbose:
        print("\nBeginning Bowtie alignment...")
    off_target_results_file = off_target_discovery(
        spacers_df=spacers_df,
        nuclease_info=nuc,
        cpus=cores,
        refgenome=reference,
        large_index_size=largeindex,
        reject=reject,
        number_mismatches_to_consider=number_mismatches_to_consider,
        verbose=verbose,
    )

    spacers_df = off_target_scoring(
        otrf=off_target_results_file,
        spacers_df=spacers_df,
        nuclease_info=nuc,
        rule_set=off_target_rule_set,
        off_target_score_threshold=off_target_score_threshold,
        off_target_count_threshold=off_target_count_threshold,
        verbose=verbose,
    )

    if paired:
        guide_library = assemble_paired_library(
            spacers=spacers_df,
            on_target_score_threshold=on_target_score_threshold,
            off_target_score_threshold=off_target_score_threshold,
            number_upstream_spacers=number_upstream_spacers,
            number_downstream_spacers=number_downstream_spacers,
        )
    else:
        guide_library = assemble_library(
            spacers=spacers_df,
            on_target_score_threshold=on_target_score_threshold,
            off_target_score_threshold=off_target_score_threshold,
            spacers_per_feature=spacers_per_feature,
        )
    guide_library.to_csv(output_library)
    print("Finished.")


class SpeciesAttribute(str, Enum):
    taxon_id = "taxon_id"
    accession = "accession"
    aliases = "aliases"
    division = "division"
    groups = "groups"
    release = "release"
    name = "name"
    strain = "strain"
    strain_collection = "strain_collection"
    display_name = "display_name"
    assembly = "assembly"
    common_name = "common_name"


@merrycrispr.command()
def new_species(
    species_value: str = typer.Argument(
        ...,
        "--species_value",
        "-s",
        help=(
            "Species value to retrieve sequences for.  Must be a value for the "
            "attribute chosen, i.e. 'canis_familiaris' if using the 'name' attribute"
        ),
    ),
    species_attribute: SpeciesAttribute = typer.Argument(
        ...,
        "--species_attribute",
        "-a",
        help=("Species attribute to search when retrieving sequences."),
    ),
    dest: Optional[Path] = typer.Option(
        None,
        "--dest",
        "-d",
        help=(
            "Location in which to place downloaded files and, if selected, new "
            "Bowtie index."
        ),
    ),
    build_bowtie: bool = typer.Option(
        False,
        "--build_bowtie",
        "-b",
        help="Build a bowtie index for the retrieved species",
    ),
    show_available_builds: bool = typer.Option(
        False,
        "--show_available_builds",
        "-w",
        help="Display a list of the species available from Ensembl",
    ),
):
    """Import a new species from Ensembl.  Optionally, build a Bowtie index for the files.
    \f
    Parameters
    -----------
    species_value : `str`
        Species value to retrieve sequences for.  Must be a value for the attribute chosen, i.e.
         'canis_familiaris' if using the 'name' attribute
    species_attribute : `str`
        Species attribute to search when retrieving sequences. Acceptable values include
        'taxon_id', 'accession', 'aliases', 'division', 'groups', 'release', 'name', 'strain',
        'strain_collection', 'display_name', 'assembly', and 'common_name'
    dest : `dest`, optional
        Directory in which to store files.  If not provided, files will be placed in a new
        subdirectory called `mc_resources`.
    build_bowtie : `bool`
        Should a Bowtie index be built?
    show_available_builds : `bool`
        Show the species currently available from Ensembl

    Return
    -------
    """
    message = ""
    if show_available_builds:
        available_species()

    else:
        if not dest:
            dest = Path().cwd() / "mc_resources"
        if not dest.exists():
            dest.mkdir()
        gtf, fasta = get_resources(
            species_value=species_value,
            species_attribute=species_attribute,
            resource_folder=dest,
        )
        message += f"GTF located at {gtf}, FASTA located at {fasta}"
        if build_bowtie:
            bowtie_location = build_bowtie_index(fasta=fasta, dest=dest)
            message += f" Bowtie index located at {bowtie_location}"


if __name__ == "__main__":
    merrycrispr()
