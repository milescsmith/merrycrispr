#!/usr/bin/env python3

import os
from typing import Optional

import click
import numpy as np
import pandas as pd
import pyfaidx
from pkg_resources import resource_filename

from merrycrispr.find_spacers import find_spacers
from merrycrispr.library_assembly import assemble_paired_library, assemble_library
from merrycrispr.off_target_scoring import off_target_discovery, off_target_scoring
from merrycrispr.on_target_scoring import on_target_scoring
from merrycrispr.seqextractor import (
    extract,
    extract_for_tss_adjacent,
    display_gtf_features,
    display_gtf_genes,
    display_gtf_geneids,
)
from merrycrispr.species_getter import (
    get_resources,
    build_bowtie_index,
    available_species,
)


@click.group()
def main():
    pass


@main.command()
@click.option(
    "--library_type",
    "-l",
    help="target library type. Accepted values are 'knockout', "
    "'repressor', 'activator', 'excision', or 'Cas13'",
    default=None,
    required=False,
    type=click.Choice(("knockout", "repressor", "activator", "excision", "Cas13")),
)
@click.option(
    "--gtf", "-g", help="input GTF/GFF file", default=None, required=False, type=str
)
@click.option("--fasta", "-f", help="FASTA sequence file", default=None, type=str)
@click.option("--output", "-o", help="output file", default=None, type=str)
@click.option("--feature_type", "-t", help="feature type", default=None)
@click.option("--gene_name", "-n", help="gene(s) to extract", default=None, type=str)
@click.option(
    "--bound",
    help="Retrieve a given number of bases on either side of the feature "
    "instead of the sequence corresponding to a feature",
    default=0,
    type=int,
)
@click.option(
    "--show_features",
    help="Scan an annotation file to identify the features present",
    is_flag=True,
)
@click.option(
    "--show_genes",
    help="Scan an annotation file to identify the genes present",
    default=False,
    is_flag=True,
)
@click.option(
    "--show_geneids",
    help="Scan an annotation file to identify the geneIDs present",
    default=False,
    is_flag=True,
)
@click.help_option()
def prep_sequences(
    gtf: str,
    library_type: str,
    fasta: str,
    output: str,
    feature_type: Optional[str] = None,
    gene_name: Optional[str] = None,
    bound: int = 0,
    show_features: Optional[bool] = False,
    show_genes: Optional[bool] = False,
    show_geneids: Optional[bool] = False,
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
    },
    skip_blank_lines=True,
)

AVAILABLE_NUCLEASES = ", ".join(NUCLEASES["nuclease"])


@main.command()
@click.option(
    "--input_sequences",
    "-i",
    help=f"Input FASTA file containing sequences to target.",
    default=None,
    required=False,
    type=str,
)
@click.option(
    "--output_library",
    "-p",
    help="Name of file to write library to (in CSV format).",
    default=None,
    required=False,
    type=str,
)
@click.option(
    "--nuclease",
    "-n",
    help=f"Cas nuclease to design for.  Current options include {AVAILABLE_NUCLEASES}",
    default="SpCas9",
    required=False,
    type=str,
)
@click.option(
    "--reference",
    "-r",
    help="Path to the directory containing the appropriate Bowtie index.",
    default=None,
    type=str,
)
@click.option(
    "--largeindex", help="", default=False, required=False, type=bool, is_flag=True
)
@click.option(
    "--rule_set",
    help="On-target score rule set to use. Current options include '1', '2', and 'Azimuth'",
    default="Azimuth",
    type=str,
)
@click.option(
    "--on_target_score_threshold",
    "-on",
    help="Spacers with an on-target score below this will be ignored.",
    default=0,
    required=False,
    type=int,
)
@click.option(
    "--off_target_score_threshold",
    "-off",
    help="Spacers with an off-target score below this will be ignored.",
    default=0,
    required=False,
    type=int,
)
@click.option(
    "--off_target_count_threshold",
    help="Spacers with more than this many off-targets will be ignored.",
    default=100,
    required=False,
    type=int,
)
@click.option(
    "--number_mismatches_to_consider",
    help="Number of mismatches to allow in potential off-targets.  A number between 0 and 3.  Without setting `off_target_count_threshold` appropriately, higher values passed to this option may result in extremely long run times.",
    default=3,
    type=int,
)
@click.option(
    "--spacers_per_feature",
    help="Number of spacers to find for each feature. Use '0' to return all spacers.",
    default=6,
    type=int,
)
@click.option(
    "--paired",
    help="Should spacers be designed to work as pairs (e.g. for excision)?",
    default=False,
    type=bool,
    is_flag=True,
)
@click.option(
    "--number_upstream_spacers",
    help=f"If designing paired spacers, number of spacers to design that target upstream of the "
    f"feature.",
    default=3,
    type=int,
)
@click.option(
    "--number_downstream_spacers",
    help=f"If designing paired spacers, number of spacers to "
    f"design that target downstream of the feature.",
    default=3,
    type=int,
)
# reenable this option just as soon as I figure out how to implement it.
# @click.option(
#     "--min_paired_distance",
#     help="If designing paired spacers, minimum space required between the up- and downstream "
#     "spacers.",
#     default=0,
#     type=int,
# )
@click.option(
    "--cores",
    "-c",
    help="Number of processors to use. By default, will use all available.",
    default=None,
    type=int,
)
@click.help_option()
def create_library(
    input_sequences: str = None,
    output_library: str = None,
    reference: str = None,
    restriction_sites: str = None,
    largeindex: bool = False,
    on_target_score_threshold: int = 0,
    off_target_score_threshold: int = 0,
    off_target_count_threshold: int = 100,
    number_mismatches_to_consider: int = 3,
    nuclease: str = "SpCas9",
    spacers_per_feature: int = 9,
    reject: bool = False,
    paired: bool = False,
    rule_set: str = "Azimuth",
    number_upstream_spacers: int = 0,
    number_downstream_spacers: int = 0,
    cores: int = 0,
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

    Return
    ------
    :type reference: object
    """
    targets = pyfaidx.Fasta(input_sequences)

    global NUCLEASES
    nuc = NUCLEASES[NUCLEASES["nuclease"] == nuclease].to_dict(orient="records")[0]

    spacers_df = find_spacers(
        itemlist=targets, nuclease_info=nuc, restriction_sites=restriction_sites
    )

    initialnumber = spacers_df.shape[0]

    spacers_df = on_target_scoring(
        ruleset=rule_set,
        spacers=spacers_df,
        on_target_score_threshold=on_target_score_threshold,
    )

    if spacers_df.shape[0] == 0:
        print("Sorry, no spacers matching that criteria were found")
        exit()
    else:
        print(
            f"Finished scoring spacers. {spacers_df.shape[0]} of {initialnumber} "
            f"spacers have an on-target score above the cutoff threshold of "
            f"{on_target_score_threshold}."
            f"\nBeginning Bowtie alignment..."
        )

    spacers_df["hash"] = spacers_df.apply(lambda x: hash(tuple(x)), axis=1)
    offtarget_results_file = off_target_discovery(
        spacers_df=spacers_df,
        nuclease_info=nuc,
        cpus=cores,
        refgenome=reference,
        large_index_size=largeindex,
        reject=reject,
        number_mismatches_to_consider=number_mismatches_to_consider,
    )

    spacers_df = off_target_scoring(
        otrf=offtarget_results_file,
        spacers_df=spacers_df,
        nuclease_info=nuc,
        off_target_score_threshold=off_target_score_threshold,
        off_target_count_threshold=off_target_count_threshold,
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


@main.command()
@click.option(
    "--species_value",
    "-s",
    help="Species value to retrieve sequences for.  Must be a value for the attribute chosen, "
    "i.e. 'canis_familiaris' if using the 'name' attribute",
    default=None,
    required=False,
    type=str,
)
@click.option(
    "--species_attribute",
    "-a",
    help="Species attribute to search when retrieving sequences. Acceptable values include "
    "'taxon_id', 'accession', 'aliases', 'division', 'groups', 'release', 'name', 'strain', "
    "'strain_collection', 'display_name', 'assembly', and 'common_name'",
    default=None,
    required=False,
    type=str,
)
@click.option(
    "--dest",
    "-d",
    help="Location in which to place downloaded files and, if selected, new "
    "Bowtie index.",
    default=None,
    required=False,
    type=str,
)
@click.option(
    "--build_bowtie",
    "-b",
    help="Build a bowtie index for the retrieved species",
    is_flag=True,
)
@click.option(
    "--show_available_builds",
    "-w",
    help="Display a list of the species available from Ensembl",
    is_flag=True,
)
@click.help_option()
def new_species(
    species_value: str,
    species_attribute: str,
    dest: Optional[str] = None,
    build_bowtie: bool = False,
    show_available_builds: bool = False,
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
            dest = os.path.expandvars("$PWD") + "/mc_resources"
        if not os.path.exists(dest):
            os.mkdir(dest)
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
    main()
