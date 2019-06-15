#!/usr/bin/env python3

from typing import List, Optional

import click
import pandas as pd
import pyfaidx
from gtfparse import read_gtf


@click.command()
@click.option(
    "--library_type",
    "-l",
    help="target library type. Accepted values are 'knockout', "
    "'repressor/activator', 'excision', or 'Cas13'",
    default=None,
    required=False,
    type=str,
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
    default=None,
    type=int,
)
@click.option(
    "--show_features",
    help="Scan a GFF file to identify the features present",
    is_flag=True,
)
@click.option(
    "--show_genes",
    help="Scan a GFF file to identify the genes present",
    default=False,
    is_flag=True,
)
@click.option(
    "--show_geneids",
    help="Scan a GFF file to identify the geneIDs present",
    default=False,
    is_flag=True,
)
@click.help_option()
def main(
    gtf: str,
    library_type: str,
    fasta: str,
    output: str,
    feature_type: str,
    gene_name: str,
    bound: int,
    show_features: bool,
    show_genes: bool,
    show_geneids: bool,
) -> None:
    """A utility for extracting sequences from a FASTA file for a given GFF annotation

    Parameters
    ----------
    gtf : `str`
        GTF/GFF file that matches the input FASTA file. Preferably one from Ensembl/GENCODE. Gzipped GTF/GFF
        files are acceptable, though their use may impose a performance penality.
    library_type : `str`, optional
        Target library type.  If provided, other parameters will be set to generate an appropriate
        set of features.
    fasta : `str`
        FASTA file corresponding to the input GTF/GFF.  Preferably one from Ensembl/GENCODE. Gzipped
        FASTA files will work, but their use results in a very substantial performance loss.
    output : `str`, optional (default: `None`)
        File to write extracted sequences to.  If not provided, seqextractor will instead return a `list` of
        :class:`~pyfaidx.Sequence` objects.
    feature_type : (default: `None`)
        A list of target feature types (i.e. gene, exon, intron, etc...).  Multiple types are allowed.
    gene_name : `str`, optional (default: `None`)
        If provided, seqextractor will only extract sequences for that list of genes.
    bound : `int`, optional (default: None)
        Instead of extracting the sequences corresponding to the feature itself, extract an area equal
        to the length of `bound` on either side of the feature.
    show_features : `bool` (default: `False`)
        Instead of extracting anything, display a list of the feature types present in the GTF/GFF.
    show_genes : `bool` (default: `False`)
        Instead of extracting anything, display a list of the genes present in the GTF/GFF.
    show_geneids : `bool` (default: `False`)
        Instead of extracting anything, display a list of the geneids present in the GTF/GFF.

    Returns
    -------
    Nothing

    """

    if gene_name:
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
        elif library_type == "repressor/activator":
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
            gene_name=gene_name.split(),
            boundary=int(bound),
        )
    return None


def extract(
    gtffile: str,
    fastafile: str,
    feature_type: str,
    outfile: str,
    gene_name: List[str],
    boundary: int = 0,
    **kwargs,
) -> Optional[List[pyfaidx.Sequence]]:
    """Match the features present in a GTF/GFF and use it to divide the sequences of a FASTA file into portions
    corresponding to those features.

    Parameters
    ----------
    gtffile : `str`
        GTF/GFF file that matches the input FASTA file.  Preferably one from Ensembl/GENCODE.  Gzipped GTF/GFF
        files are acceptable, though their use may impose a performance penality.
    fastafile :`str`
        FASTA file corresponding to the input GTF/GFF.  Preferably one from Ensembl/GENCODE.    Gzipped
        FASTA files will work, but their use results in a very substantial performance loss.
    feature_type : (default: `None`)
        A list of target feature types (i.e. gene, exon, intron, etc...).  Multiple types are allowed.
    outfile : `str`, optional (default: `None`)
        File to write extracted sequences to.
    gene_name : `str`, optional (default: `None`)
        If provided, seqextractor will only extract sequences for that list of genes.
    boundary : `int`, optional (default: None)
        Instead of extracting the sequences corresponding to the feature itself, extract an area equal
        to the length of boundary on either side of the feature.
    kwargs :
        Additional optional arguments. Specifically, "exon_number" is used to make a guess at the location
        of transcriptional start sites.

    Returns
    -------
    Optional.  If no output file is specified, then a list of :class:`~pyfaidx.Sequence` objects.
    """

    print("Parsing GTF/GFF file.")
    records = read_gtf(filepath_or_buffer=gtffile)
    if gene_name:
        records = records[records["gene_name"].isin(gene_name)]
    if feature_type:
        records = records[records["feature"] == feature_type]
        records = records[
            [
                "seqname",
                "feature",
                "start",
                "end",
                "strand",
                "frame",
                "gene_name",
                f"{feature_type}_id",
            ]
        ].drop_duplicates()
    else:
        records = records[
            ["seqname", "feature", "start", "end", "strand", "frame", "gene_name"]
        ].drop_duplicates()

    if "exon_number" in kwargs.keys():
        records = records[records["exon_number"] == kwargs["exon_number"]]

    print(f"{len(records)} total records found.")

    print(
        f"Loading the sequences in {fastafile}."
        f"Note: if this is the first time opening this file, "
        "it may take a few moments as an index is built."
    )
    sequences = pyfaidx.Fasta(fastafile)
    print(f"Finished loading {fastafile}")

    records["seq_hash"] = records.apply(lambda x: hash(tuple(x)), axis=1)
    if boundary > 0:
        records = split_records(records, boundary)
    final_list = [
        match_seq(_, sequences)
        for _ in [pd.Series(records.loc[_, :]) for _ in records.index]
    ]

    if outfile:
        with open(outfile, "w") as o_file:
            for entry in final_list:
                o_file.writelines(f"> {entry.fancy_name}\n{entry.seq}\n")
    else:
        return final_list


def display_gtf_features(gtffile: str) -> None:
    """Display the features present in a GTF/GFF

    Parameters
    ----------
    gtffile : `str`
        GTF/GFF file that matches the input FASTA file.  Preferably one from Ensembl/GENCODE.  Gzipped GTF/GFF
        files are acceptable, though their use may impose a performance penality.

    Returns
    -------
    `None`
    """

    gtf = read_gtf(gtffile)
    feature_set = set(gtf["feature"])

    print(f"{len(feature_set)} features found.  These include:")
    [print(_) for _ in feature_set]


def display_gtf_genes(gtffile: str, feature_type: list = None) -> None:
    """Display the genes present in a GTF/GFF

    Parameters
    ----------
    gtffile : `str`
        GTF/GFF file that matches the input FASTA file.  Preferably one from Ensembl/GENCODE.  Gzipped GTF/GFF
        files are acceptable, though their use may impose a performance penality.
    feature_type :

    Returns
    -------
    `None`
    """

    gtf = read_gtf(gtffile)

    if feature_type is not None:
        gtf = gtf[gtf["feature"] == feature_type]

    gene_set = set(gtf["gene_name"])

    print(f"{len(gene_set)} genes found.  These include:")
    [print(_) for _ in gene_set]
    return None


def display_gtf_geneids(gtffile: str, feature_type: list = None) -> None:
    """Display the geneids present in a GTF/GFF

    Parameters
    ----------
    gtffile : `str`
        GTF/GFF file that matches the input FASTA file.  Preferably one from Ensembl/GENCODE.  Gzipped GTF/GFF
        files are acceptable, though their use may impose a performance penality.
    feature_type :

    Returns
    -------
    `None`
    """

    gtf = read_gtf(gtffile)

    if feature_type is not None:
        gtf = gtf[gtf.feature == feature_type]

    gene_set = set(gtf["gene_id"])

    print(f"{len(gene_set)} genes found.  These include:")
    [print(_) for _ in gene_set]
    return None


def extract_for_tss_adjacent(
    gtffile: str, fastafile: str, outfile: str, gene_name: str, boundary: int = 100
) -> Optional[List[pyfaidx.Sequence]]:
    """Match the features present in a GTF/GFF and use it to divide the sequences of a FASTA file into portions
    corresponding to those features.  For a given feature, `extract_for_tss_adjacent` only extracts an area 5' to
    the start of exon 1

    Parameters
    ----------
    gtffile : `str`
        GTF/GFF file that matches the input FASTA file.  Preferably one from Ensembl/GENCODE.  Gzipped GTF/GFF
        files are acceptable, though their use may impose a performance penality.
    fastafile :`str`
        FASTA file corresponding to the input GTF/GFF.  Preferably one from Ensembl/GENCODE.    Gzipped
        FASTA files will work, but their use results in a very substantial performance loss.
    outfile : `str`, optional (default: `None`)
        File to write extracted sequences to.
    gene_name : `str`, optional (default: `None`)
        If provided, seqextractor will only extract sequences for that list of genes.
    boundary : `int`, optional (default: None)
        Instead of extracting the sequences corresponding to the feature itself, extract an area equal
        to the length of boundary on either side of the feature.

    Returns
    -------
    `None` if an output file is spcified, else then a list of Optional[List[:class:`~pyfaidx.Sequence`]] objects.
    """

    # read and parse in GTF
    print("Parsing GTF/GFF file.")
    records = read_gtf(gtffile)
    records = records[records["feature"] == "exon"]
    if gene_name:
        records = records[records["gene_name"].isin(gene_name)]

    # break up the genome into forward and reverse strands
    # for forward strand genes, we want the lowest coordinate for exon 1 of each gene
    # for reverse strand genes, we want the highest
    forward_starts = (
        records[records["strand"] == "+"]
        .groupby("gene_name")
        .apply(lambda x: x[pd.to_numeric(x["exon_number"]) == 1])
        .reset_index(drop=True)
        .groupby("gene_name")
        .apply(lambda y: y.nsmallest(1, "start"))
        .reset_index(drop=True)
    )

    reverse_starts = (
        records[records["strand"] == "-"]
        .groupby("gene_name")
        .apply(lambda x: x[pd.to_numeric(x["exon_number"]) == 1])
        .reset_index(drop=True)
        .groupby("gene_name")
        .apply(lambda y: y.nlargest(1, "start"))
        .reset_index(drop=True)
    )

    predicted_tss = pd.concat([forward_starts, reverse_starts])
    predicted_tss = predicted_tss[
        ["seqname", "feature", "start", "end", "strand", "frame", "gene_name"]
    ].drop_duplicates()

    print(f"{len(predicted_tss)} total records found.")

    print(
        f"Loading the sequences in {fastafile}."
        f"Note: if this is the first time opening this file, "
        "it may take a few moments as an index is built."
    )
    sequences = pyfaidx.Fasta(fastafile)
    print(f"Finished loading {fastafile}")

    # for our list of predicted start sites, extract +/- an interval surrounding the TSS
    # what we want is:
    #  -boundary --- start --- +boundary
    # easiest way to do this is to make our end coordinate the start coordinate plus the boundary
    # value and make the new start the original start - boundary
    predicted_tss["gene_name"] += "_TSS"
    predicted_tss["end"] = predicted_tss["start"] + boundary
    predicted_tss["start"] = predicted_tss["start"] - boundary
    final_list = [
        match_seq(_, sequences)
        for _ in [pd.Series(predicted_tss.loc[_, :]) for _ in predicted_tss.index]
    ]

    if outfile:
        with open(outfile, "w") as o_file:
            for entry in final_list:
                o_file.writelines(f"> {entry.fancy_name}\n{entry.seq}\n")
    else:
        return final_list


def match_seq(rec: pd.Series, sequences: pyfaidx.Fasta) -> pyfaidx.Sequence:
    """Given a feature in a GTF/GFF read in by gtfparse, match_seq() will extract the corresponding DNA sequence
    and create a new pyfaidx.Sequence object

    Parameters
    ----------
    rec: :class:`~pandas.Series`
        Information for a feature (i.e. gene, exon, etc...). Requires the following indices: strand, gene_name,
        feature, strand, start, end, seq_hash
    sequences: :class:`~pyfaidx.Sequence`
        Object containing sequences to match against the positions in the index.

    Returns
    -------
    :class:`~pyfaidx.Sequence object` with annotation from `rec` and sequence information from `sequences`.
    """

    try:
        rev = True if rec["strand"] == "-" else False

        seq = pyfaidx.Sequence(
            name=f"{rec['gene_name']}_"
            f"{rec['feature']}_"
            f"{rec['strand']}_"
            f"{rec['start']}_"
            f"{rec['end']}_"
            f"{rec['seq_hash']}",
            seq=sequences.get_seq(
                name=rec["seqname"], start=rec["start"], end=rec["end"], rc=rev
            ).seq,
        )
        return seq
    except ValueError:
        print(
            f"problem with {rec['gene_name']} {rec['start']} "
            f"{rec['end']} {rec['seqname']} {rec['strand']}"
        )


def split_records(rec: pd.DataFrame, padding: int) -> pd.DataFrame:
    """    Mutates a DataFrame of features from a GTF/GFF file into features corresponding to the upstream
    and downstream regions of those features.

    Parameters
    ----------
    rec: :class:`~Pandas.DataFrame`
        A dataframe with all features of interest.
    padding: `int`
        The amount of space on either size of the feature to return.

    Returns
    -------
    A new :class:`~Pandas.DataFrame` with twice as many rows as the input, with
        each row corresponding to an upstream and downstream window of sequence
    """

    rec_pos_upstream = rec[rec["strand"] == "+"].reset_index()
    rec_pos_downstream = rec[rec["strand"] == "+"].reset_index()
    rec_neg_upstream = rec[rec["strand"] == "-"].reset_index()
    rec_neg_downstream = rec[rec["strand"] == "-"].reset_index()

    rec_pos_upstream["end"] = rec_pos_upstream["start"]
    rec_pos_upstream["start"] -= padding
    rec_pos_upstream["gene_name"] += "-upstream"
    # we cannot use an underscore because that is used later in the
    # pyfaidx.Sequence objects' fancy name
    rec_pos_downstream["start"] = rec_pos_downstream["end"]
    rec_pos_downstream["end"] += padding
    rec_pos_downstream["gene_name"] += "-downstream"

    rec_neg_downstream["end"] = rec_neg_downstream["start"]
    rec_neg_downstream["start"] -= padding
    rec_neg_downstream["gene_name"] += "-downstream"
    # we cannot use an underscore because that is used later in the
    # pyfaidx.Sequence objects' fancy name
    rec_neg_upstream["start"] = rec_neg_upstream["end"]
    rec_neg_upstream["end"] += padding
    rec_neg_upstream["gene_name"] += "-upstream"

    return pd.concat(
        [rec_pos_upstream, rec_pos_downstream, rec_neg_upstream, rec_neg_downstream]
    ).reset_index()


if __name__ == "__main__":
    main()
