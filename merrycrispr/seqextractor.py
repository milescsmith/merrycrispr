#!/usr/bin/env python3

from typing import List, Optional

import pandas as pd
import pyfaidx
from gtfparse import read_gtf
from tqdm import tqdm


def extract(
    gtffile: str,
    fastafile: str,
    feature_type: str,
    outfile: str,
    gene_name: List[str],
    boundary: int = 0,
    **kwargs,
) -> Optional[List[pyfaidx.Sequence]]:
    """Match the features present in a GTF/GFF and use it to divide the sequences of a FASTA file
    into portions corresponding to those features.

    Parameters
    ----------
    gtffile : `str`
        GTF/GFF file that matches the input FASTA file.  Preferably one from Ensembl/GENCODE.
        Gzipped GTF/GFF files are acceptable, though their use may impose a performance penality.
    fastafile :`str`
        FASTA file corresponding to the input GTF/GFF.  Preferably one from Ensembl/GENCODE. Gzipped
        FASTA files will work, but their use results in a very substantial performance loss.
    feature_type : (default: `None`)
        A list of target feature types (i.e. gene, exon, intron, etc...).  Multiple types are
        allowed.
    outfile : `str`, optional (default: `None`)
        File to write extracted sequences to.
    gene_name : `str`, optional (default: `None`)
        If provided, seqextractor will only extract sequences for that list of genes.
    boundary : `int`, optional (default: None)
        Instead of extracting the sequences corresponding to the feature itself, extract an area
        equal to the length of boundary on either side of the feature.
    kwargs :
        Additional optional arguments. Specifically, "exon_number" is used to make a guess at the
        location of transcriptional start sites.

    Returns
    -------
    Optional.  If no output file is specified, then a list of :class:`~pyfaidx.Sequence` objects.
    """

    print("Parsing GTF/GFF file.")
    records = read_gtf(filepath_or_buffer=gtffile, chunksize=1024)
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

    tqdm.pandas(desc="Adding hash")
    records["seq_hash"] = records.progress_apply(lambda x: hash(tuple(x)), axis=1)
    if boundary > 0:
        records = split_records(records, boundary)
    final_list = [
        match_seq(_, sequences)
        for _ in tqdm(
            [pd.Series(records.loc[_, :]) for _ in records.index],
            total=len(records),
            unit="items",
        )
    ]

    if outfile:
        with open(outfile, "w") as o_file:
            for entry in final_list:
                o_file.writelines(f">{entry.fancy_name}\n{entry.seq}\n")
    else:
        return final_list


def extract_for_tss_adjacent(
    gtffile: str, fastafile: str, outfile: str, gene_name: str, boundary: int = 100
) -> Optional[List[pyfaidx.Sequence]]:
    """Match the features present in a GTF/GFF and use it to divide the sequences of a FASTA file
    into portions corresponding to those features.  For a given feature, `extract_for_tss_adjacent`
    only extracts an area 5' to the start of exon 1

    Parameters
    ----------
    gtffile : `str`
        GTF/GFF file that matches the input FASTA file.  Preferably one from Ensembl/GENCODE.
        Gzipped GTF/GFF files are acceptable, though their use may impose a performance penality.
    fastafile :`str`
        FASTA file corresponding to the input GTF/GFF.  Preferably one from Ensembl/GENCODE. Gzipped
        FASTA files will work, but their use results in a very substantial performance loss.
    outfile : `str`, optional (default: `None`)
        File to write extracted sequences to.
    gene_name : `str`, optional (default: `None`)
        If provided, seqextractor will only extract sequences for that list of genes.
    boundary : `int`, optional (default: None)
        Instead of extracting the sequences corresponding to the feature itself, extract an area
        equal to the length of boundary on either side of the feature.

    Returns
    -------
    `None` if an output file is spcified, else then a list of
    Optional[List[:class:`~pyfaidx.Sequence`]] objects.
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
        for _ in tqdm(
            [pd.Series(predicted_tss.loc[_, :]) for _ in predicted_tss.index],
            total=len(predicted_tss),
        )
    ]

    if outfile:
        with open(outfile, "w") as o_file:
            for entry in final_list:
                o_file.writelines(f">{entry.fancy_name}\n{entry.seq}\n")
    else:
        return final_list


def match_seq(rec: pd.Series, sequences: pyfaidx.Fasta) -> pyfaidx.Sequence:
    """Given a feature in a GTF/GFF read in by gtfparse, match_seq() will extract the corresponding
    DNA sequence and create a new pyfaidx.Sequence object

    Parameters
    ----------
    rec : :class:`~pandas.Series`
        Information for a feature (i.e. gene, exon, etc...). Requires the following indices: strand,
        gene_name, feature, strand, start, end, seq_hash
    sequences : :class:`~pyfaidx.Sequence`
        Object containing sequences to match against the positions in the index.

    Returns
    -------
    :class:`~pyfaidx.Sequence object` with annotation from `rec` and sequence information from
    `sequences`.
    """

    try:
        rev: bool = bool(rec["strand"] == "-")

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
    """Mutates a DataFrame of features from a GTF/GFF file into features corresponding to the
    upstream and downstream regions of those features.

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


def display_gtf_features(gtffile: str) -> None:
    """Display the features present in a GTF/GFF

    Parameters
    ----------
    gtffile : `str`
        GTF/GFF file that matches the input FASTA file.  Preferably one from Ensembl/GENCODE.
        Gzipped GTF/GFF files are acceptable, though their use may impose a performance penality.

    Returns
    -------
    `None`
    """

    gtf = read_gtf(gtffile)
    feature_set = set(gtf["feature"])

    print(f"{len(feature_set)} features found.  These include:")
    for _ in feature_set:
        print(_)


def display_gtf_genes(gtffile: str, feature_type: Optional[List[str]] = None):
    """Display the genes present in a GTF/GFF

    Parameters
    ----------
    gtffile : `str`
        GTF/GFF file that matches the input FASTA file.  Preferably one from Ensembl/GENCODE.
        Gzipped GTF/GFF files are acceptable, though their use may impose a performance penality.
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
    for _ in gene_set:
        print(_)


def display_gtf_geneids(gtffile: str, feature_type: Optional[List[str]] = None):
    """Display the geneids present in a GTF/GFF

    Parameters
    ----------
    gtffile : `str`
        GTF/GFF file that matches the input FASTA file.  Preferably one from Ensembl/GENCODE.
        Gzipped GTF/GFF files are acceptable, though their use may impose a performance penality.
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
    for _ in gene_set:
        print(_)
