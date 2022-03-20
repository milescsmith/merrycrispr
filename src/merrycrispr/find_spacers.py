#!/usr/bin/env python3

import re
from functools import partial
from typing import List, Optional

import numpy as np
import pandas as pd
import pyfaidx
import regex
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from tqdm.auto import tqdm


def find_spacers(
    itemlist: pyfaidx.Fasta,
    nuclease_info: dict,
    restriction_sites: Optional[List[str]] = None,
    chunks: int = 8,
) -> pd.DataFrame:
    """Find protospacers in a sequence for a given nuclease.  The search region can be more
    expansive than just the PAM (for scoring purposes) and strand can be taken into account.
    `find_spacers()` will ignore sequences that have a poly(T) sequence, high GC content, or
    a motif matching given restriction nuclease sequence.

    Parameters
    ----------
    itemlist : :class:`~pyfaidx.Fasta`
        Parsed FASTA with sequences to examine for spacers
    nuclease_info : dict
        Information for the nuclease to use.  Required keys include `start`, `end`, `strand`,
        `pam`, `spacer_regex`
    restriction_sites : `List[str]`, optional (default: `None`)
        For a found spacer, ignore it if it contains the sequence for recognition by these
        restriction endonucleases.
    chunks : `int`, optional (default: 8)
        Number of pieces to divide the spacer dataframe into.  Higher number
        means less memory used at a time, but may result in slower processing

    Return
    ------
    :class:`~pandas.DataFrame`
    """
    spacer_regex = regex.compile(nuclease_info["spacer_regex"])
    spacer_start: int = nuclease_info["start"]
    spacer_end: int = nuclease_info["end"]

    # Set the restriction sites that we are going to make sure are not in our
    # spacers
    if restriction_sites:
        rsb = RestrictionBatch(restriction_sites)
    else:
        rsb = None

    # For each entry in the file (i.e. exonic sequence), find all of the
    # potential protospacer sequences.
    spacers_df = fasta_to_df(itemlist)

    tqdm.pandas(desc="finding forward spacers", unit="sequences")
    spacers_df["forward_spacers"] = spacers_df["sequence"].progress_apply(
        spacer_regex.findall
    )
    tqdm.pandas(desc="finding reverse spacers", unit="sequences")
    spacers_df["reverse_spacers"] = spacers_df["reverse_complement"].progress_apply(
        spacer_regex.findall
    )
    spacers_df = spacers_df.drop(columns=["sequence", "reverse_complement"])

    chunked_spacer_dfs = np.array_split(spacers_df, chunks)
    pivot_partial = partial(
        pivot_spacers,
        spacer_start=spacer_start,
        spacer_end=spacer_end,
        restriction_sites=rsb,
    )
    spacers_df = pd.concat(map(pivot_partial, chunked_spacer_dfs))

    # duplicates were sneaking in.
    spacers_df = spacers_df.groupby("spacer").first().reset_index()

    return spacers_df


def fasta_to_df(fasta: pyfaidx.Fasta) -> pd.DataFrame:
    """Convert the fasta file from seqextractor to a Pandas DataFrame

    Parameters
    ----------
    fasta : :class:`pyfaidx.Fasta`
        Parsed FASTA with sequences to examine for spacers that needs to be
        converted to a Pandas DataFrame

    Results
    ----------
    :class:`pd.DataFrame`
    """
    df = pd.DataFrame(
        [fasta[_].name.split("_") for _ in fasta.keys()],
        columns=["gene_name", "feature_id", "strand", "start", "stop", "seq_hash"],
    )

    df = df.astype(
        {
            "feature_id": "category",
            "gene_name": "category",
            "strand": "category",
            "start": np.uint32,
            "stop": np.uint32,
            "seq_hash": np.int32,
        },
        copy=False,
    )

    df["sequence"] = pd.Series([fasta[_][:].seq for _ in fasta.keys()])
    df["reverse_complement"] = pd.Series(
        [fasta[_][:].reverse.complement.seq for _ in fasta.keys()]
    )
    return df


def pivot_spacers(
    wide_df: pd.DataFrame,
    spacer_start: int,
    spacer_end: int,
    restriction_sites: Optional[RestrictionBatch] = None,
) -> pd.DataFrame:
    """`Regex.findall` returns a list of matches, with the entirety of the
    list being placed in a single cell of the `wide_df` for each sequence. This
    function makes a row for each spacer in those lists, filtering undesireable
    sequences (high/low GC content, presence of BsmBI sites, multiples), and
    scores them with the chose on-target algorithm.

    Parameters
    ----------
    wide_df : :class:`pandas.DataFrame`
    spacer_start : `int`
    spacer_end : `int`
    restriction_sites : :class:`Bio.Restriction.RestrictionBatch`, optional

    Returns
    -------
    :class:`pandas.DataFrame`
    """
    BsmBI_fwd = "GAGACG"
    BsmBI_rev = "CGTCTC"

    print("separating spacer lists")
    tqdm.pandas(desc="forward", unit="sequences")
    long_df = (
        wide_df.progress_apply(
            lambda x: pd.Series(x["forward_spacers"], dtype="str"), axis=1
        )
        .stack()
        .reset_index(level=1)
        .drop(columns=["level_1"])
        .rename(columns={0: "spacer"})
    )

    tqdm.pandas(desc="reverse", unit="sequences")
    long_df_rev = (
        wide_df.progress_apply(
            lambda x: pd.Series(x["reverse_spacers"], dtype="str"), axis=1
        )
        .stack()
        .reset_index(level=1)
        .drop(columns=["level_1"])
        .rename(columns={0: "spacer"})
    )

    print("appending lists")
    long_df = long_df.append(long_df_rev)

    print("eliminating duplicates")
    # eliminate anything that appears more than once
    # will probably have to run this again for the entire list, not just a chunk
    long_df_counts = (
        long_df["spacer"]
        .value_counts()
        .reset_index()
        .rename(columns={"index": "spacer", "spacer": "count"})
    )
    long_df = long_df.reset_index().merge(long_df_counts, on="spacer", left_index=True)
    long_df = long_df[long_df["count"] == 1].drop(columns="count").set_index("index")
    if restriction_sites:
        long_df = long_df[
            ~long_df["spacer"].apply(restriction_sites_present, rsb=restriction_sites)
        ]

    print("eliminating undesirable sequences")
    # eliminate those with a polyT or a BsmBI site within the spacer
    long_df = long_df[
        long_df["spacer"]
        .str[spacer_start:spacer_end]
        .str.match(f"^((?!T{{4,}}|{BsmBI_fwd}|{BsmBI_rev}).)*$")
    ]  # match antthing not followed by a polyT or a BsmBI recognition site in the forward seq or a BsmBI recognition site in the reverse seq

    print("eliminating high/low GC sequences")
    long_df = long_df[long_df["spacer"].apply(GC).between(20, 80)]

    print("merging with spacer info")
    wide_df = wide_df.drop(columns=["forward_spacers", "reverse_spacers"])
    long_df = wide_df.merge(long_df, how="inner", left_index=True, right_index=True)

    return long_df


def GC(seq: str) -> float:
    """Calculate the QC content of a string corresponding to a genetic sequence\f

    Parameters
    ----------
    seq : `str`
        DNA sequence to examine for GC content

    Returns
    -------
    `float`
    """
    gc = len(re.findall(string=seq, pattern="[GgCc]"))
    try:
        return gc * 100.0 / len(seq)
    except ZeroDivisionError:
        return 0.0


def restriction_sites_present(spacer: str, rsb: RestrictionBatch) -> List[int]:
    """Determine if and where a set of restriction sites are present in a
    sequence\f

    Parameters
    ----------
    spacer : `str`
        Spacer sequence to examine for restriction sites.

    Returns
    -------
    :class:`typing.List`[`int`]
    """

    sites = bool([_ for results in rsb.search(Seq(spacer)).values() for _ in results])
    return sites
