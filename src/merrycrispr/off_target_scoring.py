#!/usr/bin/env python3

import sys
from distutils.spawn import find_executable
from multiprocessing import cpu_count
from subprocess import check_call
from tempfile import mkstemp
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import regex
from tqdm.auto import tqdm


def hsu_offtarget_score(mismatched_positions: List[int]) -> float:
    """Calculate the likelihood a Cas9 protospacer will cut at a particular off-target site
    Equation from http://crispr.mit.edu/about

    The mismatch scoring algorithm from the Zhang group has three terms:
    1) the effect of a mismatch at a particular position
    2) the product of the mismatch scores, weighted by the mean distance between each mismatch
    3) a penalty for the number of mismatches
    Score is from 0 to 1, with higher scores indicating a higher likelihood the off-target will be
    cut. e.g. the score for mismatches at [15,16,17,18,19] is infinitesimally small, indicating
    that those mismatches are highly destablizing
    \f
    Parameters
    ----------
    mismatched_positions : :class:`~typing.List`[`int`]

    Return
    ------
    `float`
    """

    mismatched_positions = [int(_) for _ in mismatched_positions]

    # remember that the for on target scoring, we have 4N-spacer-NGG-3N
    # for Cas9 off-target scoring, we only care about the spacer portion
    # also, Python uses 0-indexed arrays, so we also have to subtract 1
    # from the position reported by Bowtie
    if not mismatched_positions:
        score = 1
    else:
        # experimentally determined weighting penality for mismatch at each
        # position
        M = [
            0,
            0,
            0.014,
            0,
            0,
            0.395,
            0.317,
            0,
            0.389,
            0.079,
            0.445,
            0.508,
            0.613,
            0.851,
            0.731,
            0.828,
            0.615,
            0.804,
            0.685,
            0.583,
        ]
        # if there is only one mismatch, we should ignore the second two terms.
        if len(mismatched_positions) == 1:
            try:
                score = 1.0 - M[mismatched_positions[0]]
            except Exception as error:
                print(error)
        else:
            nmm = len(mismatched_positions)
            mean_distance = (
                (max(mismatched_positions) - min(mismatched_positions))
                * 1.0
                / (nmm - 1)
            )
            term_2 = 1 / ((((19 - mean_distance) / 19) * 4) + 1)
            term_3 = 1.0 / (nmm**2)
            term_1 = 1
            for n in mismatched_positions:
                try:
                    term_1 *= 1 - M[n - 20]
                except Exception as error:
                    print(error)
                    print(f"M: {M}, n: {n}\n")
                    sys.exit(1)
            score = np.prod([term_1, term_2, term_3])
    return score


def sumofftargets(offtargets: List[List[int]], rule_set: Optional[str] = None) -> float:
    """Add all of the potential off-target scores together so that the higher
    the offtarget score, the more desirable the spacer
    \f
    Parameters
    ----------
    offtargets : :class:`~typing.List`[:class:`~typing.List`[`int`]]
    start : `int`
    end : `int`

    Return
    ------
    `float`
    """
    if rule_set is None:
        return 0
    elif isinstance(rule_set, str):
        if rule_set.lower() == "hsu":
            sum_score = np.sum([hsu_offtarget_score(x) for x in offtargets])
        elif rule_set.lower() == "none":
            sum_score = 0

        if sum_score == 0:
            return 100
        else:
            return (1.0 / sum_score) * 100


def off_target_discovery(
    spacers_df: pd.DataFrame,
    nuclease_info: dict,
    cpus: int = 0,
    refgenome: Optional[str] = None,
    large_index_size: bool = False,
    reject: Optional[int] = None,
    number_mismatches_to_consider: int = 3,
    verbose: bool = False,
) -> str:
    """Identify potential protospacer off-targets using Bowtie.
    \f
    Parameters
    ----------
    spacers_df : :class:`~pandas.DataFrame`
    cpus : `int`
    refgenome : :class:`~typing.Optional`[`str`]
    large_index_size : `bool`
    reject : :class:`~typing.Optional`[`int`]
        Passed to the `-m` Bowtie parameter.  Suppress all alignments for a
        particular spacer if more than this number of reportable alignments
        exist for it.
    number_mismatches_to_consider : `int`
        An integer between 0 and 3.  Passed to the `-v` Bowtie parameter.
        Report alignments with at most this number of mismatches.

    Return
    ------
    `str`
    """
    if refgenome is None:
        raise ValueError("No reference Bowtie index provided")
    # keep only the two columns necessary right now
    spacers_df = spacers_df.loc[:, ["hash", "spacer"]].drop_duplicates()
    if cpus == 0:
        cpus = cpu_count()
    program = find_executable("bowtie")

    _, spacers_to_score = mkstemp()
    _, off_target_results = mkstemp()
    with open(spacers_to_score, "w") as tempfile:
        for entry in spacers_df.iterrows():
            tempfile.writelines(
                f">{entry[1]['hash']}\n{entry[1]['spacer'][nuclease_info['start']-1:nuclease_info['end']-1]}\n"
            )
    print(f"targets are in {spacers_to_score}")
    command = f"{program} -a -p {cpus} -v {number_mismatches_to_consider}"
    if reject:
        command += f" -m {reject}"

    if large_index_size:
        command += f" --large-index {refgenome}"
    else:
        command += f" {refgenome}"

    command += f" -f {spacers_to_score} {off_target_results}"

    try:
        if verbose:
            print(f"bowtie command: {command.split()}")
        check_call(command.split())
    except BaseException:
        raise SystemExit("Bowtie encountered an error. Please check the log file.")

    if verbose:
        print(f"Bowtie finished. Results at {off_target_results}")

    return off_target_results


def off_target_scoring(
    otrf: str,
    spacers_df: pd.DataFrame,
    nuclease_info: Dict[str, Any],
    rule_set: Optional[str] = None,
    off_target_score_threshold: int = 0,
    off_target_count_threshold: Optional[int] = 100,
    verbose: bool = False,
) -> object:
    """Calculate a cumulative off-target score for a protospacer
    \f
    Parameters
    -----------
    otrf : `str`
        Path to the results from Bowtie
    spacers_df : :class:`~pandas.DataFrame`
        Dataframe containing spacers.  Format should be `{'gene_name',
        'feature_id', 'start','stop','strand','spacer'}`
    nuclease_info : `str`
        dictionary series with nuclease characteristics from nuclease_list.csv
    off_target_score_threshold : `int`
        Total off-target score threshold beyond which a spacer is rejected.
        Ranges from 0 to 100.
    off_target_count_threshold : `int`, default: 100
        Number of potential mismatches that should be tolerated.  Spacers
        exceeding the threshold will be discarded
    verbose : `bool`

    Return
    -------
    :class:`~pandas.DataFrame` matching the one passed to spacers_df containing
    off-target scores
    """

    bowtie_results = pd.read_csv(
        otrf,
        header=None,
        names=[
            "hash",
            "strand",
            "refseq",
            "position",
            "seq",
            "readquality",
            "aligncount",
            "mismatches",
        ],
        usecols=["hash", "mismatches"],
        dtype={"hash": "int64", "mismatches": "str"},
        na_filter=False,
        skip_blank_lines=True,
        sep="\t",
        memory_map=True,
    )

    if verbose:
        print(f"Total alignments from Bowtie: {bowtie_results.shape[0]}")

    # We need to reduce the number of spacers we examine.  For the most part,
    # those with a lot of potential off-targets (>1000?) have really low
    # scores and are worthless.  Some have >10,000 (!) potential off-targets
    # and should just be thrown out.
    results_count = bowtie_results.groupby("hash").agg("count").reset_index()
    filtered_results = bowtie_results[
        bowtie_results["hash"].isin(
            results_count[results_count["mismatches"] < off_target_count_threshold][
                "hash"
            ]
        )
    ]

    # Keep only those spacers that have fewer than our cutoff
    spacers_df = spacers_df[spacers_df["hash"].isin(filtered_results["hash"])]

    mmpos = regex.compile("[0-9]{1,}")
    tqdm.pandas(desc="converting mismatches", unit="spacers")
    filtered_results["locations"] = filtered_results["mismatches"].progress_apply(
        mmpos.findall
    )

    tqdm.pandas(desc="collapsing mismatches", unit="spacers")
    collapsed_results = (
        filtered_results.groupby("hash")
        .progress_apply(lambda x: x["locations"].values)
        .reset_index()
        .rename(index=str, columns={0: "locations"})
    )

    tqdm.pandas(desc="scoring mismatches", unit="spacers")
    collapsed_results["off_target_score"] = collapsed_results.apply(
        lambda x: sumofftargets(x["locations"], rule_set=rule_set), axis=1
    )
    spacers_df = spacers_df.merge(collapsed_results, on="hash")

    tqdm.pandas("counting off-targets", unit="spacers")
    spacers_df["off_targets"] = spacers_df.progress_apply(
        lambda x: len(x["locations"]) - 1, axis=1
    )
    spacers_df = spacers_df.drop(columns=["locations"])

    spacers_df = spacers_df[spacers_df["off_target_score"] > off_target_score_threshold]
    return spacers_df
