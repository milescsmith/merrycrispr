from distutils.spawn import find_executable
from multiprocessing import cpu_count
from subprocess import check_call
from sys import stderr, exit
from typing import Union

import numpy as np
import pandas as pd
from regex import compile
from serializable.primitive_types import NoneType


def scoreCas9offtarget(mismatched_positions: list,
                       start: int,
                       end: int) -> float:
    """ Calculate the likelihood a Cas9 protospacer will cut at a particular off-target site
    Equation from http://crispr.mit.edu/about

    The mismatch scoring algorithm from the Zhang group has three terms:
    1) the effect of a mismatch at a particular position
    2) the product of the mismatch scores, weighted by the mean distance between each mismatch
    3) a penalty for the number of mismatches
    Score is from 0 to 1, with higher scores indicating a higher likelihood the off-target will be cut.
    e.g. the score for mismatches at [15,16,17,18,19] is infinitesimally small, indicating that those
    mismatches are highly destablizing
    """
    search_region = range(start,end+1)
    mismatched_positions = [int(_)-1-4 for _ in mismatched_positions if int(_) in search_region]
    # remember that the for on target scoring, we have 4N-spacer-NGG-3N
    # for Cas9 off-target scoring, we only care about the spacer portion
    # also, Python uses 0-indexed arrays, so we also have to subtract 1
    # from the position reported by Bowtie
    if not mismatched_positions:
        score = 1
    else:
        # experimentally determined weighting penality for mismatch at each position
        M = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613,
             0.851, 0.731, 0.828, 0.615, 0.804, 0.685, 0.583]
        if len(mismatched_positions) == 1:  # if there is only one mismatch, we should ignore the second two terms.
            score = 1 - M[mismatched_positions[0]]
        else:
            nmm = len(mismatched_positions)
            mean_distance = (max(mismatched_positions) - min(mismatched_positions)) * 1.0 / (nmm - 1)
            term_2 = (1 / ((((19 - mean_distance) / 19) * 4) + 1))
            term_3 = 1.0 / (nmm ** 2)
            term_1 = 1
            for n in mismatched_positions:
                term_1 *= 1 - M[n - 20]
            score = term_1 * term_2 * term_3
    return score


def sumofftargets(offtargets: list,
                  start: int,
                  end: int) -> float:
    """
    Add all of the potential off-target scores together so that the higher the offtarget score, the
    more desirable the spacer
    :param offtargets:
    :param start:
    :param end:
    :return:
    """
    sum_score = sum(scoreCas9offtarget(x, start, end) for x in offtargets)
    if sum_score == 0:
        return 100
    else:
        return (1.0 / sum_score) * 100


def off_target_discovery(spacers_df: pd.DataFrame,
                         cpus: int = 0,
                         refgenome: str = "",
                         large_index_size: bool = False,
                         reject: Union[bool, int] = False) -> str:

    spacers_df = spacers_df[spacers_df["spacer"].isin(spacers_df["spacer"].unique())]
    spacers_df['hash'] = spacers_df.apply(lambda x: hash(tuple(x)), axis=1)
    if cpus is 0:
        cpus = cpu_count()
    program = find_executable("bowtie")

    with open('temp.fa', 'w') as tempfile:
        for entry in spacers_df.iterrows():
            tempfile.writelines(f">{entry[1]['hash']}\n{entry[1]['spacer']}\n")

    command = f"{program} -a -p {cpus}"
    if reject:
        command = command + f" -m {reject}"

    if large_index_size:
        command = command + f" --large-index {refgenome}"
    else:
        command = command + f" {refgenome}"

    command = command + " -f temp.fa offtargets.fa"

    try:
        check_call(command.split())
    except:
        stderr.write("Bowtie encountered an error. Please check the log file.")
        exit(-1)

    print("Bowtie finished.")
    bowtie_results_file = 'offtargets.fa'

    return bowtie_results_file


def off_target_scoring(otrf: str,
                       spacers_df: pd.DataFrame,
                       nuclease_info: pd.Series,
                       offtarget_score_threshold: int,
                       offtarget_count_threshold: Union[int, NoneType]) -> object:
    """

    :param offtarget_score_threshold:
    :param otrf: path to the results from Bowtie
    :param spacers_df: Pandas dataframe containing spacers.  Format should be
        {'gene_name','feature_id','start','stop','strand','spacer'}
    :param nuclease_info: Pandas series with nuclease characteristics from nuclease_list.csv
    :param offtarget_count_threshold: Integer value of the number of potential mismatches that should be tolerated.
    :return: A Pandas dataframe matching the one passed to spacers_df containing off-target scores
    """

    mmpos_re = compile('[0-9]{1,}')

    bowtie_results = pd.read_csv(otrf,
                                 header=None,
                                 names=["hash",
                                        "strand",
                                        "refseq",
                                        "position",
                                        "seq",
                                        "readquality",
                                        "aligncount",
                                        "mismatches"],
                                 sep="\t")

    print(f"Total alignments from Bowtie: {bowtie_results.shape[0]}")
    bowtie_results = bowtie_results.fillna(0)

    spacers_df['offtarget_score'] = np.repeat(0, spacers_df.shape[0])
    spacers_df['number_matching'] = np.repeat(0, spacers_df.shape[0])

    # for each spacer
    for i in spacer_df["hash"].unique():
        # get everything for that ["hash"]
        matching_locations = bowtie_results[bowtie_results["hash"] == i]

        # if the number of mismatches is above a threshold, remove the spacer
        # if there are more than one perfect matches
        if offtarget_count_threshold & matching_locations.shape[0] > offtarget_count_threshold or \
                len(matching_locations[matching_locations["mismatches"] == 0].index) > 1:
            score = 0
        # if there is only one entry - no offtargets, assign a score of 0
        elif matching_locations.shape[0] == 1:
            score = 100
        # elif there are mismatch positions, get the positions, make a list holding lists of those positions, and score
        else:
            bounds = range(int(spacers_df.loc[spacers_df["hash"] == i, "start"]), \
                           int(spacers_df.loc[spacers_df["hash"] == i, "stop"]))  # ideally, this would take refseq into consideration
            matching_locations = matching_locations.drop(
                matching_locations[matching_locations["position"].isin(bounds)].index)
            mmpos = [mmpos_re.findall(str(_[1]['mismatches'])) for _ in matching_locations.iterrows()]
            score = sumofftargets(mmpos,
                                  start=nuclease_info["start"][0],
                                  end=nuclease_info["end"][0])
        spacers_df.loc[spacers_df["hash"] == i, "offtarget_score"] = score
        spacers_df.loc[spacers_df["hash"] == i, "number_matching"] = matching_locations.shape[0]

    if offtarget_score_threshold
    return spacer_df
