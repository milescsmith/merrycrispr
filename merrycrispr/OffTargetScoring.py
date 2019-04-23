import pandas as pd
import numpy as np
from regex import compile
from distutils.spawn import find_executable
from subprocess import check_call
from multiprocessing import cpu_count
from sys import stderr, exit


def off_target_discovery(spacer_df: pd.DataFrame,
                         cpus: int = 0,
                         refgenome: str = "",
                         reject: bool = False,
                         large_index_size: bool = False) -> str:

    spacer_df = spacer_df[spacer_df["spacer"].isin(spacer_df["spacer"].unique())]  # in case there are duplicates
    spacer_df['hash'] = spacer_df.apply(lambda x: hash(tuple(x)), axis=1)

    with open('temp.fa', 'w') as tempfile:
        for entry in spacer_df.iterrows():
            tempfile.writelines(f">{entry[1]['hash']}\n{entry[1]['spacer']}\n")

    program = find_executable("bowtie")
    if cpus is 0:
        cpus = cpu_count()

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


def off_target_scoring(bowtie_results_file: str,
                       nuclease_info: pd.Series,
                       spacer_df: pd.DataFrame,
                       count_threshold: int = 0) -> pd.DataFrame:

    mismatch_spots_re = compile('[0-9]{1,}')

    bowtie_results = pd.read_csv(bowtie_results_file,
                                 header=None,
                                 names=["hash",
                                        "strand",
                                        "refseq",
                                        "position",
                                        "seq",
                                        "readquality",
                                        "aligncount",
                                        "mismatches"],
                                 sep="\t").fillna(0)

    print(f"Total alignments from Bowtie: {bowtie_results.shape[0]}")

    spacer_df['offtarget_score'] = np.repeat(0, spacer_df.shape[0])
    spacer_df['number_matching'] = np.repeat(0, spacer_df.shape[0])

    # for each spacer
    for i in spacer_df["hash"].unique():
        # get everything for that ["hash"]
        matching_locations = bowtie_results[bowtie_results["hash"] == i]

        # if the number of mismatches is above a threshold, remove the spacer
        # if there are more than one perfect matches
        if matching_locations.shape[0] > count_threshold or \
                len(matching_locations[matching_locations["mismatches"] == 0].index) > 1:
            score = 0
        # if there is only one entry - no offtargets, assign a score of 0
        elif matching_locations.shape[0] == 1:
            score = 100
        # elif there are mismatch positions, get the positions
        # and make a list holding lists of those positions and score
        else:
            bounds = int(nuclease_info['start'][0]) - int(nuclease_info["end"][0])
            matching_locations = matching_locations.drop(
                matching_locations[matching_locations["position"].isin(bounds)].index).head()
            mmpos = [mismatch_spots_re.findall(str(_[1]['mismatches'])) for _ in matching_locations.iterrows()]
            score = sumofftargets(mmpos,
                                  start=nuclease_info["start"][0],
                                  end=nuclease_info["end"][0])
        spacer_df.loc[spacer_df["hash"] == i, "offtarget_score"] = score
        spacer_df.loc[spacer_df["hash"] == i, "number_matching"] = matching_locations.shape[0]

    return spacer_df


def scoreCas9offtarget(mismatched_positions: list,
                       start: int,
                       end: int) -> float:

    search_region = range(start, end+1)
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


def sumofftargets(offtargetlist: list,
                  start: int,
                  end: int) -> float:

    sum_score = sum(scoreCas9offtarget(x, start, end) for x in offtargetlist)
    if sum_score == 0:
        return 100
    else:
        return (1.0 / sum_score) * 100
