# TODO: add ability to limit guides to the template strand
# TODO: add ability to cut with Cas13a, including limiting guides to those that would bind
# mRNA and exonic only regions, maybe even exon-exon boundaries

__author__ = "milescsmith"
__email__ = "mileschristiansmith@gmail.com"

from functools import partial
from multiprocessing import Manager, Pool
from typing import Union

import click
import numpy as np
import pandas as pd
import progressbar
import pyfaidx
from Bio.Alphabet import IUPAC
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from azimuth import model_comparison
from regex import compile

from on_target_score_calculator import calc_score
from .libraryassembly import assemble_paired_library, assemble_library
from .offtargetscoring import off_target_discovery, off_target_scoring


# from subprocess import check_output


@click.command()
@click.option("--input", "-i",
              help="Input FASTA file containing sequences to target.",
              default=None,
              required=True,
              type=str)
@click.option("--output", "-p",
              help="Name of file to write library to (in CSV format).",
              default=None,
              required=True,
              type=str)
@click.option("--nuclease", "-n",
              help="Cas nuclease to design for.  Current options include 'SpCas9', 'SaCas9', 'Cpf1', and 'Cas13'",
              default="SpCas9",
              required=False,
              type=str)
@click.option("--reference", "-r",
              help="Path to the directory containing the appropriate Bowtie index.",
              default=None,
              type=str)
@click.option("--largeindex",
              help="",
              default=False,
              required=False,
              type=bool)
@click.option("--rule_set",
              help="On-target score rule set to use.  Current options include '1', '2', and 'Azimuth'",
              default="Azimuth",
              type=str)
@click.option("--ontarget_score_threshold", "-on",
              help="Spacers with an on-target score below this will be ignored.",
              default=0,
              required=False,
              type=int)
@click.option("--offtarget_score_threshold", "-off",
              help="Spacers with an off-target score below this will be ignored.",
              default=0,
              required=False,
              type=int)
@click.option("--offtarget_count_threshold", help="Spacers with more than this many off-targets will be ignored.",
              default=100,
              required=False,
              type=int)
@click.option("--spacers_per_feature",
              help="Number of spacers to find for each feature.",
              default=6,
              type=int)
@click.option("--paired",
              help="Should spacers be designed to work as pairs (e.g. for excision)?",
              default=False,
              type=bool)
@click.option("--number_upstream_spacers",
              help="If designing paired spacers, number of spacers to design that target upstream of the feature.",
              default=3,
              type=int)
@click.option("--number_downstream_spacers",
              help="If designing paired spacers, number of spacers to design that target downstream of the feature.",
              default=3,
              type=int)
@click.option("--min_paired_distance",
              help="If designing paired spacers, minimum space required between the up- and downstream spacers.",
              default=0,
              type=int)
@click.option("--cores", "-c",
              help="Number of processors to use. By default, will use all available.",
              default=None,
              type=int)
@click.help_option()
def main():
    pass


def build_library(input_sequences: str = None,
                  outfile: str = None,
                  refgenome: str = None,
                  restriction_sites: list = None,
                  largeindex: bool = False,
                  on_target_score_threshold: int = 0,
                  off_target_score_threshold: Union[int, NoneType] = None,
                  off_target_count_threshold: Union[int, NoneType] = None,
                  nuclease: str = "SpCas9",
                  spacers_per_feature: int = 9,
                  reject: bool = False,
                  paired: bool = False,
                  rules: str = "Azimuth",
                  number_upstream_spacers: int=0,
                  number_downstream_spacers: int=0,
                  numcores: int=0):

    targets = pyfaidx.Fasta(input_sequences)
    nucleases = pd.read_csv("data/nuclease_list.csv",
                            dtype={"nuclease": str,
                                   "pam": str,
                                   "spacer_regex": str,
                                   "start": np.int8,
                                   "end": np.int8},
                            skip_blank_lines=True)
    ni = nucleases[nucleases["nuclease"] == nuclease]

    spacers_df = find_guides(itemlist=targets,
                             nuclease_info=ni,
                             restriction_sites=restriction_sites)

    initialnumber = spacers_df.shape[0]

    spacers_df = on_target_scoring(ruleset=rules,
                                   spacers=spacers_df,
                                   on_target_score_threshold=on_target_score_threshold)

    if spacers_df.shape[0] == 0:
        print("Sorry, no spacers matching that criteria were found")
        return 0
    else:
        print(f"Finished scoring spacers. {spacers_df.shape[0]} of {initialnumber} "
              f"spacers have an on-target score above the cutoff threshold of {on_target_score_threshold}."
              f"\nBeginning Bowtie alignment...")

    offtarget_results_file = off_target_discovery(spacers_df=spacers_df,
                                                  cpus = numcores,
                                                  refgenome=refgenome,
                                                  large_index_size=largeindex,
                                                  reject=reject)

    spacers_df = off_target_scoring(otrf=offtarget_results_file,
                                    spacers_df=spacers_df,
                                    nuclease_info=ni,
                                    offtarget_score_threshold=off_target_score_threshold,
                                    offtarget_count_threshold=off_target_count_threshold)

    if paired:
        guide_library = assemble_paired_library(spacers=spacers_df,
                                                on_target_threshold=on_target_score_threshold,
                                                off_target_threshold=off_target_score_threshold,
                                                number_upstream_spacers=number_upstream_spacers,
                                                number_downstream_spacers=number_downstream_spacers)
    else:
        guide_library = assemble_library(spacers=spacers_df,
                                         on_target_threshold=on_target_score_threshold,
                                         off_target_score_threshold=off_target_score_threshold,
                                         spacers_per_feature=spacers_per_feature)
    guide_library.to_csv(outfile)
    print("Finished.")


# TODO need to keep strandedness in mind when searching for Cas13 guides
def find_guides(itemlist: pyfaidx.Fasta,
                nuclease_info: pd.DataFrame,
                restriction_sites: list) -> pd.DataFrame:

    spacer_regex = compile(nuclease_info["spacer_regex"].item())
    spacer_start = nuclease_info["start"].item()
    spacer_end = nuclease_info["end"].item()

    print(f"{len(itemlist.keys())} sequences to search for spacers.")
    widgets = ["Examining sequence: ",
               progressbar.Counter(),
               " ",
               progressbar.Percentage(),
               " ",
               progressbar.Bar(),
               progressbar.Timer()]
    progress = progressbar.ProgressBar(widgets=widgets, maxval=len(itemlist)).start()

    # Set the restriction sites that we are going to make sure are not in our spacers
    rsb = RestrictionBatch(restriction_sites)

    # For each entry in the file (i.e. exonic sequence), find all of the potential protospacer sequences.
    # Return a list
    spacers_df = pd.DataFrame(columns=["gene_name","feature_id","start","stop","strand","spacer","misc"])
    for item in itemlist.keys():
        # have to use the alternative Regex module instead of Re so that findall can detect overlapping
        # sequences
        spacers = (spacer_regex.findall(itemlist[item][:].seq, overlapped=True) +
                   spacer_regex.findall(itemlist[item][:].reverse.complement.seq, overlapped=True))

        info = dict(zip(["gene_name", "feature_id", "strand", "start", "end", "misc"], item.split("_")))

        for ps in spacers:
            # Note that ps[4:24] is the actual protospacer.  I need the rest of the sequence for scoring
            ps_seq = Seq(ps[spacer_start:spacer_end], IUPAC.unambiguous_dna)
            ps_full_seq = Seq(ps, IUPAC.unambiguous_dna)

            # Get rid of anything with T(4+) as those act as RNAPIII terminators
            if "TTTT" in ps:
                # TODO Should this also eliminate anything with G(4)?
                pass
            # Get rid of anything that has the verboten restriction sites
            elif bool([y for y in rsb.search(ps_full_seq).values() if y != []]):
                pass
            # BsmBI/Esp3I is used in most of the new CRISPR vectors, especially for library construction.
            # Biopython misses potential restriction sites as it tries to match GAGACGN(5), whereas we need to find
            # matches of just the GAGACG core.  The next four lines take care of that.
            elif "GAGACG" in ps[spacer_start:spacer_end]:
                pass
            elif "CGTCTC" in ps[spacer_start:spacer_end]:
                pass
            # Eliminate potentials with a GC content <20 or >80%
            elif GC(ps_seq) <= 20 or GC(ps_seq) >= 80:
                pass
            else:
                ps_start = itemlist[item][:].seq.find(ps) + int(info["start"])
                spacer_data = {"gene_name": [info["gene_name"]],
                               "feature_id": [info["feature_id"]],
                               "start": [ps_start],
                               "stop": [ps_start + len(ps)],
                               "strand": [info["strand"]],
                               "spacer": [ps],
                               "misc": [info["misc"]]}
                _ = pd.DataFrame.from_dict(spacer_data)
                spacers_df = pd.concat([spacers_df, _])

    progress.finish()
    return spacers_df


def on_target_scoring(ruleset: str,
                      spacers: pd.DataFrame,
                      on_target_score_threshold: float) -> pd.DataFrame:

    if ruleset == 1:
        spacerlist = spacers["spacer"].tolist()
        initialnumber = len(spacers)
        print(f"Found {initialnumber} potential spacers.  Now scoring")
        sublist = []
        queue = Manager().Queue()
        pool = Pool()
        func = partial(score_entry, method=calc_score, place=queue, cutoff=on_target_score_threshold)
        mapObj = pool.map_async(func, spacerlist, callback=sublist.append)
        scoring_widgets = ["Scoring sequence: ", progressbar.Counter(), " ", progressbar.Percentage(), " ",
                           progressbar.Bar(), progressbar.Timer()]
        scoring_progress = progressbar.ProgressBar(widgets=scoring_widgets, maxval=len(spacers)).start()
        # Initialize progress
        score_position = 0
        # While the pool has not finished its task
        while not mapObj.ready():
            # Get the report from the process queue on how many spacers they have scored since we last looked
            for _ in range(queue.qsize()):
                score_position += queue.get()
                queue.task_done()
            # Give that number over to the progressbar
            scoring_progress.update(score_position)
        mapObj.wait()
        spacerscores = np.asarray([x for x in sublist[0]])
        spacers["score"] = spacerscores
        scoring_progress.finish()
    elif ruleset == "Azimuth":
        spacers["score"] = model_comparison.predict(spacers["spacer"].values)
    spacers = spacers[spacers["score"] > on_target_score_threshold]
    return spacers


def score_entry(scoring_entry: dict,
                **kwargs) -> dict:
    calc_method = kwargs["method"]
    place = kwargs["place"]
    cutoff = float(kwargs["cutoff"])
    scoring_entry["score"] = calc_method(scoring_entry["sequence"])
    # print("new entry: {}".format(sub[0]))
    # print("old entry: {}".format(scoring_entry))
    place.put(1)
    if scoring_entry["score"] > cutoff:
        return scoring_entry

if __name__ == "__main__":
    nucleases = pd.read_csv("data/nuclease_list.csv",
                            dtype={"nuclease": str,
                                   "pam": str,
                                   "spacer_regex": str,
                                   "start": np.int8,
                                   "end": np.int8},
                            skip_blank_lines=True)
    main()