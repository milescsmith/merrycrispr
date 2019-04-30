# TODO: add ability to limit guides to the template strand
# TODO: add ability to cut with Cas13a, including limiting guides to those that would bind
# mRNA and exonic only regions, maybe even exon-exon boundaries
from merrycrispr.OffTargetScoring import off_target_scoring, off_target_discovery

__author__ = 'milescsmith'
__email__ = 'mileschristiansmith@gmail.com'

import click
import gc
import pyfaidx
import pandas as pd
import numpy as np
from Bio import SeqIO
from multiprocessing import cpu_count, Manager, Pool
from typing import Union

import progressbar
from Bio.Alphabet import IUPAC, single_letter_alphabet
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from csv import writer
from functools import partial
from operator import attrgetter
from regex import compile
from subprocess import check_output

from azimuth import model_comparison
from on_target_score_calculator import calc_score


def main():
    pass

def build_library(input_sequences: str=None,
                  outfile: str=None,
                  refgenome: str=None,
                  restriction_sites: list=None,
                  largeindex: bool=False,
                  ontarget_score_threshold: int=0,
                  offtarget_score_threshold: Union[int, NoneType]=None,
                  offtarget_count_threshold: Union[int, NoneType]=None,
                  nuclease: str='Cas9',
                  return_limit: str=9,
                  reject: bool=False,
                  paired: bool=False,
                  rules: str="Azimuth",
                  numcores: int=0,
                  number_upstream_spacers: int=0,
                  number_downstream_spacers: int=0):

    targets = pyfaidx.Fasta(input_sequences)
    nucleases = pd.read_csv('data/nuclease_list.csv',
                            dtype={'nuclease': str,
                                   'pam': str,
                                   'spacer_regex': str,
                                   'start': np.int8,
                                   'end': np.int8},
                            skip_blank_lines=True)
    ni = nucleases[nucleases['nuclease'] == nuclease]

    spacers_df = find_guides(itemlist=targets,
                             nuclease_info=ni,
                             restriction_sites=restriction_sites)

    initialnumber = spacers_df.shape[0]

    spacers_df = on_target_scoring(ruleset=rules,
                                   spacers=spacers_df,
                                   ontarget_score_threshold=ontarget_score_threshold)

    if spacers_df.shape[0] == 0:
        print("Sorry, no spacers matching that criteria were found")
        return 0
    else:
        print(f"Finished scoring spacers. {spacers_df.shape[0]} of {initialnumber} "
              f"spacers have an on-target score above the cutoff threshold of {ontarget_score_threshold}."
              f"\nBeginning Bowtie alignment...")

    offtarget_results_file = off_target_discovery(spacers_df=spacers_df,
                                                  cpus = numcores,
                                                  refgenome=refgenome,
                                                  large_index_size=largeindex,
                                                  reject=reject)

    spacers_df = off_target_scoring(otrf=offtarget_results_file,
                                    spacers_df=spacers_df,
                                    nuclease_info=ni,
                                    offtarget_score_threshold=offtarget_score_threshold,
                                    offtarget_count_threshold=offtarget_count_threshold)

    ###--- Library assembly         ---###

    guides = assemble_guide_list(gene_name=z,
                                 spacer_list=formatted_spacers,
                                 paired=paired,
                                 number_upstream_spacers=number_upstream_spacers,
                                 number_downstream_spacers=number_downstream_spacers,
                                 return_limit=return_limit)



            toplist.append(guides['toplist'])
            not_found_list.append(guides['not_found'])


            print("There were {} spacers that matched your parameters.".format(len(olist)))

        print("Writing file.")
        try:
            with open(outfile, 'w') as ofile:
                output = writer(ofile, dialect='excel')
                output.writerows([headerlist])
                output.writerows(olist)
            if not_found_list:
                print(
                    "There were {} entries for which I could not find paired spacers given the current settings. They "
                    "are listed in the no_spacers_found.txt file".format(len(not_found_list)))
                with open("no_spacers_found.txt", 'w') as problem_file:
                    for problem in not_found_list:
                        problem_file.writelines(problem + '\n')
        except IOError:
            print("There is trouble writing to the file.  Perhaps it is open in another application?")
            choice = input("Would you like to try again? [y/n]")
            if choice == 'y' or choice == 'Y':
                try:
                    with open(outfile, 'w') as ofile:
                        output = writer(ofile, dialect='excel')
                        output.writerows([headerlist])
                        output.writerows(olist)
                except:
                    print("Sorry, still was unable to write to the file")
            else:
                return 0

        print("Finished.")

    except StopIteration:
        print("No spacers were found that correspond to the parameters you set.")


# TODO need to keep strandedness in mind when searching for Cas13 guides
def find_guides(itemlist: pyfaidx.Fasta,
                nuclease_info: pd.DataFrame,
                restriction_sites: list) -> pd.DataFrame:

    spacer_regex = compile(nuclease_info['spacer_regex'].item())
    spacer_start = nuclease_info['start'].item()
    spacer_end = nuclease_info['end'].item()

    print(f'{len(itemlist.keys())} sequences to search for spacers.')
    widgets = ['Examining sequence: ',
               progressbar.Counter(),
               ' ',
               progressbar.Percentage(),
               ' ',
               progressbar.Bar(),
               progressbar.Timer()]
    progress = progressbar.ProgressBar(widgets=widgets, maxval=len(itemlist)).start()

    # Set the restriction sites that we are going to make sure are not in our spacers
    rsb = RestrictionBatch(restriction_sites)

    # For each entry in the file (i.e. exonic sequence), find all of the potential protospacer sequences.
    # Return a list
    spacers_df = pd.DataFrame(columns=['gene_name','feature_id','start','stop','strand','spacer'])
    for item in itemlist.keys():
        # have to use the alternative Regex module instead of Re so that findall can detect overlapping
        # sequences
        spacers = (spacer_regex.findall(itemlist[item][:].seq, overlapped=True) +
                   spacer_regex.findall(itemlist[item][:].reverse.complement.seq, overlapped=True))

        info = dict(zip(['gene_name', 'feature_id', 'strand', 'start', 'end'], item.split("_")))

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
            elif 'GAGACG' in ps[spacer_start:spacer_end]:
                pass
            elif 'CGTCTC' in ps[spacer_start:spacer_end]:
                pass
            # Eliminate potentials with a GC content <20 or >80%
            elif GC(ps_seq) <= 20 or GC(ps_seq) >= 80:
                pass
            else:
                ps_start = itemlist[item][:].seq.find(ps) + int(info['start'])
                spacer_data = {'gene_name': [info['gene_name']],
                               'feature_id': [info['feature_id']],
                               'start': [ps_start],
                               'stop': [ps_start + len(ps)],
                               'strand': [info['strand']],
                               'spacer': [ps]}
                _ = pd.DataFrame.from_dict(spacer_data)
                spacer_df = pd.concat([spacer_df, _])

    progress.finish()
    gc.enable()
    gc.collect()

    return spacer_df


def on_target_scoring(ruleset: str,
                      spacers: pd.DataFrame,
                      ontarget_score_threshold: float) -> pd.DataFrame:

    if ruleset == 1:
        spacerlist = spacers['spacer'].tolist()
        initialnumber = len(spacers)
        print(f"Found {initialnumber} potential spacers.  Now scoring")
        sublist = []
        queue = Manager().Queue()
        pool = Pool()
        func = partial(score_entry, method=calc_score, place=queue, cutoff=ontarget_score_threshold)
        mapObj = pool.map_async(func, spacerlist, callback=sublist.append)
        scoring_widgets = ['Scoring sequence: ', progressbar.Counter(), ' ', progressbar.Percentage(), ' ',
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
        spacers['score'] = spacerscores
        scoring_progress.finish()
    elif ruleset == "Azimuth":
        spacers['score'] = model_comparison.predict(spacers['spacer'].values)
    spacers = spacers[spacers_df['score'] >  ontarget_score_threshold]
    return spacers


def score_entry(scoring_entry: dict,
                **kwargs) -> dict:
    calc_method = kwargs['method']
    place = kwargs['place']
    cutoff = float(kwargs['cutoff'])
    scoring_entry['score'] = calc_method(scoring_entry['sequence'])
    # print("new entry: {}".format(sub[0]))
    # print("old entry: {}".format(scoring_entry))
    place.put(1)
    if scoring_entry['score'] > cutoff:
        return scoring_entry


def assemble_guide_list(gene_name: str,
                        spacer_df: pd.DataFrame,
                        paired: bool=False,
                        number_upstream_spacers: int=3,
                        number_downstream_spacers: int=3,
                        return_limit: str = 9) -> dict:

    # constants for use when making paired guides
    bsmbi_arm_5 = 'aaaAgcaCGAGACG'
    right_extra_spacer = 'GGTTCTATGC'
    direct_repeat = 'GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC'
    left_extra_spacer = 'GATAGTTGCC'
    bsmbi_arm_3 = 'CGTCTCGTTTTaaaa'

    all_spacers_for_gene = spacer_df[spacer_df['gene_name']== gene_name]
    toplist = []
    not_found_list =[]

    ###--- Paired guide assembly block ---###

    if paired:  # if we are making the paired list
        # divide list into the upstream list and the downstream list
        upstream_list = [b for b in all_spacers_for_gene if 'upstream' in b.GeneName]
        downstream_list = [c for c in all_spacers_for_gene if 'downstream' in c.GeneName]
        # it is possible that, given the cutoff values used, we will find no spacers for either the
        # upstream or downstream portions.  Make sure that we aren't going to try to make combinations with
        # an empty set
        if (len(upstream_list) != 0) and (len(downstream_list) != 0):
            # sort the lists, first by on-target then off-target
            ranked_upstream = sorted(upstream_list, key=attrgetter('score', 'offtargetscore'), reverse=True)
            ranked_downstream = sorted(downstream_list, key=attrgetter('score', 'offtargetscore'),
                                       reverse=True)
            # and then grab the top user-defined number from each
            if len(ranked_upstream) > number_upstream_spacers:
                top_upstream_list = [ranked_upstream[i] for i in range(0, number_upstream_spacers)]
            else:
                top_upstream_list = ranked_upstream
            if len(ranked_downstream) > number_downstream_spacers:
                top_downstream_list = [ranked_downstream[j] for j in range(0, number_downstream_spacers)]
            else:
                top_downstream_list = ranked_downstream
            # make permutations of the top upstream * downstream
            permutations = [[x.seq, y.seq] for x in top_upstream_list for y in top_downstream_list]
            # then we add on the extra bases that are necessary for Cas9 to process the pseudo-array and to
            # anneal into our vector
            combinations = [bsmbi_arm_5 + right_extra_spacer + str(
                combo[0]).lower() + direct_repeat + left_extra_spacer + str(combo[1]).lower()
                            + bsmbi_arm_3 for combo in permutations]
            # each combination is a new entity entirely so (and essentially just a string representing the
            # combined spacer sequence so we have to create new entries in the type FormattedResult for each
            # paired spacer here so that the existing file writing routine doesn't shit itself
            # First, we have to make a fake header
            tempName = upstream_list[0].GeneName.split()[0]
            # print("tempName = {}".format(tempName))
            tempID = upstream_list[0].GeneID
            # print("tempID = {}".format(tempID))
            description = f"{tempID}|none|{tempName}|||0|0|0|0"
            # we need to make a dict in the form of ['description','position','score','spacer','offtargetscore']
            tempKeys = ['description', 'position', 'score', 'spacer', 'offtargetscore']
            for w in combinations:
                tempValues = [description, 0, 0, w, 0]
                v = dict(zip(tempKeys, tempValues))
                print(w)
                print(FormattedResult(v))
                toplist.append(FormattedResult(v))
        # If there are no upstream or downstream spacers, we should let the user know for what genes this
        # was encountered
        else:
            not_found_list.append(all_spacers_for_gene[0].GeneName.split()[0])

    ###--- End paired guide assembly block ---###

    ###--- Non-paired guide assembly block ---###
    else:  # for everything else...
        # sort the spacers, first by the on-target score then by the off-target score
        ranked_spacers = sorted(all_spacers_for_gene, key=attrgetter('score', 'offtargetscore'),
                                reverse=True)
        # if step:
        #     pass
        if return_limit == 'all' or len(ranked_spacers) <= int(return_limit):
            # if we have fewer spacers than the return limit, we return everything
            for w in ranked_spacers: toplist.append(w)
        else:
            # otherwise, return a number of the top spacers than correspond to the return_limit
            for w in range(0, int(return_limit)): toplist.append(ranked_spacers[w])

    guide_dict = {'toplist': toplist, 'not_found': not_found_list}

    return guide_dict
