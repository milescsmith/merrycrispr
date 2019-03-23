# TODO: add ability to limit guides to the template strand
# TODO: add ability to cut with Cas13a, including limiting guides to those that would bind mRNA and exonic only regions, maybe even exon-exon boundaries

__author__ = 'milescsmith'
__email__ = 'mileschristiansmith@gmail.com'

import click
import gc
import pyfaidx
import pandas as pd
import numpy as np
from Bio import SeqIO
from multiprocessing import cpu_count, Manager, Pool

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


class FormattedResult(object):
    """ Class to automate the formatting of the final protospacers
    """
    __slots__ = ('GeneID', 'ExonRank', 'GeneName', 'seq', 'GC', 'position', 'score', 'offtargetscore')

    # noinspection PyPep8
    def __init__(self, recorddict):
        # ['description','position','score','spacer','offtargetscore']
        # description = 100767553|Unknown|Man1a1|||||1838438|1838994 1838822 0.111881537716
        # GENEID | TRANSCRIPTID | GENENAME | EXON RANK | CONSTITUTIVE EXON | 5' UTR END | 3' UTR STOP | EXON START | EXON END
        # GENEID | GENENAME | EXON RANK | SEQ_START | SEQ_END
        tempdescript = recorddict['description'].split('|')
        # self.GeneID = tempdescript[0]
        # self.ExonRank = tempdescript[3]  # Ensembl ExonIDs are nearly useless for sorting as they not in numerical order
        # self.GeneName = tempdescript[2]
        # self.seq = str(recorddict['spacer'])
        # self.GC = GC(self.seq)
        # self.position = recorddict['position']
        # self.score = recorddict['score']
        # self.offtargetscore = recorddict['offtargetscore']
        self.GeneID = tempdescript[0]
        self.GeneName = tempdescript[1]
        self.seq = str(recorddict['spacer'])
        self.GC = GC(self.seq)
        self.position = recorddict['position']
        self.score = recorddict['score']
        self.offtargetscore = recorddict['offtargetscore']

    def __repr__(self):
        return 'GeneID: {self.GeneID}, GeneName: {self.GeneName}, ExonRank: {self.ExonRank}, Sequence: {self.seq}, ' \
               'Position: {self.position}, GC%: {self.GC}, score: {self.score}, ' \
               'off-target score: {self.offtargetscore}'.format(self=self)
#
#
# class m_SeqRecord(SeqRecord):
#     """An altered form of SeqRecord that reads in and parses a header that contains the 5' UTR end and 3' UTR start
#     information and uses that to trim the sequence
#
#     From BioMart, you need to select the following Attributes in this order:
#     Ensembl Gene ID
#     Ensembl Transcript ID
#     Exon sequences
#     Associated Gene Name
#     Exon Rank in Transcript
#     Constitutive Exon
#     5' UTR End
#     3' UTR Start
#     Exon Chr Start (bp)
#     Exon Chr End (bp)
#
#     If you are using GFFextractor, everything should already be formatted correctly
#     """
#
#     __slots__ = ('id', 'tempseq', 'GeneID', 'TranscriptID', 'GeneName', 'ExonRank', 'ConstExon',
#                  'FPUTRend', 'TPUTRstart', 'exonStart', 'exonEnd')
#
#     def __init__(self, seq, id="<unknown id>", name="<unknown name>", description="<unknown description>"):
#         # Essentially, all the information we need other than the sequence is passed as the 'id',
#         # with each portion separated by a '|', in the order defined above
#         self.id = id
#         self.tempseq = seq
#         temp = self.id.split('|')
#         self.GeneID = temp[0]
#         self.TranscriptID = temp[1]
#         self.GeneName = temp[2]
#         self.ExonRank = temp[3]
#         if temp[4] == '':
#             self.ConstExon = False
#         else:
#             self.ConstExon = True
#         try:
#             self.FPUTRend = int(temp[5])
#         except:
#             self.FPUTRend = ''
#         try:
#             self.TPUTRstart = int(temp[6])
#         except:
#             self.TPUTRstart = ''
#         try:
#             self.exonStart = int(temp[7])
#         except:
#             self.exonStart = ''
#         try:
#             self.exonEnd = int(temp[8])
#         except:
#             self.exonEnd = ''
#
#         # If there is a 5'UTR
#         if self.FPUTRend != '':
#             # Cut off an amount of the sequence that is equal to the length of the UTR
#             sequence = str(self.tempseq)[(self.FPUTRend - self.exonStart):(self.exonEnd - self.exonStart)]
#         # The same, but do it for a 3'UTR
#         elif self.TPUTRstart != '':
#             sequence = str(self.tempseq)[0:(self.TPUTRstart - self.exonStart)]
#         # Otherwise, leave the sequence alone
#         else:
#             sequence = str(self.tempseq)
#
#         SeqRecord.__init__(self, Seq(sequence, IUPAC.unambiguous_dna), id, name, description, dbxrefs=None,
#                            features=None, annotations=None, letter_annotations=None)
#
#
# def m_FastaIterator(handle, alphabet=single_letter_alphabet):
#     """ FastaIterator from the Bio.SeqIO.FastaIO module altered to use the m_SeqRecord class
#       Generator function to iterate over Fasta records (as SeqRecord objects).
#
#       Arguments:
#        - handle - input file
#        - alphabet - optional alphabet
#        - title2ids - A function that, when given the title of the FASTA
#          file (without the beginning >), will return the id, name and
#          description (in that order) for the record as a tuple of strings.
#          If this is not given, then the entire title line will be used
#          as the description, and the first word as the id and name.
#          :param handle:
#          :param alphabet:
#       """
#
#     from Bio.SeqIO.FastaIO import SimpleFastaParser
#     for title, sequence in SimpleFastaParser(handle):
#         try:
#             first_word = title.split(None, 1)[0]
#         except IndexError:
#             assert not title, repr(title)
#             first_word = ""
#         yield m_SeqRecord(Seq(sequence, alphabet), id=first_word, name=first_word, description=title)


def main():


def build_library(input_sequences: str=None,
                  outfile: str=None,
                  refgenome: str=None,
                  restriction_sites: list=None,
                  largeindex: bool=False,
                  cutoff: int=0,
                  offtargetcutoff: int=0,
                  nuclease: str='Cas9',
                  return_limit: str=9,
                  reject: bool=False,
                  paired: bool=False,
                  rules: int=2,
                  numcores: int=0,
                  number_upstream_spacers: int=0,
                  number_downstream_spacers: int=0):

    targets = pyfaidx.Fasta(input_sequences)
    nucleases = pd.read_csv('data/nuclease_list.csv')
    ni = nucleases[nucleases['nuclease'] == nuclease]

    # Setup the regular expression for the pam site and spacer
    # if nuclease == 'Cas9':
    #     # Matches a 20N-NGG target sequence plus the -4->-3 and +1->+3 nucleotides needed for scoring
    #     pam = r'(?i)[ACGT]{25}[G]{2}[ACGT]{3}'
    # elif nuclease == 'Cpf1':
    #     # Matches a 30N-NGG-NNN target
    #     pam = r'(?i)[T]{2,}[A-Z]{25}'
    # elif (nuclease is 'Cas13a') or (nuclease is 'Csc2'):
    #     pam = r'[ATGC]{21}[ATC]{3}' # Setting this at 21 for the moment as the Abudayyeh et al. Science 550 (2018) paper
    #                                 # showed efficient cutting down to 20 nt.
    #                                 # TODO: add option to make this longer?

    ###--- Find potential guides block      ---###

    spacerlist = find_guides(itemlist=targets,
                             nuclease_info=ni,
                             restriction_sites=restriction_sites)

    ###--- End potential guide search block ---###
    ###--- On-target scoring block          ---###

    spacerlist = on_target_scoring(ruleset=rules, spacers=spacerlist, cutoff=cutoff)

    ###--- End scoring block ---###

    if len(spacerlist) == 0:
        print("Sorry, no spacers matching that criteria were found")
        return 0
    else:
        print("Finished scoring spacers.  {} of {} spacers have an on-target score above the cutoff threshold of "
              "{}.\nBeginning Bowtie alignment...".format(len(spacerlist), initialnumber, cutoff))

    ###--- Off-target discovery ---###

    offtarget_results_file = off_target_discovery(spacerlist=spacerlist,
                                                  cpus = numcores,
                                                  refgenome=refgenome,
                                                  large_index_size=largeindex,
                                                  reject=reject)

    ###--- End off-target discovery ---###
    ###--- Off-target scoring       ---###

    prunedlist = off_target_scoring(otrf=offtarget_results_file,
                                    nuclease=nuclease,
                                    offtargetcutoff=offtargetcutoff,
                                    spacerlist=spacerlist,
                                    reject=reject)

    ###--- End off-target scoring   ---###
    ###--- Library assembly         ---###
    try:
        not_found_list = []
        if len(prunedlist[0]['description'].split('|')) == 4:  # need to make sure the header format is correct
            formatted_spacers = [FormattedResult(entry) for entry in prunedlist]
            # Create a set of all the GeneIDs in our list
            geneset = {y.GeneName for y in formatted_spacers}
            toplist = []

            for z in geneset:
                guides = assemble_guide_list(gene_name=z,
                                             spacer_list=formatted_spacers,
                                             paired=paired,
                                             number_upstream_spacers=number_upstream_spacers,
                                             number_downstream_spacers=number_downstream_spacers,
                                             return_limit=return_limit)

            toplist.append(guides['toplist'])
            not_found_list.append(guides['not_found'])


            # print('toplist is {} long'.format(len(toplist)))
            olist = [[entry.GeneID, entry.GeneName, entry.seq, entry.GC, entry.position,
                      entry.score, entry.offtargetscore] for entry in toplist]
            headerlist = ['GeneID', 'GeneName', 'seq', '%GC', 'position', 'score', 'off-target score']
            print("There were {} spacers that matched your parameters.".format(len(olist)))
        else:  # if the header format it isn't in the correct format (though I don't know how that would happen),
            # just dump all the spacers we found into a file
            finallist = spacerlist
            olist = [[entry.id.split(' ')[0], entry.seq, GC(entry.seq), entry.position.split(' ')[1],
                      entry.score.split(' ')[2]] for entry in finallist]
            headerlist = ['ID', 'seq', '%GC', 'position', 'score']
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
    spacer_df = pd.DataFrame(columns=['gene_name','feature_id','start','stop','strand','spacer'])
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
                      cutoff: float) -> pd.DataFrame:

    if ruleset == 1:
        spacerlist = spacers['spacer'].tolist()
        initialnumber = len(spacers)
        print(f"Found {initialnumber} potential spacers.  Now scoring")
        sublist = []
        queue = Manager().Queue()
        pool = Pool()
        func = partial(score_entry, method=calc_score, place=queue, cutoff=cutoff)
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
        spacers['score'] = model_comparison.predict(spacers['spacer'].as_matrix())
    spacers = spacers[spacer_df['score'] > cutoff]
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

def off_target_discovery(spacerlist: list,
                         cpus: int = 0,
                         refgenome: str = "",
                         large_index_size: bool = False,
                         reject: bool = False) -> str:

    if cpus is 0:
        cpus = cpu_count()

    # The below would be great if it actually worked.  Unfortunately, Bowtie will not allow us to set right-most
    # nucleotides as an absolute requirement, only left-most.  I think I tried doing something where I flipped
    # everthing around so that we could do that, but something horrible happened.
    '''
    # add the 'NGG' to the end of the spacer so we can restrict off-target analysis to those only those sites that
    # have a pam
    with open('temp.fa', 'w') as tempfile:
        for entry in spacerlist:
            tempfile.writelines(">%s %s %s\n%s\n" % (entry['description'], entry['position'],
                                                     entry['score'], (entry['spacer']+'NGG')))
    '''
    # TODO: this has to be changed to work with the longspacer above
    with open('temp.fa', 'w') as tempfile:
        for entry in spacerlist:
            tempfile.writelines(f">{entry['description']} {entry['position']} {entry['score']}\n{entry['spacer']}\n"
                                )

    # delete lists we are not going to use in the future to save on some memory


    # off-target identification with Bowtie
    # Use Bowtie to find if this particular sequence has any potential off targets (i. e. two or fewer mismatches)
    # At current, Bowtie is set to return everything that a particular spacer matches within the reference genome
    # with two or fewer mismatches.
    # TODO switch to Bowtie2 so we can account for gaps in mismatches
    # The problem with Bowtie2 is that it outputs in the SAM format.  Biopython does not have a SAM parser, though
    # there is another library, pysam, that will read it.  This, of course, complicates things.  Bowtie, on the
    # otherhand, will output to a plaintext file we can read.
    # TODO modify this setup to account for NAG PAMs or do we even need to? (i. e. are we already?)
    program = 'bowtie'
    cpus = "-p" + str(cpus)
    # So, if you go to tinker with this, note that there can be NO spaces in any of the strings you pass.  Maybe if
    # you pass a variable (like 'refgenome'). The '--suppress' below is a good example - that is one argument passed to
    # Bowtie, "--suppress 3,5,6,7", but everything goes all to hell unless you pass those as two separate arguments.
    if reject:
        if large_index_size:
            check_output([program, '-a', "--suppress", "3,5,6,7", cpus, '-m', reject, '--large-index',
                          refgenome, '-f', 'temp.fa', 'offtargets.fa'])
        else:
            check_output([program, '-a', "--suppress", "3,5,6,7", cpus, '-m', reject, refgenome, '-f',
                          'temp.fa', 'offtargets.fa'])
    else:
        if large_index_size:
            check_output([program, '-a', "--suppress", "3,5,6,7", cpus, '--large-index',
                          refgenome, '-f', 'temp.fa', 'offtargets.fa'])
        else:
            check_output([program, '-a', "--suppress", "3,5,6,7", cpus, refgenome, '-f',
                          'temp.fa', 'offtargets.fa'])

    print("Bowtie finished.  Begining offtarget analysis...")
    bowtie_results_file = 'offtargets.fa'

    return bowtie_results_file


def scoreofftarget(mismatched_positions: list) -> float:
    """ Calculate the likelihood a protospacer will cut at a particular off-target site
    Equation from http://crispr.mit.edu/about

    The mismatch scoring algorithm from the Zhang group has three terms:
    1) the effect of a mismatch at a particular position
    2) the product of the mismatch scores, weighted by the mean distance between each mismatch
    3) a penalty for the number of mismatches
    Score is from 0 to 1, with higher scores being good.
    """

    # experimentally determined weighting penality for mismatch at each position
    M = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0, 0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.731, 0.828,
         0.615, 0.804, 0.685, 0.583]
    if len(mismatched_positions) == 1:  # if there is only one mismatch, we should ignore the second two terms.
        score = 1 - M[mismatched_positions[0]]
        return score
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


def sumofftargets(offtargetlist: list) -> float:

    sum_score = sum((scoreofftarget(x) for x in offtargetlist))
    return 1.0 / (1 + sum_score) * 100


def off_target_scoring(otrf: str,
                       nuclease: str = 'Cas9',
                       offtargetcutoff: int = 0,
                       spacerlist: list = None,
                       reject: bool = False) -> list:

    if nuclease is 'Cas9':
        # TODO: adjust this to account for other PAMs
        mmpos = '[0-9]{1,}'
        mmpos_re = compile(mmpos)
        # these are so we ignore the pam when looking for mismatches
        reject_pos = [20, 21, 22]

    total_lines = check_output(['wc', '-l', otrf])
    print(f"Total alignments from Bowtie: {total_lines}")

    oftcount = 1
    with open(otrf) as offtargetsfile:  # parse all the bowtie results into something we can use
        # Read in the first line and parse.
        keys = ['readname', 'strand', 'position', 'mmpositions']
        lines = (line.strip('\n').split('\t') for line in offtargetsfile)
        listofspacers = {f'{spacer["description"]} {spacer["position"]} {str(spacer["score"])}': spacer
                         for spacer in spacerlist}

        bowtie_widgets = ['Processing off-targets. Examining: ', progressbar.Counter(), ' ', progressbar.Percentage(),
                          ' ', progressbar.Bar(), progressbar.Timer()]
        bowtie_progress = progressbar.ProgressBar(widgets=bowtie_widgets, maxval=total_lines).start()

        for line in lines:
            bowtie_progress.update(oftcount)  # update the position on the progressbar
            oftcount += 1  # update the count used by the progressbar
            entry = dict(zip(keys, line))
            # process the mismatch positions
            perfectmatch = 0  # variable that keeps track of how many mismatches are between the guide and whatever
            # Bowtie found
            pos = mmpos_re.findall(entry['mmpositions'])
            mmlist = None
            if len(pos) == 0:  # Bowtie returns a match for the spacer itself in the genome
                # if we have two such matches, we should indicate that there is a perfect off-target
                perfectmatch += 1
            else:
                # Make a list of where the mismatches are, as long as they are in the first 20 nucleotides.
                # We don't care what is in the mismatch position, only WHERE they are.
                temp = [int(w) for w in pos if int(w) not in reject_pos]
                if temp:
                    mmlist = temp
            # find the entry in list of spacers to which it belongs
            # add the list of mismatches to that entry
            if perfectmatch > 1:
                listofspacers[entry['readname']]['otpositions'].append('discard')
            elif mmlist:
                listofspacers[entry['readname']]['otpositions'].append(mmlist)
        bowtie_progress.finish()

        matching_widgets = ['Matching off-targets. Examining: ', progressbar.Counter(), ' ', progressbar.Percentage(),
                            ' ', progressbar.Bar(), progressbar.Timer()]
        matching_progress = progressbar.ProgressBar(widgets=matching_widgets, maxval=len(listofspacers)).start()
        matching_pos = 0

        for key in listofspacers:
            matching_pos += 1
            matching_progress.update(matching_pos)
            discard = False
            # Tally up the offtarget score for the set of off targets and set the spacer's off-target score
            for x in listofspacers[key]['otpositions']:
                if 'discard' in x:  # just skip over anything we found has more than one other perfect match.
                    discard = True
            if len(listofspacers[key]['otpositions']) != 0 or discard:
                if discard:
                    listofspacers[key]['offtargetscore'] = 0
                # Speed this up by rejecting things with over a certain set of matching off-targets
                elif reject and len(listofspacers[key]['otpositions']) > reject:
                    listofspacers[key]['offtargetscore'] = 0
                else:
                    listofspacers[key]['offtargetscore'] = sumofftargets(listofspacers[key]['otpositions'])

        matching_progress.finish()

    print("\nFinished scoring off-targets.")

    # There is no need to process spacers that have an off-target score below the cutoff value
    # so we make a new list with which we can work that is populated with the good stuff
    prunedlist = [listofspacers[spacer] for spacer in listofspacers if
                  float(listofspacers[spacer]['offtargetscore']) >= float(offtargetcutoff)]

    return prunedlist


def assemble_guide_list(gene_name: str,
                        spacer_list: list,
                        paired: bool,
                        number_upstream_spacers: int = 3,
                        number_downstream_spacers: int = 3,
                        return_limit: str = 9) -> dict:

    # constants for use when making paired guides
    bsmbi_arm_5 = 'aaaAgcaCGAGACG'
    right_extra_spacer = 'GGTTCTATGC'
    direct_repeat = 'GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC'
    left_extra_spacer = 'GATAGTTGCC'
    bsmbi_arm_3 = 'CGTCTCGTTTTaaaa'

    all_spacers_for_gene = [_ for _ in spacer_list if _.GeneName == gene_name]
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


def parse_item_info(item: str) -> dict:

    item = dict(zip(['gene_name','feature_id','strand'_'start','end'], item.split("_")))
