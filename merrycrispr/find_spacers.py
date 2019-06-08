import pandas as pd
import progressbar
import pyfaidx

from Bio.Restriction import RestrictionBatch
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqUtils import GC

# TODO need to keep strandedness in mind when searching for Cas13 guides
def find_spacers(itemlist: pyfaidx.Fasta,
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