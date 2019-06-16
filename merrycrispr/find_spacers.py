from typing import Optional, List

import pandas as pd
import pyfaidx
import regex
from Bio.Alphabet import IUPAC
from Bio.Restriction import RestrictionBatch
from Bio.Seq import Seq
from Bio.SeqUtils import GC


def find_spacers(
    itemlist: pyfaidx.Fasta,
    nuclease_info: dict,
    restriction_sites: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Find protospacers in a sequence for a given nuclease.  The search region can be more
    expansive than just the PAM (for scoring purposes) and strand can be taken into account.
    `find_spacers()` will ignore sequences that have a poly(T) sequence, high GC content, or
    a motif matching given restriction nuclease sequence.

    Parameters
    __________
    itemlist : :class:`~pyfaidx.Fasta`
        Parsed FASTA with sequences to examine for spacers
    nuclease_info : dict
        Information for the nuclease to use.  Required keys include `start`, `end`, `strand`,
        `pam`, `spacer_regex`
    restriction_sites : `List[str]`, optional (default: `None`)
        For a found spacer, ignore it if it contains the sequence for recognition by these
        restriction endonucleases.


    Return
    ______
    :class:`~pandas.DataFrame`

    """
    spacer_regex: regex.Regex = regex.compile(nuclease_info["spacer_regex"])
    spacer_start: int = nuclease_info["start"]
    spacer_end: int = nuclease_info["end"]

    # Set the restriction sites that we are going to make sure are not in our
    # spacers
    if restriction_sites:
        rsb = RestrictionBatch(restriction_sites)

    # For each entry in the file (i.e. exonic sequence), find all of the
    # potential protospacer sequences.
    spacers_df = pd.DataFrame(
        columns=[
            "gene_name",
            "feature_id",
            "strand",
            "start",
            "end",
            "spacer",
            "seq_hash",
        ]
    )
    for item in itemlist.keys():
        # have to use the alternative Regex module instead of Re so that findall
        # can detect overlapping sequences
        # Since Cas13/Csc2 recognizes mRNA, we need to keep strandedness in mind
        # go ahead and check for all three possibilites to future proof for
        # other nucleases.
        try:
            if nuclease_info["strand"] == "b":
                spacers = spacer_regex.findall(
                    itemlist[item][:].seq, overlapped=True
                ) + spacer_regex.findall(
                    itemlist[item][:].reverse.complement.seq, overlapped=True
                )
            elif nuclease_info["strand"] == "p":
                spacers = spacer_regex.findall(itemlist[item][:].seq, overlapped=True)
            elif nuclease_info["strand"] == "n":
                spacers = spacer_regex.findall(
                    itemlist[item][:].reverse.complement.seq, overlapped=True
                )
        except ValueError:
            print(
                f"The current nuclease is set to recognize {nuclease_info['strand']} strand, "
                f"which is not one I recognize"
            )

        info = dict(
            zip(
                ["gene_name", "feature_id", "strand", "start", "end", "seq_hash"],
                item.split("_"),
            )
        )

        for ps in spacers:
            # Note that ps[4:24] is the actual protospacer.  I need the rest of
            # the sequence for scoring
            ps_seq = Seq(ps[spacer_start:spacer_end], IUPAC.unambiguous_dna)
            ps_full_seq = Seq(ps, IUPAC.unambiguous_dna)

            # Get rid of anything with T(4+) as those act as RNAPIII
            # terminators
            if "TTTT" in ps:
                # TODO Should this also eliminate anything with G(4)?
                pass
            # Get rid of anything that has the verboten restriction sites
            elif bool([y for y in rsb.search(ps_full_seq).values() if y != []]):
                pass
            # BsmBI/Esp3I is used in most of the new CRISPR vectors, especially for
            # library construction. Biopython misses potential restriction sites as it tries
            # to match GAGACGN(5), whereas we need to find matches of just the GAGACG core.
            # The next four lines take care of that.
            elif "GAGACG" in ps[spacer_start:spacer_end]:
                pass
            elif "CGTCTC" in ps[spacer_start:spacer_end]:
                pass
            # Eliminate potentials with a GC content <20 or >80%
            elif GC(ps_seq) <= 20 or GC(ps_seq) >= 80:
                pass
            else:
                ps_start = itemlist[item][:].seq.find(ps) + int(info["start"])
                spacer_data = {
                    "gene_name": [info["gene_name"]],
                    "feature_id": [info["feature_id"]],
                    "strand": [info["strand"]],
                    "start": [ps_start],
                    "end": [ps_start + len(ps)],
                    "spacer": [ps],
                    "seq_hash": [info["seq_hash"]],
                }
                _ = pd.DataFrame.from_dict(spacer_data)
                spacers_df = pd.concat([spacers_df, _])

    return spacers_df
