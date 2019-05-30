#!/usr/bin/env python3
import pandas as pd


def assemble_library(spacers: pd.DataFrame,
                     on_target_threshold: int=100,
                     off_target_threshold: int=100,
                     spacers_per_feature: int=6) -> pd.DataFrame:

    spacers = spacers[spacers['score'] > on_target_threshold]
    spacers = spacers[spacers['offtarget_score'] > off_target_threshold]
    grouped = spacers.groupby('gene_name').apply(lambda x: x.nlargest(spacers_per_feature, 'score')).reset_index(drop=True)

    return grouped


def assemble_paired_library(spacers: pd.DataFrame,
                            on_target_threshold: int=100,
                            off_target_threshold: int=100,
                            number_upstream_spacers: int=6,
                            number_downstream_spacers: int=6) -> pd.DataFrame:

    spacers = spacers[spacers['score'] > on_target_threshold]
    spacers = spacers[spacers['offtarget_score'] > off_target_threshold]
    grouped = spacers.groupby('gene_name').apply(lambda x: x.nlargest(spacers_per_feature, 'score')).reset_index(
        drop=True)
    return grouped

def assemble_guide_list(gene_name: str,
                        spacer_df: pd.DataFrame,
                        paired: bool=False,
                        number_upstream_spacers: int=3,
                        number_downstream_spacers: int=3,
                        return_limit: int=9) -> pd.DataFrame:

    # constants for use when making paired guides
    bsmbi_arm_5 = "aaaAgcaCGAGACG"
    right_extra_spacer = "GGTTCTATGC"
    direct_repeat = "GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC"
    left_extra_spacer = "GATAGTTGCC"
    bsmbi_arm_3 = "CGTCTCGTTTTaaaa"

    all_spacers_for_gene = spacer_df[spacer_df["gene_name"] == gene_name]
    toplist = []
    not_found_list =[]

    ###--- Paired guide assembly block ---###

    if paired:  # if we are making the paired list
        # divide list into the upstream list and the downstream list
        upstream_list = [b for b in all_spacers_for_gene if "upstream" in b.GeneName]
        downstream_list = [c for c in all_spacers_for_gene if "downstream" in c.GeneName]
        # it is possible that, given the cutoff values used, we will find no spacers for either the
        # upstream or downstream portions.  Make sure that we aren"t going to try to make combinations with
        # an empty set
        if (len(upstream_list) != 0) and (len(downstream_list) != 0):
            # sort the lists, first by on-target then off-target
            ranked_upstream = sorted(upstream_list, key=attrgetter("score", "offtargetscore"), reverse=True)
            ranked_downstream = sorted(downstream_list, key=attrgetter("score", "offtargetscore"),
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
            # paired spacer here so that the existing file writing routine doesn"t shit itself
            # First, we have to make a fake header
            tempName = upstream_list[0].GeneName.split()[0]
            # print("tempName = {}".format(tempName))
            tempID = upstream_list[0].GeneID
            # print("tempID = {}".format(tempID))
            description = f"{tempID}|none|{tempName}|||0|0|0|0"
            # we need to make a dict in the form of ["description","position","score","spacer","offtargetscore"]
            tempKeys = ["description", "position", "score", "spacer", "offtargetscore"]
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
        ranked_spacers = sorted(all_spacers_for_gene, key=attrgetter("score", "offtargetscore"),
                                reverse=True)
        # if step:
        #     pass
        if return_limit == "all" or len(ranked_spacers) <= int(return_limit):
            # if we have fewer spacers than the return limit, we return everything
            for w in ranked_spacers: toplist.append(w)
        else:
            # otherwise, return a number of the top spacers than correspond to the return_limit
            for w in range(0, int(return_limit)): toplist.append(ranked_spacers[w])

    guide_dict = {"toplist": toplist, "not_found": not_found_list}

    return guide_dict