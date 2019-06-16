#!/usr/bin/env python3
from itertools import product

import pandas as pd

# the capitals below are to make it easier to visually discern the arm/extra spacer boundaries
BSMBI_ARM_5 = "aaaagcacgagacG"
RIGHT_EXTRA_SPACER = "ggttctatgc"
DIRECT_REPEAT = "gttttagagctatgctgttttgaatggtcccaaaac"
LEFT_EXTRA_SPACER = "gatagttgcc"
BSMBI_ARM_3 = "Cgtctcgttttaaaa"


def assemble_library(
    spacers: pd.DataFrame,
    on_target_score_threshold: int = 100,
    off_target_score_threshold: int = 100,
    spacers_per_feature: int = 6,
) -> pd.DataFrame:
    """Creates a final list of protospacers for synthesis

    Parameters
    __________
    spacers : :class:`~pd.DataFrame`
        Dataframe with all spacers found by :module:`~find_spacers.find_spacers`,
        scores added by :module:`~on_target_scoring.on_target_scoring` and
        :module:`~off_target_scoring.off_target_scoring`
    on_target_score_threshold : int, optional (default: 100)
        Spacers with an on-target score below this threshold will be removed
    off_target_score_threshold : int, optional (default: 100)
        Spacers with an off-target score below this threshold will be removed
    spacers_per_feature : int, optional (default: 6)
        The number of spacers to return for each gene

    Return
    ______
    :class:`~pd.DataFrame` with the final spacer sequences for synthesis
    """

    spacers = spacers[spacers["on_target_score"] > on_target_score_threshold]
    spacers = spacers[spacers["off_target_score"] > off_target_score_threshold]
    grouped = (
        spacers.groupby("gene_name")
        .apply(lambda x: x.nlargest(spacers_per_feature, "on_target_score"))
        .reset_index(drop=True)
    )

    return grouped


def assemble_paired_library(
    spacers: pd.DataFrame,
    on_target_score_threshold: int = 100,
    off_target_score_threshold: int = 100,
    number_upstream_spacers: int = 3,
    number_downstream_spacers: int = 3,
    mix_and_match: bool = True,
) -> pd.DataFrame:
    """Creates a final list of protospacers for synthesis.  Used to create excision libraries,
    where two spacers are necessary to cause cuts at either side of a feature.
    `assemble_paired_library()` will take a set of upstream and set of downstream spacers,
    generate all permutations for those originating for the same feature, and assemble them in a
    synthetic SpCas9 spacer array

    Parameters
    __________
     : :class:`~pd.DataFrame`
        Dataframe with all spacers found by :module:`~find_spacers.find_spacers`,
        scores added by :module:`~on_target_scoring.on_target_scoring` and
        :module:`~off_target_scoring.off_target_scoring`
    on_target_score_threshold : int, optional (default: 100)
        Spacers with an on-target score below this threshold will be removed
    off_target_score_threshold : int, optional (default: 100)
        Spacers with an off-target score below this threshold will be removed
    number_upstream_spacers : int, optional (default: 3)
        Number of spacers upstream of a gene to use
    number_downstream_spacers : int, optional (default: 3)
        Number of spacers upstream of a gene to use
    mix_and_match : bool, optional (default: True)
        If `True`, permutations of the final upstream and downstream spacers will be assebled
        into a larger synthetic spacer array construct.

    Return
    ______
    :class:`~pd.DataFrame` with the final spacer sequences for synthesis.  If `mix_and_match` is
    `True`, then this will correspond to the spacer arrays; if `False`, then this will be a
    listing of the final upstream and downstream spacers.
    """

    spacers = spacers[spacers["on_target_score"] > on_target_score_threshold]
    spacers = spacers[spacers["off_target_score"] > off_target_score_threshold]
    upstream_spacers = spacers[spacers["gene_name"].str.contains("upstream")]
    downstream_spacers = spacers[spacers["gene_name"].str.contains("downstream")]

    grouped_upstream = (
        upstream_spacers.groupby("seq_hash")
        .apply(lambda x: x.nlargest(number_upstream_spacers, "on_target_score"))
        .reset_index(drop=True)
    )
    grouped_downstream = (
        downstream_spacers.groupby("seq_hash")
        .apply(lambda x: x.nlargest(number_downstream_spacers, "on_target_score"))
        .reset_index(drop=True)
    )

    if mix_and_match:
        original_targets = spacers["seq_hash"].drop_duplicates().values

        combo_df = pd.DataFrame(
            columns=[
                "gene_name",
                "feature_id",
                "strand",
                "spacer",
                "upstream_on_target_score",
                "downstream_on_target_score",
                "upstream_off_target_score",
                "downstream_off_target_score",
                "seq_hash",
                "upstream_hash",
                "downstream_hash"
            ]
        )

        for _ in original_targets:
            tmp_upstream_spacers = grouped_upstream[grouped_upstream["seq_hash"] == _]
            tmp_downstream_spacers = grouped_downstream[
                grouped_downstream["seq_hash"] == _
            ]

            for permuted_indices in product(
                tmp_upstream_spacers.index, tmp_downstream_spacers.index
            ):
                upstream_index, downstream_index = permuted_indices
                instance_df = pd.DataFrame(
                    {
                        "gene_name": tmp_upstream_spacers["gene_name"]
                        .drop_duplicates()
                        .item()
                        .strip("-upstream"),
                        "feature_id": tmp_upstream_spacers["feature_id"],
                        "strand": tmp_upstream_spacers["strand"],
                        "spacer": "".join(
                            [
                                BSMBI_ARM_5,
                                RIGHT_EXTRA_SPACER,
                                tmp_upstream_spacers.loc[upstream_index, "spacer"],
                                DIRECT_REPEAT,
                                LEFT_EXTRA_SPACER,
                                tmp_downstream_spacers.loc[downstream_index, "spacer"],
                                BSMBI_ARM_3,
                            ]
                        ),
                        "upstream_on_target_score": tmp_upstream_spacers.loc[
                            upstream_index, "on_target_score"
                        ],
                        "downstream_on_target_score": tmp_downstream_spacers.loc[
                            downstream_index, "on_target_score"
                        ],
                        "upstream_off_target_score": tmp_upstream_spacers.loc[
                            upstream_index, "off_target_score"
                        ],
                        "downstream_off_target_score": tmp_downstream_spacers.loc[
                            downstream_index, "off_target_score"
                        ],
                        "seq_hash": tmp_upstream_spacers.loc[
                            upstream_index, "seq_hash"
                        ],
                        "upstream_hash": tmp_upstream_spacers.loc[
                            upstream_index, "hash"
                        ],
                        "downstream_hash": tmp_downstream_spacers.loc[
                            downstream_index, "hash"
                        ],
                    }
                )
                combo_df = pd.concat([combo_df, instance_df]).drop_duplicates()
        return combo_df
    else:
        return pd.concat([grouped_upstream, grouped_downstream])
