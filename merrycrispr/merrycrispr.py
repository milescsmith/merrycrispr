#!/usr/bin/env python3

# TODO: add ability to limit guides to the template strand
# TODO: add ability to cut with Cas13a, including limiting guides to those that would bind
# mRNA and exonic only regions, maybe even exon-exon boundaries

__author__ = "milescsmith"
__email__ = "mileschristiansmith@gmail.com"

from typing import Optional, List

import click
import numpy as np
import pandas as pd
import pyfaidx

from .find_spacers import find_spacers
from .library_assembly import assemble_paired_library, assemble_library
from .off_target_scoring import off_target_discovery, off_target_scoring
from .on_target_scoring import on_target_scoring

NUCLEASES = pd.read_csv(
    "data/nuclease_list.csv",
    dtype={
        "nuclease": str,
        "pam": str,
        "spacer_regex": str,
        "start": np.int8,
        "end": np.int8,
        "strand": str,
    },
    skip_blank_lines=True,
)

AVAILABLE_NUCLEASES = ", ".join(NUCLEASES["nuclease"])


@click.command()
@click.option(
    "--input",
    "-i",
    help=f"Input FASTA file containing sequences to target.",
    default=None,
    required=True,
    type=str,
)
@click.option(
    "--output",
    "-p",
    help="Name of file to write library to (in CSV format).",
    default=None,
    required=True,
    type=str,
)
@click.option(
    "--nuclease",
    "-n",
    help=f"Cas nuclease to design for.  Current options include {AVAILABLE_NUCLEASES}",
    default="SpCas9",
    required=False,
    type=str,
)
@click.option(
    "--reference",
    "-r",
    help="Path to the directory containing the appropriate Bowtie index.",
    default=None,
    type=str,
)
@click.option("--largeindex", help="", default=False, required=False, type=bool)
@click.option(
    "--rule_set",
    help="On-target score rule set to use. Current options include '1', '2', and 'Azimuth'",
    default="Azimuth",
    type=str,
)
@click.option(
    "--ontarget_score_threshold",
    "-on",
    help="Spacers with an on-target score below this will be ignored.",
    default=0,
    required=False,
    type=int,
)
@click.option(
    "--offtarget_score_threshold",
    "-off",
    help="Spacers with an off-target score below this will be ignored.",
    default=0,
    required=False,
    type=int,
)
@click.option(
    "--offtarget_count_threshold",
    help="Spacers with more than this many off-targets will be ignored.",
    default=100,
    required=False,
    type=int,
)
@click.option(
    "--spacers_per_feature",
    help="Number of spacers to find for each feature.",
    default=6,
    type=int,
)
@click.option(
    "--paired",
    help="Should spacers be designed to work as pairs (e.g. for excision)?",
    default=False,
    type=bool,
)
@click.option(
    "--number_upstream_spacers",
    help="If designing paired spacers, number of spacers to design that target upstream of the feature.",
    default=3,
    type=int,
)
@click.option(
    "--number_downstream_spacers",
    help="If designing paired spacers, number of spacers to design that target downstream of the feature.",
    default=3,
    type=int,
)
@click.option(
    "--min_paired_distance",
    help="If designing paired spacers, minimum space required between the up- and downstream spacers.",
    default=0,
    type=int,
)
@click.option(
    "--cores",
    "-c",
    help="Number of processors to use. By default, will use all available.",
    default=None,
    type=int,
)
@click.help_option()
def main(
    input_sequences: str = None,
    outfile: str = None,
    refgenome: str = None,
    restriction_sites: Optional[List[str]] = None,
    largeindex: bool = False,
    on_target_score_threshold: int = 0,
    off_target_score_threshold: int = 0,
    off_target_count_threshold: int = 0,
    nuclease: str = "SpCas9",
    spacers_per_feature: int = 9,
    reject: bool = False,
    paired: bool = False,
    rules: str = "Azimuth",
    number_upstream_spacers: int = 0,
    number_downstream_spacers: int = 0,
    numcores: int = 0,
) -> None:

    targets = pyfaidx.Fasta(input_sequences)

    global NUCLEASES
    nuc = NUCLEASES[NUCLEASES["nuclease"] == nuclease].to_dict(orient="records")[0]

    spacers_df = find_spacers(
        itemlist=targets, nuclease_info=nuc, restriction_sites=restriction_sites
    )

    initialnumber = spacers_df.shape[0]

    spacers_df = on_target_scoring(
        ruleset=rules,
        spacers=spacers_df,
        on_target_score_threshold=on_target_score_threshold,
    )

    if spacers_df.shape[0] == 0:
        print("Sorry, no spacers matching that criteria were found")
        return 0
    else:
        print(
            f"Finished scoring spacers. {spacers_df.shape[0]} of {initialnumber} "
            f"spacers have an on-target score above the cutoff threshold of {on_target_score_threshold}."
            f"\nBeginning Bowtie alignment..."
        )

    spacers_df["hash"] = spacers_df.apply(lambda x: hash(tuple(x)), axis=1)
    offtarget_results_file = off_target_discovery(
        spacers_df=spacers_df,
        cpus=numcores,
        refgenome=refgenome,
        large_index_size=largeindex,
        reject=reject,
    )

    spacers_df = off_target_scoring(
        otrf=offtarget_results_file,
        spacers_df=spacers_df,
        nuclease_info=nuc,
        off_target_score_threshold=off_target_score_threshold,
        off_target_count_threshold=off_target_count_threshold,
    )

    if paired:
        guide_library = assemble_paired_library(
            spacers=spacers_df,
            on_target_score_threshold=on_target_score_threshold,
            off_target_score_threshold=off_target_score_threshold,
            number_upstream_spacers=number_upstream_spacers,
            number_downstream_spacers=number_downstream_spacers,
        )
    else:
        guide_library = assemble_library(
            spacers=spacers_df,
            on_target_score_threshold=on_target_score_threshold,
            off_target_score_threshold=off_target_score_threshold,
            spacers_per_feature=spacers_per_feature,
        )
    guide_library.to_csv(outfile)
    print("Finished.")


if __name__ == "__main__":
    main()
