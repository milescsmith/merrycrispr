# src/merrycrispr/__init__.py
from importlib.metadata import metadata, version

try:
    __author__ = metadata(__name__)["Author"]
except KeyError:
    __author__ = "unknown"

try:
    __email__ = metadata(__name__)["Author-email"]
except KeyError:  # pragma: no cover
    __email__ = "unknown"

try:
    __version__ = version(__name__)
except KeyError:  # pragma: no cover
    __version__ = "unknown"

from .find_spacers import find_spacers
from .library_assembly import assemble_library, assemble_paired_library
from .off_target_scoring import (hsu_offtarget_score, off_target_discovery,
                                 off_target_scoring, sumofftargets)
from .on_target_scoring import on_target_scoring, score_entry
from .rule_set_one import calc_score
from .seqextractor import (display_gtf_features, display_gtf_geneids,
                           display_gtf_genes, extract,
                           extract_for_tss_adjacent, match_seq, split_records)
from .species_getter import (EnsemblRestClient, available_species,
                             build_bowtie_index, get_resources)


import typer
from rich.console import Console

console = Console()


def version_callback(value: bool):
    """Prints the version of the package."""
    if value:
        console.print(
            f"[yellow]plinkliftover[/] version: [bold blue]{__version__}[/]"
        )
        raise typer.Exit()


merrycrispr = typer.Typer(
    name="merrycrispr",
    help="Generate CRISPR libraries suitable for screening assays",
    add_completion=False,
)

__doc__ = """\
merrycrispr
------------

Main function used to create library:

.. autosummary::
   :toctree: .

   main

species_getter
---------------

Retrieve annotation and sequence files from Ensembl

.. autosummary::
   :toctree: .

   EnsemblRestClient
   available_species
   get_resources

Prepare a new species for use in MerryCRISPR

.. autosummary::
   :toctree: .

   build_bowtie_index


seqextractor
-------------

Scan GTF for information.

.. autosummary::
   :toctree: .

   display_gtf_features
   display_gtf_genes
   display_gtf_geneids

Create FASTAs to search for spacers.

.. autosummary::
   :toctree: .

   extract
   extract_for_tss_adjacent

Utilities for parsing FASTAs or matching annotations

.. autosummary::
   :toctree: .

   match_seq
   split_records

find_spacers
------------

Description

.. autosummary::
   :toctree: .

   find_spacers


on_target_scoring
-----------------

Description

.. autosummary::
   :toctree: .

   on_target_scoring
   score_entry

off_target_scoring
------------------

Find and score potential off-targets

.. autosummary::
   :toctree: .

   hsu_offtarget_score
   sumofftargets
   off_target_discovery
   off_target_scoring


library_assembly
----------------

Description

.. autosummary::
   :toctree: .

   assemble_library
   assemble_paired_library
"""
