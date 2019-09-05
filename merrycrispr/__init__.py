from merrycrispr.__main__ import main
from merrycrispr._version import get_versions
from merrycrispr.find_spacers import find_spacers
from merrycrispr.library_assembly import assemble_library, assemble_paired_library
from merrycrispr.off_target_scoring import (
    hsu_offtarget_score,
    sumofftargets,
    off_target_discovery,
    off_target_scoring,
)
from merrycrispr.on_target_scoring import on_target_scoring, score_entry
from merrycrispr.rule_set_one import calc_score
from merrycrispr.seqextractor import (
    extract,
    extract_for_tss_adjacent,
    match_seq,
    split_records,
    display_gtf_features,
    display_gtf_genes,
    display_gtf_geneids,
)
from merrycrispr.species_getter import (
    EnsemblRestClient,
    available_species,
    build_bowtie_index,
    get_resources,
)

__author__ = ("Miles Smith",)
__email__ = "mileschristiansmith@gmail.com"

__version__ = get_versions()["version"]

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
