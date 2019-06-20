from .merrycrispr import main
from .find_spacers import find_spacers
from .on_target_scoring import on_target_scoring, score_entry
from .off_target_scoring import scoreCas9offtarget, sumofftargets, off_target_discovery, off_target_scoring
from .library_assembly import assemble_library, assemble_paired_library
from .seqextractor import extract, extract_for_tss_adjacent, match_seq, split_records, display_gtf_features, display_gtf_genes, display_gtf_geneids

__doc__ = """\
merrycrispr
------------

Main function used to create library:

.. autosummary::
   :toctree: .

   main

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

Description

.. autosummary::
   :toctree: .

   scoreCas9offtarget
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

__author__ = 'Miles Smith',
__email__ = "mileschristiansmith@gmail.com"

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions