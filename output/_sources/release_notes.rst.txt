.. role:: small
.. role:: smaller
.. role:: noteversion


Version 2.0 :small:`June 17, 2019`
------------------------------------
- Initial release of major overhaul of MerryCRISPR:
   - Replaced custom GTF/GFF parsing functions with `gtfparse <https://github.com/openvax/gtfparse>`_
   - Handling of FASTA files by `pyfaidx <https://github.com/mdshw5/pyfaidx>`_
   - Removed hardcoding of nuclease information to a csv file, so new nucleases can be easily added.
   - Almost total refactoring of code.