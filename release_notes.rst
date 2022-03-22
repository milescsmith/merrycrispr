.. role:: small
.. role:: smaller
.. role:: noteversion


Version 4.1: :small:`March, 21 2022`
-----------------------------------
- Replace `click` with `typer` for cli

Version 3.0: :small:`July,18 2021`
------------------------------------
- Use Poetry for PEP517 compliance
- Add use of pre-commit tools, linters, type checkers, etc....

Version 2.2: :small:`July, 21 2019`
------------------------------------
- Major overhaul of off-target scoring
    - Replaced `for`-loops with Pandas `apply`
    - ~40-fold speed-up


Version 2.1: :small:`July, 5 2019`
------------------------------------
- Introduced `new-species` subcommand
    - Able to fetch annotations and sequences from Ensembl
    - Can build a Bowtie index from fetched sequences
- Added Documentation


Version 2.0: :small:`June, 17 2019`
------------------------------------
- Initial release of major overhaul of MerryCRISPR:
    - Replaced custom GTF/GFF parsing functions with `gtfparse <https://github.com/openvax/gtfparse>`_
    - Handling of FASTA files by `pyfaidx <https://github.com/mdshw5/pyfaidx>`_
    - Removed hardcoding of nuclease information to a csv file, so new nucleases can be easily added.
    - Almost total refactoring of code.
        - Replaced use of lists with Pandas
    - Replaced old Rule Set 2 with a version of Azimuth overhauled for Python 3 compatibility