Installation
------------

There are three potential methods one can use to install MerryCRISPR, in order of easiest to most difficult: `Docker`_, `Miniconda`_, and `PyPI`_.

Docker
~~~~~~~
`Docker <https://www.docker.com/>`_ is a tool for running software as self-contained images, packinging the necessary operating system and libraries together.
This means that even though MerryCRISPR was written on MacOS and works most reliably when run on an Unix-style OS, the MerryCRISPR container should work well when run on a Windows 10-based system.
Also, no other software needs to be installed once the MerryCRISPR container has been.
(A form of Docker *is* available for Win7, but it is more difficult to setup and use due to some necessary features that were not introduced until Win10).

Installation of Docker is beyond the scope of this guide.  Once it is installed and working, the MerryCRISPR container and be installed from within the merrycrispr directory by running::

    docker build . -t merrycrispr:latest

You can check that the container has been successfully installed and is working by running::

    docker run merrycrispr merrycrispr --help

    ├── Dockerfile
    ├── LICENSE
    ├── README.md
    ├── README.rst
    ├── docker_excision_test.fa
    ├── docs
    ├── ensembl
    ├── environment.yml
    ├── merrycrispr
    ├── output
    ├── requirements.txt
    ├── setup.cfg
    ├── setup.py
    └── versioneer.py


Miniconda
~~~~~~~~~~~~~
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ is a package manager, initially focusing on Python packages by now encompasing
much more, especially bioinfomatics software.  Conda makes installing sofware in descrete "environments" - where all software is automatically
installed with all its dependencies, where only compatible versions of software and libraries are installed, and everthing is kept seperate from
outside packages or those in other environments.  If Python is not already installed on your, Conda is the best option to use.

PyPI
~~~~~
More things.

