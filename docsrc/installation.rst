Installation
------------

There are three potential methods one can use to install MerryCRISPR, in order
of easiest to most difficult: `Docker`_, `Miniconda`_, and `PyPI`_.

Docker
~~~~~~~
`Docker <https://www.docker.com/>`_ is a tool for running software as
self-contained images, packinging the necessary operating system and libraries
together. This means that even though MerryCRISPR was written on MacOS and
works most reliably when run on an Unix-style OS, the MerryCRISPR container
should work well when run on a Windows 10-based system. Also, no other software
needs to be installed once the MerryCRISPR container has been. (A form of
Docker *is* available for Win7, but it is more difficult to setup and use due
to some necessary features that were not introduced until Win10).

`Installation <https://docs.docker.com/install/>`_ of Docker is beyond the
scope of this guide, but should be relatively straightforward.  Once Docker is
installed and working, the MerryCRISPR container can be `pulled
<https://gitlab.com/milothepsychic/merrycrispr/container_registry>`_ from the
project's GitLab repository::

    docker login registry.gitlab.com
    docker pull registry.gitlab.com/milothepsychic/merrycrispr

(note that currently, you will need both a GitLab account and permissions to
the project in order to pull the image).

Alternatively, the docker image can be built locally and installed by runnning
the following command from within the base merrycrispr directory::

    docker build . -t merrycrispr:latest

You can check that the container has been successfully installed and is working
by running::

    docker run merrycrispr merrycrispr --help

For use, see the :ref:`usage:Docker` section in :ref:`usage:Usage`.


Miniconda
~~~~~~~~~~~~~
`Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ is a package
manager, initially focusing on Python packages by now covering much more,
especially bioinfomatics software.  Conda makes it easy to install sofware in
descrete "environments", where all software is automatically installed with all
its dependencies, where only compatible versions of software and libraries are
installed, and everthing is kept seperate from outside packages or those in
other environments.  If Python is not already installed on your, Conda is the
best option to use.
After downloading and installing Conda, create an en

PyPI
~~~~~
More things.

