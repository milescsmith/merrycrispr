#!/usr/bin/env python
import sys

if sys.version_info < (3,):
    sys.exit("MerryCRISPR requires Python >= 3.6")
from pathlib import Path
from setuptools import setup, find_packages
import versioneer

try:
    from merrycrispr import __author__, __email__
except ImportError:  # Deps not yet installed
    __author__ = __email__ = ""

setup(
    name="MerryCRISPR",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Generate gene and genome-wide CRISPR/Cas9 knockout and excision libraries",
    long_description=Path("README.rst").read_text("utf-8"),
    url="https://github.com/milescsmith/merrycrispr",
    author=__author__,
    author_email=__email__,
    license="LGPL-3.0-or-later",
    python_requires=">=3.6",
    install_requires=[
        l.strip() for l in Path("requirements.txt").read_text("utf-8").splitlines()
    ],
    classifiers=[
        "Development Status :: 3 - Production",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    extras_require=dict(doc=["sphinx", "sphinx_rtd_theme", "sphinx_autodoc_typehints"]),
    packages=find_packages(),
    include_package_data=True,
    entry_points={"console_scripts": ["merrycrispr = merrycrispr.__main__:main"]},
    keywords="CRISPR",
    package_dir={"merrycrispr": "merrycrispr"},
    package_data={"merrycrispr": ["data/*.*"]},
)
