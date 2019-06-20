#!/usr/bin/env python
import sys
if sys.version_info < (3,):
    sys.exit('MerryCRISPR requires Python >= 3.6')
from pathlib import Path
from setuptools import setup, find_packages
import versioneer

try:
    from merrycrispr import __author__, __email__
except ImportError:  # Deps not yet installed
    __author__ = __email__ = ''

setup(
    name="MerryCRISPR",
    version=versioneer.get_version(),
    description="Machine Learning-Based Predictive Modelling of CRISPR/Cas9 guide efficiency",
    long_description=Path('README.rst').read_text('utf-8'),
    url="https://github.com/milescsmith/merrycrispr",
    author=__author__,
    author_email=__email__,
    license="Proprietary",
    python_requires=">=3.6",
      install_requires=[
        l.strip() for l in
        Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    classifiers=[
        "Development Status :: 3 - Production",
        'Environment :: Console',
        "Intended Audience :: Developers",
        'Programming Language :: Python :: 3',
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    extras_require=dict(
        doc=['sphinx', 'sphinx_rtd_theme', 'sphinx_autodoc_typehints'],
    ),
    packages=find_packages(),
    include_package_data=True,

    keywords="CRISPR",
    package_dir={"merrycrispr": "merrycrispr"},
    package_data={"merrycrispr": ["data/*.*"]},
    dependency_links=[
        "https://github.com/lmjohns3/theanets/tarball/master#egg=theanets-0.8.0rc0",
        "https://github.com/milescsmith/Azimuth/tarball/master#egg=azimuth-3.0",
    ],
)
