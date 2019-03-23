#!/usr/bin/env python

from setuptools import setup

setup(name='MerryCRISPR',
      version='2.0',
      author='Miles Smith',
      author_email="mileschristiansmith@gmail.com",
      description="Machine Learning-Based Predictive Modelling of CRISPR/Cas9 guide efficiency",
      url='https://gitlab.com/milothepsychic/merrycrispr',
      license='Proprietary',
      classifiers=['Development Status :: 3 - Production',
                   'Intended Audience :: Developers',
                   'Programming Language :: Python :: 3.6', 'Programming Language :: Python :: 3.7'],
      packages=["merrycrispr"],
      keywords='CRISPR',
      python_requires='>=3.6',
      package_dir={
            'merrycrispr': 'merrycrispr'
      },
      package_data={
            'merrycrispr': ['data/*.*']
      },
      install_requires=['biopython', 'progressbar2', 'regex', 'Azimuth', 'pyfaidx', 'click', 'gtfparse', 'pandas', 'numpy'],
      dependency_links=['https://github.com/lmjohns3/theanets/tarball/master#egg=theanets-0.8.0rc0']
      )