#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name = 'ppdx',
    version = '1.0',
    description = 'Build protein-protein complexes and compute binding descriptors',
    long_description = 'Build protein-protein complexes and compute binding descriptors',
    url = 'https://github.com/SimoneCnt/ppdx',
    author = 'Simone Conti',
    author_email = 'simonecnt@gmail.com',
    license = 'GPLv3',
    packages = find_packages(),
    package_data = {'ppdg': ['data/*']},
    install_requires = ['numpy', 'matplotlib', 'biopython', 'joblib', 'parsl', 'mdtraj', 'sklearn']
)

