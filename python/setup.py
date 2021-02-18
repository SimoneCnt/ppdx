#!/usr/bin/env python

from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name = 'ppdg',
    version = '0.2',
    description = 'Build protein-protein complexes and calculate binding affinity',
    long_description = readme(),
    url = 'https://github.com/SimoneCnt/ppdg',
    author = 'Simone Conti',
    author_email = 'simonecnt@gmail.com',
    license = 'All rights reserved',
    packages = find_packages(),
    package_data = {'ppdg': ['data/*']},
    #install_requires=[],
)

