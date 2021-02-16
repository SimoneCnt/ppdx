#!/usr/bin/env python3

from .config import readconfig, printconfig
from .pdb import Pdb
from .get_descriptors import get_descriptors, get_descriptors_average, clean, eval_pkl
from .mmgbsa import mmgbsa
from .runomm import runomm

from . import tools
from . import makemodel
from . import scoring


# Manage scripts and data in data/
import os
import pkg_resources
def get_data_fname(fname):
    return pkg_resources.resource_filename(__package__, os.path.join('data', fname))
def link_data(fname, overwrite=False):
    if os.path.isfile(fname):
        if overwrite:
            os.remove(fname)
        else:
            return
    source = pkg_resources.resource_filename(__package__, os.path.join('data', fname))
    if not os.path.isfile(source):
        raise ValueError('File %s does not exist! Cannot link it...' % (source))
    os.symlink(source, fname)

