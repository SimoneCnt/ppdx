#!/usr/bin/env python3

from .pdb import Pdb
from .compute import eval_descriptors, get_descriptors, eval_pkl, compute_core, save_descriptors_json, clean

from . import config
from . import tools
from . import makemodel
from . import scoring

# Manage scripts and data in data/
import os
import pkg_resources
def get_data_fname(fname):
    return pkg_resources.resource_filename(__package__, os.path.join('data', fname))
def link_data(fname, overwrite=True):
    if os.path.lexists(fname):
        os.remove(fname)
    source = pkg_resources.resource_filename(__package__, os.path.join('data', fname))
    if not os.path.isfile(source):
        raise ValueError('File %s does not exist! Cannot link it...' % (source))
    os.symlink(source, fname)

