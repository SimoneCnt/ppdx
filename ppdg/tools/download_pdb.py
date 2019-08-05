#!/usr/bin/env python3

import os
import urllib.request
import logging
log = logging.getLogger(__name__)


def download_pdb(pdbid, saveto=None, overwrite=False):
    '''
        Download a PDB file from the RSCB Protein Data Bank
    '''
    url = 'https://files.rcsb.org/download/%s.pdb' % (pdbid)
    if not saveto:
        saveto = '%s.pdb' % (pdbid)
    if os.path.isfile(saveto) and not overwrite:
        return
    log.info('Downloading PDB %s and saving as %s' % (pdbid, saveto))
    response = urllib.request.urlopen(url)
    with open(saveto, 'w') as fp:
        for line in response.readlines():
            line = line.decode('utf-8')
            fp.write(line)

