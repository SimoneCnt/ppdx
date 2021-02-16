#!/usr/bin/env python3

import ppdg
from timeit import default_timer as timer
import os, sys, random
import Bio.pairwise2
import Bio.SubsMat.MatrixInfo
import logging
log = logging.getLogger(__name__)

def rosetta1(target_sequence, template, wrkdir):
    """
        Build one model for the complex. Save it in the wrkdir folder.
    """

    # Make working directory and enter it
    if os.path.isdir(wrkdir):
        if os.path.isfile(os.path.join(wrkdir, 'model.pdb')):
            return
    else:
        os.makedirs(wrkdir)
    time_start = timer()
    log.info("Creating model %s/model.pdb..." % (wrkdir))
    basepath = os.getcwd()
    os.chdir(wrkdir)

    # Create fasta file with target sequence
    target_sequence_file = 'target_sequence.fasta'
    with open(target_sequence_file, 'w') as fp:
        fp.write(">target\n%s\n" % (target_sequence))

    # Prepare template
    if not os.path.isfile('PDBID.pdb'):
        if isinstance(template, ppdg.Pdb):
            template.write('PDBID.pdb')
        else:
            template_path = os.path.join(basepath, template)
            if os.path.isfile(template_path):
                os.symlink(os.path.join(basepath, template), 'PDBID.pdb')
            else:
                raise ValueError('Impossible to find template file %s' % (template))

    # Remove / because Bio do not recognize them and rosetta thread neither...
    target_sequence = target_sequence.replace('/','')
    template_sequence = ppdg.Pdb("PDBID.pdb").get_sequence().replace('/','')

    # Make alignment
    alignment = Bio.pairwise2.align.globalds(target_sequence, template_sequence, Bio.SubsMat.MatrixInfo.blosum62, -10, -0.5)
    s1 = alignment[0][0]
    s2 = alignment[0][1]

    # Write sequence alignemnt in grishnan (?) format
    with open('alignment.gri', 'w') as fp:
        fp.write("## target PDBID_thread\n")
        fp.write("# \n")
        fp.write("scores_from_program: 0\n");
        fp.write("0 " + s1 + "\n")
        fp.write("0 " + s2 + "\n")
        fp.write("--\n")

    log.info('Running rosetta thread...')
    rthread = ppdg.ROSETTABIN+"/partial_thread.static.linuxgccrelease " + \
                " -database " + ppdg.ROSETTA + "/main/database" + \
                " -in:file:fasta " + target_sequence_file + \
                " -in:file:alignment alignment.gri" + \
                " -in:file:template_pdb PDBID.pdb" + \
                " -ignore_unrecognized_res " + \
                " >rosetta_thread.out 2>&1 "
    ppdg.tools.execute(rthread);
    if not os.path.isfile('PDBID_thread.pdb'):
        log.error('Rosetta thread failed! Look at the output file rosetta_thread.out')

    log.info("Running Rosetta script...")
    ppdg.link_data('hybridize.xml')
    rscript = ppdg.ROSETTABIN + "/rosetta_scripts.static.linuxgccrelease " + \
                " -database " + ppdg.ROSETTA + "/main/database " + \
                " -in:file:fasta " + target_sequence_file + \
                " -parser:protocol hybridize.xml " + \
                " -default_max_cycles 200 " + \
                " -dualspace " + \
                " -restore_talaris_behavior " + \
                " -score:set_weights pro_close 0 " + \
                " >rosetta_script.out 2>&1 "
    ppdg.tools.execute(rscript);
    if not os.path.isfile('S_0001.pdb'):
        log.error('Rosetta hybridize failed! Look at the output file rosetta_script.out')

    # Done!
    if os.path.lexists('model.pdb'):
        os.remove('model.pdb')
    os.symlink("S_0001.pdb", "model.pdb")
    os.chdir(basepath)
    log.info("Model %s/model.pdb created!" % (wrkdir))
    time_end = timer()
    return time_end - time_start

