#!/usr/bin/env python3

import ppdg
from timeit import default_timer as timer
import os, sys, random
import modeller, modeller.automodel
import logging
log = logging.getLogger(__name__)

def modeller_veryfast(sequence, template, wrkdir):
    return modeller_generic(sequence, template, wrkdir, "veryfast")

def modeller_fast(sequence, template, wrkdir):
    return modeller_generic(sequence, template, wrkdir, "fast")

def modeller_slow(sequence, template, wrkdir):
    return modeller_generic(sequence, template, wrkdir, "slow")

def modeller_generic(sequence, template, wrkdir, fast):
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
    # Create the sequence file
    with open("sequence.seq", "w") as fp:
        fp.write(">P1;complex_seq\n")
        fp.write("sequence::.:.:.:.::::\n")
        fp.write(sequence + "*\n")
        fp.write(">P1;template.pdb\n")
        fp.write("structure:template.pdb:FIRST:@:LAST:.::::\n")
        fp.write("*\n")
    # Prepare template
    if not os.path.isfile('template.pdb'):
        if isinstance(template, ppdg.Pdb):
            template.write('template.pdb')
        else:
            template_path = os.path.join(basepath, template)
            if os.path.isfile(template_path):
                os.symlink(os.path.join(basepath, template), 'template.pdb')
            else:
                raise ValueError('Impossible to find template file %s' % (template))
    # Run modeller
    with open("modeller.out", "w") as fp:
        _stdout = sys.stdout
        sys.stdout = fp
        env = modeller.environ(rand_seed=-random.randrange(50000))
        #a = modeller.automodel.allhmodel(env, alnfile='sequence.seq', knowns='template.pdb', sequence='complex_seq')
        a = modeller.automodel.automodel(env, alnfile='sequence.seq', knowns='template.pdb', sequence='complex_seq')
        a.starting_model = 0
        a.ending_model = 0
        a.auto_align()
        if fast=='veryfast':
            a.md_level = None
        elif fast=='fast':
            a.md_level = modeller.automodel.refine.fast
        elif fast=='slow':
            a.md_level = modeller.automodel.refine.slow
            a.repeat_optimization = 2
        else:
            log.error('You are expected to specify fast=[veryfast,fast,slow]')
            quit()
        a.make()
        sys.stdout.flush()
        sys.stdout = _stdout
    # Done!
    if os.path.lexists('model.pdb'):
        os.remove('model.pdb')
    os.symlink("complex_seq.B99990000.pdb", "model.pdb")
    os.chdir(basepath)
    log.info("Model %s/model.pdb created!" % (wrkdir))
    time_end = timer()
    return time_end - time_start

