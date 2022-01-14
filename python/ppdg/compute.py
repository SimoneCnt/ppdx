#!/usr/bin/env python3

import os, json
import joblib
import numpy as np
import ppdg

import parsl
from parsl import python_app

import logging
log = logging.getLogger(__name__)


def eval_descriptors(protocol, desclist, inputs, nmodels=12, config=None):
    """
        Main function!
        Given:
            - The protocol to build the atomistic model of the protein-protein complex
            - A PDB file to use as template for the modelling
            - The sequence of the whole protein-protein complex
            - A tuple indicating how many chains are in the receptor and in the ligand [ nchains = ( nchains_rec, nchains_lig) ]
            - Which descriptors do you want to compute
            and
            - How many models you want to make (to make nice averages)
        This function
            - Makes the models
            - Computes the descriptors
        recycling all possible info reading the descriptors.json file in base_wrkdir (if available).
    """
    new_inputs = list()
    for name, cpxseq, nchains, template in inputs:
        for nmodel in range(nmodels):
            data = [name, protocol, template, cpxseq, nchains, desclist, nmodel, False]
            new_inputs.append(data)
    _run(new_inputs, config)


def get_descriptors(name, protocol, desc_list=None, nmodels=None, median=False):
    """
        Return the average and standard error for the descriptors in 
        desc_list for the first nmodels.
        All descriptors have to be already computed by get_descriptors,
        otherwise an error is thrown.
    """
    wrkdir = os.path.join(ppdg.WRKDIR, name)
    if not desc_list:
        desc_list = ppdg.scoring.all_descriptors()
    descfile = os.path.join(wrkdir, 'descriptors.json')
    if not os.path.isfile(descfile):
        raise ValueError("Missing file descriptors.json in directory %s\nBe sure to have called get_descriptors before!" % (wrkdir))
    log.info("Reading descriptors from %s" % (descfile))
    with open(descfile, 'r') as fp:
        alldesc = json.load(fp)
    if protocol in alldesc:
        descriptors = alldesc[protocol]
    else:
        raise ValueError("Could not find protocol %s in %s\nBe sure to have called get_descriptors before!" % (protocol, wrkdir))
    scores = dict()
    for desc in desc_list:
        if nmodels:
            lst = list()
            for i in range(nmodels):
                lst.append(descriptors[desc][str(i)])
        else:
            lst = np.array([ v for k,v in descriptors[desc].items() ])
        if median:
            avg = np.median(lst)
        else:
            avg = np.mean(lst)
        std = np.std(lst)
        err = std/np.sqrt(len(lst))
        scores[desc] = [avg, std, err]
    return scores


def eval_pkl(pkl, inputs, nmodels=12, config=None):

    unpk = joblib.load(pkl)
    if len(unpk)==4:
        protocol, desclist, scaler, reg = unpk
    elif len(unpk)==6:
        protocol, desclist, scaler, reg, x, y = unpk
        yc = reg.predict(scaler.transform(x))
        rmsd = np.sqrt(np.average(np.square(y-yc)))
        if rmsd>0.001:
            raise ValueError("Checking the pkl failed! RMSE for test set is %lf " % (rmsd))
    else:
        raise ValueError("Error unpacking the pkl file. Should have 4 or 6 elements, but found %d" % (len(unpk)))

    # Compute all descriptors
    new_inputs = list()
    for name, cpxseq, nchains, template in inputs:
        for nmodel in range(nmodels):
            data = [name, protocol, template, cpxseq, nchains, desclist, nmodel, False]
            new_inputs.append(data)
    _run(new_inputs, config)

    # Compute the pkl
    dg = dict()
    for name,_,_,_ in inputs:
        scores = get_descriptors(name, protocol, desc_list=desclist, nmodels=nmodels, median=True)
        desc = np.array([ scores[d][0] for d in desclist ]).reshape(1, -1)
        if scaler:
            desc = scaler.transform(desc)
        dg[name] = reg.predict(desc)[0]

    return dg


def _run(inputs, config):
    if not config:
        config = 'pool'
    if config=='serial':
        return _run_serial(inputs)
    elif config=='parsl':
        return _run_parsl(inputs)
    elif config.startswith('pool'):
        splt = config.split(':')
        if len(splt)==1:
            pool, ncores = splt[0], ""
        else:
            pool, ncores = splt[0], splt[1]
        if pool!='pool':
            raise ValueError('Could not understand string <%s>' % (config))
        if ncores.isdigit():
            ncores = int(ncores)
            if ncores<1: ncores=None
        else:
            ncores = None
        return _run_pool(inputs, ncores)
    else:
        raise ValueError('Could not understand string <%s>' % (config))


def _run_serial(inputs):
    for data in inputs:
        name, protocol, _, _, _, _, nmodel, _ = data
        key = '%s__%s__%s' % (name, protocol, str(nmodel))
        key2, scores = ppdg.compute_core(*data)
        assert (key==key2), "Something very wrong happened :("
        _save(name, protocol, nmodel, scores)


def _run_pool(inputs, ncores):
    from  multiprocessing import Pool
    results = dict()
    with Pool(ncores) as p:
        results = p.starmap(ppdg.compute_core, inputs)
    for data, res in zip(inputs, results):
        name, protocol, _, _, _, _, nmodel, _ = data
        key = '%s__%s__%s' % (name, protocol, str(nmodel))
        key2, scores = res
        assert (key==key2), "Something very wrong happened :("
        _save(name, protocol, nmodel, scores)


def _run_parsl(inputs):
    ppdg_config = ppdg.config.cget()
    futures = dict()
    for i,data in enumerate(inputs):
        name, protocol, _, _, _, _, nmodel, _ = data
        key = '%s__%s__%s' % (name, protocol, str(nmodel))
        futures[key] = _compute_parsl(data, ppdg_config)
        log.info("Task %d key %s submitted" % (i, key))
    for data in inputs:
        name, protocol, _, _, _, _, nmodel, _ = data
        key = '%s__%s__%s' % (name, protocol, str(nmodel))
        key2, scores = futures[key].result()
        assert (key==key2), "Something very wrong happened :("
        _save(name, protocol, nmodel, scores)

@python_app
def _compute_parsl(data, ppdg_config):
    import ppdg
    ppdg.config.cset(ppdg_config)
    return ppdg.compute_core(*data)


def compute_core(name, protocol, template, sequence, nchains, desc_wanted, nmodel, force_calc):

    base_wrkdir = os.path.join(ppdg.WRKDIR, name)
    key = '%s__%s__%s' % (name, protocol, str(nmodel))

    # Read the descriptors from file
    alldesc = dict()
    desc = dict()
    descfile = os.path.join(base_wrkdir, 'descriptors.json')
    if os.path.isfile(descfile):
        #log.info("Reading descriptors from %s" % (descfile))
        with open(descfile, 'r') as fp:
            alldesc = json.load(fp)
        if protocol in alldesc:
            desc = alldesc[protocol]
    desc_flat = ppdg.tools.switch_desc_format(desc)

    # Make a list of what we need to compute
    if str(nmodel) in desc_flat.keys():
        desc_have = desc_flat[str(nmodel)]
    else:
        desc_have = dict()
    desc_set = set()
    for desc in desc_wanted:
        if desc not in desc_have.keys() or force_calc:
            desc_set.add(desc)
    if not desc_set:
        return key, dict()

    log.info('We need to compute these new descriptors: %s for key %s' % (desc_set, key))

    # Check if tgz exists
    wrkdir = os.path.join(base_wrkdir, protocol+'_%s' % (str(nmodel)))
    if os.path.isfile(wrkdir+'.tgz'):
        if os.path.isdir(wrkdir):
            raise ValueError('Both tgz and dir exists for %s' % (wrkdir))
        log.info('Un-compressing %s.tgz' % (wrkdir))
        basepath = os.getcwd()
        os.chdir(base_wrkdir)
        ppdg.tools.execute('tar -xf %s_%s.tgz && rm %s_%s.tgz' % (protocol, str(nmodel), protocol, str(nmodel)))
        os.chdir(basepath)

    # Check template file exists
    if not os.path.isfile(template):
        raise ValueError('Template file %s does not exist! Working directory is %s' % (template, os.getcwd()))

    # Check nchains. A tuple. First is the number of chains in the receptor, 
    # second the number of chains in the ligand.
    if len(nchains)!=2:
        raise ValueError('Cannot understand nchains = %s. Use (nchains_receptor, nchains_ligand).' % (nchains))
    nchains_rec = nchains[0]
    nchains_lig = nchains[1]
    nchains_tot = len(sequence.split('/'))
    if nchains_rec + nchains_lig != nchains_tot:
        raise ValueError('You said the receptor has %d chains and the ligand %d, but the sequence has a total of %d chains.' 
                % (nchains_rec, nchains_lig, nchains_tot))

    # Compute
    scores = ppdg.makemodel.make_model(wrkdir, protocol, template, sequence)
    nchfile = os.path.join(wrkdir, 'nchains.dat')
    if not os.path.isfile(nchfile):
        with open(nchfile,'w') as fp:
            fp.write('%d %d' % nchains)
    if protocol in ['modeller_veryfast']:
        nsteps = 10
    else:
        nsteps = 100
    ppdg.makemodel.charmify(os.path.join(wrkdir, 'model.pdb'), nsteps=nsteps)
    ppdg.makemodel.split_complex(wrkdir, nchains)
    scores2 = ppdg.scoring.evaluate(wrkdir, desc_wanted, desc_have, force_calc)
    scores.update(scores2)
    return key, scores


def _save(name, protocol, nmodel, scores):

    if not scores:
        return

    base_wrkdir = os.path.join(ppdg.WRKDIR, name)

    # Read the descriptors from file
    alldesc = dict()
    desc = dict()
    descfile = os.path.join(base_wrkdir, 'descriptors.json')
    if os.path.isfile(descfile):
        with open(descfile, 'r') as fp:
            alldesc = json.load(fp)
        if protocol in alldesc:
            desc = alldesc[protocol]
    desc_flat = ppdg.tools.switch_desc_format(desc)
    desc_flat[str(nmodel)] = scores

    # Write descriptors to file
    desc2 = ppdg.tools.switch_desc_format(desc_flat)
    if desc2!=desc:
        alldesc[protocol] = desc2
        with open(descfile, 'w') as fp:
            json.dump(alldesc, fp, indent=4, sort_keys=True)
    return


def clean(wrkdir=None):
    # List of files generated during model construction and evaluation of the descriptors.
    # Commented file names should be kept, uncommented ones can be deleted without problems.
    # wrkdir is the same as ppdg.WRKDIR

    if not wrkdir:
        wrkdir = ppdg.WRKDIR

    log.info('Cleaning directory %s' % (wrkdir))

    filelist = [

        ## Modeller
        'complex_seq.D00000000', 'complex_seq.ini', 'complex_seq.rsr', 'complex_seq.sch',
        'complex_seq.V99990000', 'family.mat', 'modeller.out',
        #'complex_seq.B99990000.pdb', 'model.pdb', 'template.pdb', 'sequence.seq', 'sequence.seq.ali',

        ## Charmify
        'add_disulfide.str', 'buildgen.inp',
        # 'chain_a.pdb', 'chain_b.pdb', 'chain_c.pdb',
        'disu.str', 'model-chm.err', 'model-chm.out',
        #'model-chm.cor', 'model-chm.pdb', 'model-chm.psf',
        #'complex-chm.cor', 'complex-chm.pdb', 'complex-chm.psf',

        # Split complex
        'extract.inp', 'ligand-chm.err', 'ligand-chm.out', 'receptor-chm.err', 'receptor-chm.out',
        #'complex.pdb', 'complexAB.pdb', 'ligand.pdb', 'ligandB.pdb', 'receptor.pdb', 'receptorA.pdb',
        #'ligand-chm.psf', 'ligand-chm.cor', 'ligand-chm.pdb', 'receptor-chm.cor', 'receptor-chm.pdb', 'receptor-chm.psf',

        'ligand-chm-facts.pdb', 'ligand-chm-gbsw.pdb', 'ligand-chm-gbmv.pdb',
        'receptor-chm-facts.pdb', 'receptor-chm-gbsw.pdb', 'receptor-chm-gbmv.pdb',
        'complex-chm-facts.pdb', 'complex-chm-gbsw.pdb', 'complex-chm-gbmv.pdb',
    ]

    wrkdir = os.path.abspath(wrkdir)
    basepath = os.getcwd()
    sysdir = [os.path.join(wrkdir, o) for o in os.listdir(wrkdir) if os.path.isdir(os.path.join(wrkdir, o))]

    for sysname in sysdir:
        os.chdir(sysname)
        modelsdirs = [o for o in os.listdir(sysname) if os.path.isdir(o)]
        for model in modelsdirs:
            if os.path.isfile('%s.tgz' % (model)):
                raise ValueError('Both folder and tgz exist for %s/%s' % (sysname, model))
            log.info('Cleaning dir %s/%s' % (sysname, model))
            for f in filelist:
                torm = os.path.join(model, f)
                if os.path.isfile(torm):
                    os.remove(torm)
            #log.info('Compressing dir %s/%s' % (sysname, model))
            ppdg.tools.execute('tar --remove-files -czf %s.tgz %s' % (model, model))
        os.chdir(basepath)

