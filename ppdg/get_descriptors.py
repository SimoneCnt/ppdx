#!/usr/bin/env python3

import os, json
import numpy as np
import ppdg
import logging
log = logging.getLogger(__name__)

def get_descriptors_average(wrkdir, protocol, desc_list=None, nmodels=None):
    """
        Return the average and standard error for the descriptors in 
        desc_list for the first nmodels.
        All descriptors have to be already computed by get_descriptors,
        otherwise an error is thrown.
    """
    if not desc_list:
        desc_list = ppdg.scoring.all_descriptors()
    descfile = os.path.join(wrkdir, 'descriptors.json')
    if not os.path.isfile(descfile):
        raise "Missing file descriptors.json in directory %s\nBe sure to have call get_descriptors before!" % (wrkdir)
    log.info("Reading descriptors from %s" % (descfile))
    with open(descfile, 'r') as fp:
        alldesc = json.load(fp)
    if protocol in alldesc:
        descriptors = alldesc[protocol]
    else:
        raise "Could not find protocol %s in %s\nBe sure to have call get_descriptors before!" % (protocol, wrkdir)
    scores = dict()
    for desc in desc_list:
        if nmodels:
            lst = list()
            for i in range(nmodels):
                lst.append(descriptors[desc][str(i)])
        else:
            lst = np.array([ v for k,v in descriptors[desc].items() ])
        avg = np.mean(lst)
        std = np.std(lst)
        err = std/np.sqrt(len(lst))
        scores[desc] = [avg, std, err]
    return scores


def get_descriptors(base_wrkdir, protocol, template, sequence, nchains, desc_wanted, nmodels):
    """
        Main function!
        Given:
            - The working directory to use
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

    # Check template file exists
    if not os.path.isfile(template):
        raise ValueError('Template file %s does not exist!' % (template))

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

    # Read the descriptors from file
    alldesc = dict()
    desc = dict()
    descfile = os.path.join(base_wrkdir, 'descriptors.json')
    if os.path.isfile(descfile):
        log.info("Reading descriptors from %s" % (descfile))
        with open(descfile, 'r') as fp:
            alldesc = json.load(fp)
        if protocol in alldesc:
            desc = alldesc[protocol]
    desc_flat = _switch_desc_format(desc)

    # Compute what needed
    for i in range(nmodels):
        wrkdir = os.path.join(base_wrkdir, protocol+'_%d' % (i))
        if str(i) in desc_flat.keys():
            desc_have = desc_flat[str(i)]
        else:
            desc_have = dict()
        scores = _get_descriptors_core(wrkdir, protocol, template, sequence, nchains, desc_have, desc_wanted)
        desc_flat[str(i)] = scores

    # Write descriptors to file
    desc2 = _switch_desc_format(desc_flat)
    if desc2!=desc:
        alldesc[protocol] = desc2
        log.info("Writing descriptors to %s" % (descfile))
        with open(descfile, 'w') as fp:
            json.dump(alldesc, fp, indent=4, sort_keys=True)
    return alldesc


def _get_descriptors_core(wrkdir, protocol, template, sequence, nchains, desc, desc_wanted):
    """
        Core function to make one model and get the descriptors.
        Do not call this directly, but use get_descriptors.
    """
    scores = ppdg.makemodel.make_model(wrkdir, protocol, template, sequence)
    ppdg.makemodel.charmify(os.path.join(wrkdir, 'model.pdb'))
    ppdg.makemodel.split_complex(wrkdir, nchains)
    scores2 = ppdg.scoring.evaluate(wrkdir, desc_wanted, desc)
    scores.update(scores2)
    return scores


def _switch_desc_format(descriptors):
    """
        Convert dict from 
            index - descriptor - value
        to
            descriptor - index - value
        and vice versa.
        For internal use only...
    """
    desc_flat = dict()
    for key, values in descriptors.items():
        for index, value in values.items():
            if not index in desc_flat.keys():
                desc_flat[index] = dict()
            desc_flat[index][key] = float(value)
    return desc_flat


