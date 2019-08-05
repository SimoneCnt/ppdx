#!/usr/bin/env python3

import ppdg
import os, json
import logging
log = logging.getLogger(__name__)

class PPComplex():

    def __init__(self, name, template, sequence, nchains):
        """
            Create a PPComplex object. 
             name       the ID to identify it. Any string is valid.
             template   the path of PDB file used as template to build the complex.
             sequence   the amino-acid sequence of the full complex. 
                        Use / to separate chains.
             n_chains   the number of chains in the receptor and ligand as a 
                        tuple (n_chains_receptor, n_chains_ligand)
            Notes:
            - n_chains_receptor + n_chains_ligand must be equal to the number
                of chains in the sequence
             - the order of the chains in the sequence must match the sequence
                in the template.
             - clean the template as much as possible! Keep only the necessary 
                chains and put them in correct order.
             - The receptor is the first set of chains, the ligand the second.
        """

        # Name/ID for the complex
        self.name = name

        # Target sequence of the complex
        self.sequence = sequence

        # The PDB template to use to make the model of the complex
        self.template = template
        if not os.path.isfile(template):
            raise ValueError('Template file %s does not exist!' % (template))

        # A tuple. First is the number of chains in the receptor, second
        #   the number of chains in the ligand.
        if len(nchains)!=2:
            raise ValueError('Cannot understand nchains = %s. Use (nchains_receptor, nchains_ligand).' % (nchains))
        nchains_rec = nchains[0]
        nchains_lig = nchains[1]
        nchains_tot = len(sequence.split('/'))
        if nchains_rec + nchains_lig != nchains_tot:
            raise ValueError('You said the receptor has %d chains and the ligand %d, but the sequence has a total of %d chains.' 
                    % (nchains_rec, nchains_lig, nchains_tot))
        self.nchains = (nchains_rec, nchains_lig)

        # Working directory
        self.wrkdir = os.path.join(ppdg.WRKDIR, self.name.replace('/', '-').replace(' ', '-').lower())

        # Descriptors
        self.descriptors = dict()
        self.descriptors_changed = False
        self.read_descriptors()

        log.info('PPComplex object created successfully! Working directory is %s' % (self.wrkdir))

    def read_descriptors(self):
        """
            Read descriptors from file. File is stored in working directory.
        """
        scorefile = os.path.join(self.wrkdir, "descriptors.json")
        if os.path.isfile(scorefile):
            log.info("Reading descriptors from %s" % (scorefile))
            with open(scorefile, 'r') as fp:
                self.descriptors = json.load(fp)

    def save_descriptors(self):
        """
            Save descriptors to file in working directory.
            Save to disk only if some descriptor changed.
        """
        scorefile = os.path.join(self.wrkdir, "descriptors.json")
        if self.descriptors_changed:
            log.info("Writing descriptors to file %s" % (scorefile))
            with open(scorefile, 'w') as fp:
                json.dump(self.descriptors, fp, indent=4)
            self.descriptors_changed = False

    def make_model(self, protocol, index):
        """
            Build a molecular model of the protein protein complex
            using the specified protocol and the given index.
        """
        wrkdir = os.path.join(self.wrkdir, protocol+'_%d' % index)
        time = ppdg.makemodel.make_model(wrkdir, protocol, self.template, self.sequence)
        if time:
            if not protocol in self.descriptors:
                self.descriptors[protocol] = dict()
            if not '>TIME_makemodel' in self.descriptors[protocol]:
                self.descriptors[protocol]['>TIME_makemodel'] = dict()
            self.descriptors[protocol]['>TIME_makemodel'][str(index)] = time
        ppdg.makemodel.split_model(wrkdir, self.nchains)

    def eval_desc(self, protocol, index, desc_list):
        """
            Compute all descriptors in desc_list for the given protocol/index.
            If a descriptor is already in self.descriptors, do not recompute it.
        """
        desc_wanted = []
        for desc in desc_list:
            try:
                score = self.descriptors[protocol][desc][str(index)]
            except KeyError:
                desc_wanted.append(desc)
        if len(desc_wanted)==0:
            return
        log.info('Computing %s for protocol %s and index %d' % (str(desc_wanted), protocol, index))
        wrkdir = os.path.join(self.wrkdir, protocol+'_%d' % (index))
        scores = ppdg.scoring.evaluate(wrkdir, desc_wanted)
        if not protocol in self.descriptors:
            self.descriptors[protocol] = dict()
        for desc in scores.keys():
            if not desc in self.descriptors[protocol]:
                self.descriptors[protocol][desc] = dict()
        for desc, score in scores.items():
            self.descriptors[protocol][desc][str(index)] = score
        self.descriptors_changed = True

    def eval(self, protocol, desc_list, nmodels):
        """
            Evaluate all descriptors in desclist, using the given protocol
            to build the models. Make an average over nmodels.
        """
        for ndx in range(nmodels):
            self.make_model(protocol, ndx)
            self.eval_desc(protocol, ndx, desc_list)
            self.save_descriptors()
        scores = dict()
        for desc in desc_list:
            lst = np.array([ v for k,v in self.descriptors[protocol][desc].items() ])
            avg = np.mean(lst)
            std = np.std(lst)
            err = std/np.sqrt(len(lst))
            print('%-30s %8.3f %8.3f %8.3f' % (desc, avg, std, err))
            scores[desc] = [avg, std, err]
        return scores

