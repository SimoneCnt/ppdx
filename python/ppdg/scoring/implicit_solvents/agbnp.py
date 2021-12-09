#!/usr/bin/env python3

import os, sys
from timeit import default_timer as timer
import logging
log = logging.getLogger(__name__)

try:
    import modeller
    HAVE_MODELLER = True
except ImportError:
    HAVE_MODELLER = False


if HAVE_MODELLER:
    class AGBNPScorer(modeller.terms.AssessEnergyTerm):
        """Score the model using GB/SA implicit solvation."""
        name = 'AGBNP Implicit Solvent'

        def __init__(self, library='${LIB}/solv.lib', solvation_model=1,
                 cutoff=12.0, group=modeller.physical.gbsa):
            modeller.terms.energy_term.__init__(self)
            self.__library = library
            self.__solvation_model = solvation_model
            self.__cutoff = cutoff
            self._group = group

        def _add_term(self, edat, indx):
            import _modeller
            _modeller.mod_gbsa_create(edat, indx, self._group.get_type(),
                                  self.__library, self.__solvation_model,
                                  self.__cutoff)

def agbnp(wrkdir):
    """
        Get AGBNP solvation term for protein-protein as implemented in 
        modeller.
        [1] E. Gallicchio and R. M. Levy, "AGBNP: An analytic implicit solvent 
            model suitable for molecular dynamics simulations and high-resolution 
            modeling", Journal of Computational Chemistry, vol. 25, no. 4, 
            pp. 479-499, 2004.
    """
    import modeller
    time_start = timer()
    log.info("Getting AGBNP scoring...")
    cpx = os.path.join(wrkdir, 'complexAB.pdb')

    with open(os.path.join(wrkdir, "agbnp.out"), "w") as fp:
        _stdout = sys.stdout
        sys.stdout = fp
        env = modeller.environ()
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        mdl = modeller.scripts.complete_pdb(env, cpx)
        cpx = modeller.selection(mdl.chains[:])
        lig = modeller.selection(mdl.chains[0])
        rec = modeller.selection(mdl.chains[1])
        agbnp = AGBNPScorer()
        score = cpx.assess(agbnp) - rec.assess(agbnp) - lig.assess(agbnp)
        sys.stdout.flush()
        sys.stdout = _stdout

    desc = dict()
    desc['AGBNP'] = score
    time_end = timer()
    desc['>TIME_AGBNP'] = time_end - time_start
    return desc


