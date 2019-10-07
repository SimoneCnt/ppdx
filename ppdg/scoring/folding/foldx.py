
import os
from timeit import default_timer as timer
import ppdg
import logging
log = logging.getLogger(__name__)

def foldx(wrkdir, cpxname='complexAB.pdb'):
    """
        FoldX interaction energy [1].
        First run the RepairPDB utility (which takes quite a lot of time...)
        and then calculate the binding free energy using the AnalyseComplex
        utility. Return the total binding free energy as well as some of its
        components.
        [1] J. Schymkowitz, J. Borg, F. Stricher, R. Nys, F. Rousseau, and 
            L. Serrano, "The FoldX web server: an online force field", Nucleic 
            Acids Research, vol. 33, pp. W382-W388, 2005.
    """
    time_start = timer()
    basepath = os.getcwd()
    os.chdir(wrkdir)
    log.info("Getting FoldX scoring...")
    if not os.path.isfile("rotabase.txt"):
        if os.path.isfile(os.path.join(ppdg.FOLDX, 'rotabase.txt')):
            os.symlink(os.path.join(ppdg.FOLDX, 'rotabase.txt'), 'rotabase.txt')
        else:
            raise ValueError('Impossible to find file rotabase.txt in directory %s. Fix FOLDX variable in config.' % (ppdg.FOLDX))
    ret = ppdg.tools.execute(os.path.join(ppdg.FOLDX,'foldx') + " --command=RepairPDB --pdb=%s" % (cpxname) + " >foldx_repair.err 2>&1")
    repairname = cpxname[:-4]+'_Repair.pdb'
    ret = ppdg.tools.execute(os.path.join(ppdg.FOLDX,'foldx') + " --command=AnalyseComplex --pdb=%s --analyseComplexChains=A,B --output-file=foldx" % (repairname) + " >foldx_analyse.err 2>&1")
    desc = dict()
    with open('Interaction_foldx_AC.fxout','r') as fp:
        for line in fp.readlines():
            lsplit = line.split()
            if len(lsplit)==32:
                desc['FoldX'] = float(lsplit[5])
                desc['FoldX_backbone_hbond'] = float(lsplit[6])
                desc['FoldX_sidechain_hbond'] = float(lsplit[7])
                desc['FoldX_vdw'] = float(lsplit[8])
                desc['FoldX_elec'] = float(lsplit[9])
                desc['FoldX_solvation_polar'] = float(lsplit[10])
                desc['FoldX_solvation_hydrophobic'] = float(lsplit[11])
                desc['FoldX_entropy_sidechain'] = float(lsplit[13])
                desc['FoldX_entropy_mainchain'] = float(lsplit[14])
    os.chdir(basepath)
    time_end = timer()
    desc['>TIME_FoldX'] = time_end-time_start
    return desc
     
if __name__=='__main__':
    import sys, os
    ppdg.readconfig('ppdg.ini')
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    if len(sys.argv)==2:
        desc = foldx(sys.argv[1])
    elif len(sys.argv)==5:
        desc = foldx(sys.argv[1], sys.argv[2])
    else:
        print('Usage: %s wrkdir [complexAB.pdb]' % (sys.argv[0]))
        quit()
    print(desc)

