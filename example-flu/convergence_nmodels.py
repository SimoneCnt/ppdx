#!/usr/bin/env python3

import json
import math
import numpy as np
import matplotlib.pyplot as plt

def main():

    protocol = 'modeller_fast'

    with open('descriptors-all.json', 'r') as fp:
        alldata = json.load(fp)

    desclist = set([ k for name in alldata.keys() for k in alldata[name][protocol].keys() if not k.startswith(">TIME") and k!='NRES' ])

    nn = max([int(v) for v in set([k for name in alldata.keys() for desc in alldata[name][protocol].keys() for k in alldata[name][protocol][desc].keys()])])

    numd = len(desclist)
    ncols = 10
    nrows = math.ceil(numd/ncols)
    tmp = 2.1
    plt.figure(figsize=(tmp*ncols,tmp*nrows-0.5))
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    for n, desc in enumerate(sorted(desclist)):
        plt.subplot(nrows, ncols, n+1)
        if desc=='FoldX_solvation_hydrophobic':
            plt.title('FoldX_solv_hydroph')
        elif desc=='FoldX_solvation_polar':
            plt.title('FoldX_solv_polar')
        elif desc=='FoldX_sidechain_hbond':
            plt.title('FoldX_sidechan_hb')
        elif desc=='FoldX_backbone_hbond':
            plt.title('FoldX_backbone_hb')
        elif desc=='FoldX_entropy_sidechain':
            plt.title('FoldX_sidechan_S')
        elif desc=='FoldX_entropy_mainchain':
            plt.title('FoldX_mainchan_S')
        else:
            plt.title(desc)
        if desc in ['AGBNP', 'GBMV_POL', 'SOAP-Protein-OD']:
            plt.ylim(0.9, 19)
        else:
            plt.ylim(0.9, 1.9)
            plt.yticks([1.0, 1.4, 1.8])
        plt.xlim(1, nn+1)
        plt.xticks([4,8,12,16,20])
        for name, data in alldata.items():
            if desc not in data[protocol].keys(): continue
            if not name.startswith('FY'): continue
            v = data[protocol][desc]
            nmax = max([ int(k) for k in v.keys() ])+1
            vals = np.empty((nmax), dtype=float)
            for k in v.keys():
                vals[int(k)] = v[k]
            avgs = [ np.average(vals[0:n+1]) for n in range(nmax)]
            avgs = np.abs(avgs/avgs[-1] - 1.0) + 1.0
            plt.plot(range(1,nmax+1), avgs, 'k')

    plt.savefig('fig/convergence_nmodels.png', bbox_inches='tight', dpi=300)

main()

