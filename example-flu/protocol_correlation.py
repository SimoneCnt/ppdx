#!/usr/bin/env python3

import json
import itertools
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt

def get_dset2(alldata, p1, p2, descriptor, complexes):
    pairs = list()
    for cpx in complexes:
        v1 = None
        if p1 in alldata[cpx].keys():
            if descriptor in alldata[cpx][p1].keys():
                v1 = alldata[cpx][p1][descriptor]
        v2 = None
        if p2 in alldata[cpx].keys():
            if descriptor in alldata[cpx][p2].keys():
                v2 = alldata[cpx][p2][descriptor]
        if v1 and v2:
            pairs.append((v1,v2))
    return np.array(pairs, dtype=float).transpose()

def plot_crosscorr(mtx, nnamesx, nnamesy, cmap='Spectral', fname=None):
    matplotlib.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(13,13))
    ax = fig.add_subplot(111)
    cax = ax.matshow(mtx, interpolation='nearest', cmap=plt.get_cmap(cmap))#, vmin=-1, vmax=1)
    fig.colorbar(cax)
    ind_array_x = np.arange(mtx.shape[1])
    ind_array_y = np.arange(mtx.shape[0])
    x, y = np.meshgrid(ind_array_x, ind_array_y)
    #for x_val, y_val in zip(x.flatten(), y.flatten()):
    #    c = str(np.around(mtx[x_val][y_val], 2))
    #    ax.text(x_val, y_val, c, va='center', ha='center', fontsize=4)
    ax.set_xticklabels(nnamesx, rotation=90)
    ax.set_yticklabels(nnamesy)
    ax.set_xticks(ind_array_x)
    ax.set_yticks(ind_array_y)
    if fname:
        plt.savefig(fname, bbox_inches='tight', dpi=300)
    else:
        plt.show()
    plt.clf()


def main():

    with open('descriptors-avg.json', 'r') as fp:
        alldata = json.load(fp)

    cpxs = list((alldata.keys()))
    prot = set()
    desc = set()

    for cpx, protocols in alldata.items():
        for protocol, data in protocols.items():
            prot.add(protocol)
            for d, v in data.items():
                if not d.startswith('>TIME_') and d!='NRES':
                    desc.add(d)
    desc = sorted(desc)

    comb = [('modeller_fast', 'modeller_veryfast'), ('modeller_fast', 'modeller_slow'), ('modeller_fast', 'rosetta'),
            ('modeller_veryfast', 'modeller_slow'), ('modeller_veryfast', 'rosetta'), ('modeller_slow', 'rosetta')]
    cnames = ['MF_MVF', 'MF_MS', 'MF_R', 'MVF_MS', 'MVF_R', 'MS_R']

    mtx_rp = np.zeros((len(desc), len(comb)))
    mtx_slope = np.zeros((len(desc), len(comb)))
    for i, d in enumerate(desc):
        for j, c in enumerate(comb):
            v = get_dset2(alldata, c[0], c[1], d, cpxs)
            if len(v)>0:
                rp = stats.pearsonr(v[0], v[1])[0]
                A = np.vstack([v[0], np.ones(len(v[0]))]).T
                m, q = np.linalg.lstsq(A, v[1])[0]
            else:
                rp = float('nan')
                m = float('nan')
            mtx_rp[i][j] = rp
            mtx_slope[i][j] = abs(1-m)
    todel = list()
    newdesc = list()
    for i, d  in enumerate(desc):
        if np.all(np.isnan(mtx_rp[i])):
            todel.append(i)
        else:
            newdesc.append(d)
    mtx_rp = np.delete(mtx_rp, todel, axis=0)
    mtx_slope = np.delete(mtx_slope, todel, axis=0)
    plot_crosscorr(mtx_rp, cnames, newdesc, 'bwr_r', 'fig/protocol_correlation_rp.png')
    plot_crosscorr(mtx_slope, cnames, newdesc, 'bwr', 'fig/protocol_correlation_slope.png')

main()
 
