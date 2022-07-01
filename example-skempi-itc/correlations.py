#!/usr/bin/env python3

import json
import itertools
import numpy as np
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt


def get_experimental_data(ppdbf='ppdb/ppdb.txt', dG=False, dH=False, TdS=False):
    """
        Load experimental data. Only if temperature is 298K.
        Also skip 4JEU complex.
    """
    names = list()
    exps = list()
    with open(ppdbf, 'r') as fp:
        for line in fp:
            if line.startswith('#'): continue
            splt = line.split()
            name, temperature, adG, adH, aTdS = splt[0], splt[3], splt[4], splt[5], splt[7]
            if temperature!='298': continue
            if name.startswith('4JEU'): continue
            names.append(name)
            if dG: exps.append(float(adG))
            if dH: exps.append(float(adH))
            if TdS: exps.append(float(aTdS))
    return names, exps

def get_data_from_json(jsonf, names, descriptors, protocol):
    """
        Read the computed data for the selected complexes, descriptors and 
        protocol.
    """
    with open(jsonf, 'r') as fp:
        alldata = json.load(fp)
    X = np.zeros((len(names), len(descriptors)), dtype=float)
    for i, desc in enumerate(descriptors):
        for j, name in enumerate(names):
            X[j][i] = alldata[name][protocol][desc]
    return X

def plot_crosscorr(mtx, nnamesx, nnamesy, cmap='Spectral', fname=None):
    matplotlib.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(13,13))
    ax = fig.add_subplot(111)
    cax = ax.matshow(mtx, interpolation='nearest', cmap=plt.get_cmap(cmap))#, vmin=-1, vmax=1)
    fig.colorbar(cax)
    ind_array_x = np.arange(mtx.shape[1])
    ind_array_y = np.arange(mtx.shape[0])
    x, y = np.meshgrid(ind_array_x, ind_array_y)
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

    names, dG = get_experimental_data(dG=True)
    print(len(names), len(dG))
    names, dH = get_experimental_data(dH=True)
    print(len(names), len(dH))
    names, TdS = get_experimental_data(TdS=True)
    print(len(names), len(TdS))

    cnames = ['dG', 'dH', 'TdS']
    mtx_rp = np.zeros((len(desc), len(cnames)))
    for i, d in enumerate(desc):
        v = get_data_from_json('descriptors-avg.json', names, [d], 'modeller_fast').transpose()[0]
        if np.all(np.isfinite(v)):

            rp = stats.pearsonr(v, dG)[0]
            mtx_rp[i][0] = rp
            plt.scatter(v, dG, label='rp: %.2f' % (rp))
            plt.title('%s vs dG' % (d))
            plt.legend()
            plt.savefig('img/dg_%s.png' % (d), bbox_inches='tight', dpi=300)
            plt.clf()

            rp = stats.pearsonr(v, dH)[0]
            mtx_rp[i][1] = rp
            plt.scatter(v, dH, label='rp: %.2f' % (rp))
            plt.title('%s vs dH' % (d))
            plt.legend()
            plt.savefig('img/dh_%s.png' % (d), bbox_inches='tight', dpi=300)
            plt.clf()

            rp = stats.pearsonr(v, TdS)[0]
            mtx_rp[i][2] = rp
            plt.scatter(v, TdS, label='rp: %.2f' % (rp))
            plt.title('%s vs TdS' % (d))
            plt.legend()
            plt.savefig('img/tds_%s.png' % (d), bbox_inches='tight', dpi=300)
            plt.clf()

        else:
            mtx_rp[i][0] = float('nan')
            mtx_rp[i][1] = float('nan')
            mtx_rp[i][2] = float('nan')
            print('Skip', d)

    todel = list()
    newdesc = list()
    for i, d  in enumerate(desc):
        if np.any(np.isnan(mtx_rp[i])):
            todel.append(i)
            print(d)
        else:
            newdesc.append(d)
    mtx_rp = np.delete(mtx_rp, todel, axis=0)
    plot_crosscorr(mtx_rp, cnames, newdesc, 'bwr_r', 'correlations.png')

main()
 
