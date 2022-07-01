#!/usr/bin/env python3

import math
import json
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

def expdata(antigen):
    names = list()
    dgs = list()
    with open('ppdb/ppdb.txt') as fp:
        for line in fp:
            if line[0]=='#': continue
            name, _, _, dg, _ = line.split()
            if (not antigen) or (antigen and name.startswith(antigen)):
                names.append(name)
                dgs.append(float(dg))
    return names, dgs
    
def sample(x, std=0.45):
    return np.random.normal(x, std)

def expselfcorr():
    _, g = expdata(None)
    rpl = [ scipy.stats.pearsonr(sample(g), sample(g))[0] for i in range(100000) ]
    plt.hist(rpl, bins='auto', label='rp: %.2f $\pm$ %.2f' % (np.average(rpl), np.std(rpl)))
    plt.xlabel('Pearson correlation coefficient')
    plt.ylabel('Count')
    plt.legend()
    plt.savefig('fig/experimental-self-correlation.png', bbox_inches='tight', dpi=300)

def pearson_convergence(descname, calcd, antigen, protocol, color, shift):
    cpxs, dgexp = expdata(antigen)
    calcv = list()
    for cpx in cpxs:
        calcv.append(list(calcd[cpx][protocol][descname].values()))
    calcv = np.array(calcv)
    calcv = calcv.transpose()
    nmax = calcv.shape[0]
    replica = 1000
    rpconv = list()
    for nsamples in range(1, nmax+1):
        rp = [ scipy.stats.pearsonr(np.average(calcv[idx], axis=0), sample(dgexp, std=0.0))[0] for idx in np.random.choice(nmax, size=(replica, nsamples), replace=True) ]
        rpconv.append((np.average(rp), np.std(rp)))
    rpconv = np.array(rpconv).transpose()
    plt.errorbar(np.array(range(1, nmax+1))+shift, rpconv[0], yerr=rpconv[1], color=color)

def best_correlation(descname, calcd, antigen, protocol, color):
    cpxs, dgexp = expdata(antigen)
    calcv = list()
    for cpx in cpxs:
        calcv.append(list(calcd[cpx][protocol][descname].values()))
    calcv = np.array(calcv)
    calcv = calcv.transpose()
    calc = np.average(calcv, axis=0)
    if antigen:
        plt.scatter(dgexp, calc, color=color, s=16)
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(dgexp, calc)
    #lins = np.linspace(np.min(dgexp), np.max(dgexp), 100)
    lins = np.linspace(-12, -6, 100)
    plt.plot(lins, slope*lins+intercept, '-', color=color, label='%.2f' % (r_value))

def multiplot(mode, protocol):
    desclist = list()
    desclist += ['ZRANK', 'ZRANK2', 'pyDock', 'ATTRACT', 'FireDock']
    desclist += ['RF_HA_SRS', 'RF_CB_SRS_OD', 'ipot_aace167', 'ipot_aace18', 'ipot_aace20', 'ipot_rrce20']
    desclist += ['FACTS_TOT', 'GBSW_TOT']
    desclist += ['FoldX']
    desclist += ['Rosetta_dg', 'Prodigy_IC_NIS']
    desclist += ['DOPE', 'DOPE-HR']

    with open('descriptors-all.json') as fp:
        calcd = json.load(fp)

    numd = len(desclist)
    ncols = 6
    nrows = math.ceil(numd/ncols)
    tmp = 2.1
    plt.figure(figsize=(tmp*ncols,tmp*nrows-0.5))
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    for n, desc in enumerate(sorted(desclist)):
        plt.subplot(nrows, ncols, n+1)
        plt.title(desc)
        if mode=='pearson':
            #plt.xticks([4,8,12,16,20])
            plt.ylim(-0.6,0.8)
            pearson_convergence(desc, calcd, '93th', protocol, 'blue', -0.1)
            pearson_convergence(desc, calcd, 'rsc3', protocol, 'red', +0.1)
            pearson_convergence(desc, calcd, None, protocol, 'black', 0.0)
        else:
            plt.yticks([])
            best_correlation(desc, calcd, '93th', protocol, 'blue')
            best_correlation(desc, calcd, 'rsc3', protocol, 'red')
            best_correlation(desc, calcd, None, protocol, 'black')
            plt.legend(fontsize='x-small', loc='best', ncol=1, columnspacing=1, handlelength=1)

    plt.savefig('fig/%s_%s.png' % (mode, protocol), bbox_inches='tight', dpi=300)


expselfcorr()
multiplot('pearson', 'modeller_veryfast')
multiplot('pearson', 'modeller_fast')
multiplot('pearson', 'modeller_slow')
multiplot('pearson', 'rosetta')
multiplot('correlation', 'modeller_veryfast')
multiplot('correlation', 'modeller_fast')
multiplot('correlation', 'modeller_slow')
multiplot('correlation', 'rosetta')


