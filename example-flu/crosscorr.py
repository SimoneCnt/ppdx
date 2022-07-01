#!/usr/bin/env python3

import math
import json
import numpy as np
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.preprocessing import RobustScaler
import matplotlib
import matplotlib.pyplot as plt
import sklearn.decomposition
import sklearn.metrics


def plot_crosscorr(mtx, nice_names, fname=None):
    matplotlib.rcParams.update({'font.size': 8})
    fig = plt.figure(figsize=(13,13))
    ax = fig.add_subplot(111)
    cax = ax.matshow(mtx, interpolation='nearest', cmap=plt.get_cmap('Spectral'), vmin=-1, vmax=1)
    fig.colorbar(cax)
    ind_array = np.arange(mtx.shape[0])
    x, y = np.meshgrid(ind_array, ind_array)
    #for x_val, y_val in zip(x.flatten(), y.flatten()):
    #    c = str(np.around(mtx[x_val][y_val], 2))
    #    ax.text(x_val, y_val, c, va='center', ha='center', fontsize=4)
    ax.set_xticklabels(nice_names, rotation=90)
    ax.set_yticklabels(nice_names)
    ax.set_xticks(ind_array)
    ax.set_yticks(ind_array)
    if fname:
        plt.savefig(fname, bbox_inches='tight', dpi=300)
    else:
        plt.show()
    plt.clf()


def main():

    protocol = 'modeller_fast'

    with open('descriptors-avg.json', 'r') as fp:
        alldata = json.load(fp)

    desc = set()
    for cpx, cpxdata in alldata.items():
        desc |= set(cpxdata[protocol].keys())
    desc = sorted([ d for d in desc if not d.startswith('>') and d!='NRES' and d!='AGBNP' and d!='GBMV_POL' and d!='SOAP-Protein-OD' ])
    print('Number of descriptors:', len(desc))
    cpxs = [ c for c in alldata.keys() if c.startswith('FY')]

    data = np.zeros((len(cpxs), len(desc)), dtype=float)
    for i,d in enumerate(desc):
        for j,c in enumerate(cpxs):
            data[j][i] = alldata[c][protocol][d]

    scaler = RobustScaler().fit(data)
    X = scaler.transform(data)
    X = X.transpose()

    # Histograms
    numd = len(desc)
    ncols = 10
    nrows = math.ceil(numd/ncols)
    plt.figure(figsize=(3*ncols,3*nrows-0.5))
    plt.subplots_adjust(hspace=0.4, wspace=0.3)
    for n, d in enumerate(desc):
        plt.subplot(nrows, ncols, n+1)
        plt.title(d)
        plt.hist(X[n], bins='auto')
    plt.savefig('fig/histograms.png', bbox_inches='tight', dpi=300)
    plt.clf()

    # Dendogram
    method = 'complete'   # complete or average seem better
    Z = linkage(X, method=method, metric='correlation', optimal_ordering=True)
    fig = plt.figure(figsize=(6, 10))
    dn = dendrogram(Z, orientation='right', labels=desc)
    plt.savefig('fig/dendogram-%s.png' % (method), bbox_inches='tight', dpi=300)
    plt.clf()
    
    # Reorder based on dendogram
    labels = list(reversed(dn['ivl']))
    ndx = [desc.index(l) for l in labels]
    X = X[ndx, :]

    # Cross-correlation matrix
    size = len(desc)
    mtx = np.ones((size,size), dtype=float)
    for i in range(size):
        for j in range(i+1,size):
            rp = stats.pearsonr(X[i], X[j])[0]
            rs = stats.spearmanr(X[i], X[j])[0]
            mtx[i][j] = rp
            mtx[j][i] = rs
    plot_crosscorr(mtx, labels, 'fig/crosscorr-%s.png' % (method))

    # Explained_variance
    pca = sklearn.decomposition.PCA().fit(X.transpose())
    nc = pca.n_components_
    cumul = np.zeros(nc)
    for i in range(nc):
        cumul[i] = pca.explained_variance_ratio_[i]
        if (i>0): cumul[i] += cumul[i-1]
    cut = 25
    plt.figure(figsize=(6.4, 4.8))
    plt.plot(range(1, cut+1), cumul[:cut], 'ro-', label='Cumulative explained variance')
    plt.bar(range(1, cut+1), pca.explained_variance_ratio_[:cut], label='Explained variance ratio')
    plt.ylim(-0.05, 1.05)
    plt.xlabel('Number of components')
    plt.ylabel('Explained Variance')
    plt.xticks(range(1,cut+1))
    plt.legend()
    plt.savefig('fig/explained_variance.png', bbox_inches='tight', dpi=300)
    plt.clf()


main()

