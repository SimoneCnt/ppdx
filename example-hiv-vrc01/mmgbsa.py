#!/usr/bin/env python3

import json
import math
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

import sklearn.ensemble
import sklearn.model_selection
import sklearn.preprocessing
import sklearn.metrics
import joblib


def get_experimental_data(ppdbf='ppdb/ppdb.txt'):
    """
        Load experimental data to fit. Returns the name of the protein-protein
        complexes and their log(IC50) - only if the value is exact, no "greater
        than" values.
    """
    names = list()
    exps = list()
    with open(ppdbf, 'r') as fp:
        for line in fp:
            if line.startswith('#'): continue
            splt = line.split()
            name, dg = splt[0], splt[3]
            if dg[0]=='>': continue
            names.append(name)
            exps.append(float(dg))
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


def main(solvent, protocol):

    print(solvent, protocol)

    # List of descriptors to use
    descriptors = ['%s_ELEC' % (solvent), '%s_VDW' % (solvent), '%s_GB' % (solvent), '%s_ASP' % (solvent)]

    # Names of the complexes and experimental data
    names, y = get_experimental_data()

    # Descriptors
    X = get_data_from_json('descriptors-med.json', names, descriptors, protocol)

    # Linear Regression
    reg = sklearn.linear_model.LinearRegression().fit(X, y)
    y_pred = reg.predict(X)

    # Fitted coefficients
    for n, c in zip(descriptors, reg.coef_):
        print('\t%-10s %10.8f' % (n, c))
    print('\t%-10s %10.8f' % ('Intercept', reg.intercept_))

    # Compute some statistics
    rp = scipy.stats.pearsonr(y, y_pred)[0]
    rs = scipy.stats.spearmanr(y, y_pred)[0]
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(y, y_pred)
    rmse = np.sqrt(sklearn.metrics.mean_squared_error(y, y_pred))
    mae = sklearn.metrics.mean_absolute_error(y, y_pred)
    print("\tPearson: %.3f\n\tSpearman: %.3f" % (rp, rs))
    print("\tSlope: %.3f\n\tIntercept: %.3f\n\tPearson: %.3f\n\tP-value: %.3e\n\tStdErr: %.3f" % (slope, intercept, r_value, p_value, std_err))
    print("\tRMSE: %.3f\n\tMAE: %.3f" % (rmse, mae))

    # Plot the correlation
    y1 = [ yy for yy, name in zip(y, names) if name.startswith('93th057') ]
    yp1 = [ yy for yy, name in zip(y_pred, names) if name.startswith('93th057') ]
    y2 = [ yy for yy, name in zip(y, names) if name.startswith('rsc3') ]
    yp2 = [ yy for yy, name in zip(y_pred, names) if name.startswith('rsc3') ]
    print(len(y1), len(yp1), len(y2), len(yp2))
    #plt.plot(y, y_pred, '.', label='rp: %.2f' % (rp))
    plt.plot(y1, yp1, '.', color='blue')
    plt.plot(y2, yp2, '.', color='red')
    lins = np.linspace(np.min(y), np.max(y), 100)
    plt.plot(lins, slope*lins+intercept, '-', label='rp: %.2f' % (rp))
    plt.plot(lins, lins, ":", color='grey')
    plt.legend()

    yl = list(y)
    ypl = list(y_pred)
    rpl = list()
    for i in range(len(y)):
        ym = yl[:i]+y[i+1:]
        ypm = ypl[:i]+ypl[i+1:]
        rp = scipy.stats.pearsonr(ym, ypm)[0]
        rpl.append(rp)
    print('rp:', np.mean(rpl), np.std(rpl))

def main2():
    numd = 4
    ncols = 2
    nrows = math.ceil(numd/ncols)
    tmp = 3
    plt.figure(figsize=(tmp*ncols,tmp*nrows-0.5))
    #plt.subplots_adjust(hspace=0.4, wspace=0.3)

    plt.subplot(nrows, ncols, 1)
    plt.title('FACTS')
    plt.ylabel('Rosetta')
    main('FACTS', 'rosetta')

    plt.subplot(nrows, ncols, 2)
    plt.title('GBSW')
    main('GBSW', 'rosetta')

    plt.subplot(nrows, ncols, 3)
    plt.ylabel('Modeller')
    main('FACTS', 'modeller_fast')

    plt.subplot(nrows, ncols, 4)
    main('GBSW', 'modeller_fast')

    plt.savefig('fig/mmgbsa-regression.png', bbox_inches='tight', dpi=300)
    plt.clf()

main2()

