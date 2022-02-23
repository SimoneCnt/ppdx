#!/usr/bin/env python3

import json
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


def main(ref=True):

    # List of descriptors to use
    descriptors = ['BSA', 'BSA_A', 'BSA_P', 'BSA_C', 'NIS_A', 'NIS_P', 'NIS_C']
    protocol = 'modeller_veryfast'

    # Names of the complexes and experimental data
    names, y = get_experimental_data()

    # Load computed descriptors (choose one)
    # or from reference data
    # or from your computed json
    if ref:
        X = np.loadtxt('rf-data-ref.txt')
    else:
        X = get_data_from_json('descriptors-med.json', names, descriptors, protocol)
        np.savetxt('rf-data-new.txt', X, header=" ".join(descriptors))

    # Setup the Random Forest model
    mdlfunc = sklearn.ensemble.RandomForestRegressor
    mdlparams = {"n_estimators": 51, "max_depth": 10, "min_samples_leaf" : 2}
    random_seed1 = 52871
    random_seed2 = 41242
    mdl = mdlfunc(random_state=random_seed2, **mdlparams)

    # Split train and test set, fit a scaler and fit the random forest model
    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.5, shuffle=True, random_state=random_seed1)
    scaler = sklearn.preprocessing.RobustScaler().fit(X_train)
    mdl.fit(scaler.transform(X_train), y_train)
    y_predict = mdl.predict(scaler.transform(X_test))
    y_check = mdl.predict(scaler.transform(X_train))

    # Compute some statistics
    rp = scipy.stats.pearsonr(y_train, y_check)[0]
    rs = scipy.stats.spearmanr(y_train, y_check)[0]
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(y_train, y_check)
    rmse = np.sqrt(sklearn.metrics.mean_squared_error(y_train, y_check))
    mae = sklearn.metrics.mean_absolute_error(y_train, y_check)
    print("Training set:")
    print("\tPearson: %.3f\n\tSpearman: %.3f" % (rp, rs))
    print("\tSlope: %.3f\n\tIntercept: %.3f\n\tPearson: %.3f\n\tP-value: %.3e\n\tStdErr: %.3f" % (slope, intercept, r_value, p_value, std_err))
    print("\tRMSE: %.3f\n\tMAE: %.3f" % (rmse, mae))


    # Compute some statistics
    rp = scipy.stats.pearsonr(y_test, y_predict)[0]
    rs = scipy.stats.spearmanr(y_test, y_predict)[0]
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(y_test, y_predict)
    rmse = np.sqrt(sklearn.metrics.mean_squared_error(y_test, y_predict))
    mae = sklearn.metrics.mean_absolute_error(y_test, y_predict)
    print("Validation set:")
    print("\tPearson: %.3f\n\tSpearman: %.3f" % (rp, rs))
    print("\tSlope: %.3f\n\tIntercept: %.3f\n\tPearson: %.3f\n\tP-value: %.3e\n\tStdErr: %.3f" % (slope, intercept, r_value, p_value, std_err))
    print("\tRMSE: %.3f\n\tMAE: %.3f" % (rmse, mae))

    # Plot the correlation
    plt.plot(y_test, y_predict, '.', label='Test')
    plt.plot(y_train, y_check, '.', label='Train')
    lins = np.linspace(np.min(y), np.max(y), 100)
    plt.plot(lins, slope*lins+intercept, '-')
    plt.plot(lins, lins, ":", color='grey')
    plt.legend()
    plt.xlabel('Experimental')
    plt.ylabel('Computed')
    plt.savefig('fig/rf-correlation-%s.png' % ("ref" if ref else "new"), bbox_inches='tight', dpi=300)
    plt.clf()

    # Save the python object file with all the info
    yp = mdl.predict(scaler.transform(X))
    joblib.dump([protocol, descriptors, scaler, mdl, X, yp], 'rf.joblib')


def data_compare():
    ref = np.loadtxt('rf-data-ref.txt').transpose()
    new = np.loadtxt('rf-data-new.txt').transpose()

    for i in range(6):
        plt.plot(ref[i], new[i], '.')
        amin = min(np.amin(ref[i]), np.amin(new[i]))
        amax = max(np.amax(ref[i]), np.amax(new[i]))
        l = np.linspace(amin, amax, 100)
        plt.plot(l, l)
        plt.show()

main(True)
main(False)
data_compare()

