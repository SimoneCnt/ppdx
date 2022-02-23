#!/usr/bin/env python3

import json
import random
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

import sklearn.ensemble
import sklearn.model_selection
import sklearn.preprocessing
import sklearn.metrics


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

def fit(X, y, random_seed1, random_seed2):
    mdlfunc = sklearn.ensemble.RandomForestRegressor
    mdlparams = {"n_estimators": 51, "max_depth": 10, "min_samples_leaf" : 2}
    mdl = mdlfunc(random_state=random_seed2, **mdlparams)
    X_train, X_test, y_train, y_test = sklearn.model_selection.train_test_split(X, y, test_size=0.5, shuffle=True, random_state=random_seed1)
    scaler = sklearn.preprocessing.RobustScaler().fit(X_train)
    mdl.fit(scaler.transform(X_train), y_train)
    y_predict = mdl.predict(scaler.transform(X_test))
    y_check = mdl.predict(scaler.transform(X_train))
    return np.array(y_test), np.array(y_predict), np.array(y_train), np.array(y_check)

def get_rp(y_test, y_predict, y_train, y_check):
    rp_test = scipy.stats.pearsonr(y_test, y_predict)[0]
    rp_train = scipy.stats.pearsonr(y_train, y_check)[0]
    return rp_test #, rp_train

def pearson_random_sampling(replica=1000):
    """
        Sample multiple models based on random split of train and test sets.
    """
    descriptors = ['BSA', 'BSA_A', 'BSA_P', 'BSA_C', 'NIS_A', 'NIS_P', 'NIS_C']
    protocol = 'modeller_veryfast'
    names, y = get_experimental_data()
    X = np.loadtxt('rf-data-ref.txt')
    rp = np.array([get_rp(*fit(X, y, random.randrange(100000), random.randrange(100000))) for i in range(replica) ])
    avg = np.average(rp)
    std = np.std(rp)
    plt.hist(rp, bins='auto', label='Average: %.2f\nStandard deviation: %.2f' % (avg, std))
    plt.xlabel('Pearson correlation coefficient')
    plt.ylabel('Count')
    plt.legend()
    plt.savefig('fig/rf-pearson-error-%d.png' % (replica), bbox_inches='tight', dpi=300)
    plt.clf()
    return


def bootstrap(replica=10000):
    """
        Compute the error in the pearson correlation coefficient by boostrap.
        Compute it for the "cherry-picked" model and for a random one.
    """

    descriptors = ['BSA', 'BSA_A', 'BSA_P', 'BSA_C', 'NIS_A', 'NIS_P', 'NIS_C']
    protocol = 'modeller_veryfast'
    names, y = get_experimental_data()
    X = np.loadtxt('rf-data-ref.txt')

    print("Cherry-Picked model")
    random_seed1 = 52871
    random_seed2 = 41242
    y_test, y_predict, y_train, y_check = fit(X, y, random_seed1, random_seed2)
    rp_test = scipy.stats.pearsonr(y_test, y_predict)[0]
    rp_train = scipy.stats.pearsonr(y_train, y_check)[0]
    print("rp_test", rp_test)
    print("rp_train", rp_train)
    nsamples = len(y_test)
    rp = [ scipy.stats.pearsonr(y_test[idx], y_predict[idx])[0] for idx in np.random.choice(nsamples, size=(replica,nsamples), replace=True) ]
    avg = np.mean(rp)
    std = np.std(rp)
    print("Bootstrap", avg, std)
    plt.hist(rp, bins='auto', alpha=0.7, label='Average: %.2f\nStandard deviation: %.2f' % (avg, std))

    print("Random model")
    random_seed1 = random.randrange(100000)
    random_seed2 = random.randrange(100000)
    y_test, y_predict, y_train, y_check = fit(X, y, random_seed1, random_seed2)
    rp_test = scipy.stats.pearsonr(y_test, y_predict)[0]
    rp_train = scipy.stats.pearsonr(y_train, y_check)[0]
    print("rp_test", rp_test)
    print("rp_train", rp_train)
    nsamples = len(y_test)
    rp = [ scipy.stats.pearsonr(y_test[idx], y_predict[idx])[0] for idx in np.random.choice(nsamples, size=(replica,nsamples), replace=True) ]
    avg = np.mean(rp)
    std = np.std(rp)
    print("Bootstrap", avg, std)
    plt.hist(rp, bins='auto', alpha=0.7, label='Average: %.2f\nStandard deviation: %.2f' % (avg, std))

    plt.xlabel('Pearson correlation coefficient')
    plt.ylabel('Count')
    plt.legend()
    plt.savefig('fig/rf-bootstrap-%d.png' % (replica), bbox_inches='tight', dpi=300)
    plt.clf()


pearson_random_sampling()
bootstrap()

