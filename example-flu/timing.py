#!/usr/bin/env python3

import json
import math
import numpy as np
import matplotlib.pyplot as plt

def main():

    ## Time to compute descriptors
    protocol = 'modeller_fast'
    with open('descriptors-all.json', 'r') as fp:
        alldata = json.load(fp)
    desclist = set([ k for name in alldata.keys() for k in alldata[name][protocol].keys() if k.startswith(">TIME") and k!=">TIME_makemodel" and k!='>TIME_prodigy_ic_nis' ])

    descs = list()
    vals = list()
    errs = list()
    for desc in sorted(desclist):
        v = [ np.average(list(data[protocol][desc].values())) for data in alldata.values() if desc in data[protocol].keys() ]
        descs.append(desc.replace(">TIME_", "") if desc!='>TIME_intermolecular_contacts' else 'IC')
        vals.append(np.average(v))
        errs.append(np.std(v))

    plt.bar(range(len(descs)), vals, yerr=errs, align='center', alpha=0.5, ecolor='black', capsize=5, log=True)
    plt.ylabel("Time [s]")
    plt.xticks(range(len(descs)),descs,rotation='vertical')
    plt.savefig('fig/timing.png', bbox_inches='tight', dpi=300)


    ## Time to generate models
    with open('descriptors-all.json') as fp:
        data = json.load(fp)

    for protocol in ['modeller_veryfast', 'modeller_fast', 'modeller_slow', 'rosetta']:
        timelst = list()
        for cpx, info in data.items():
            for num, time in info[protocol]['>TIME_makemodel'].items():
                timelst.append(time)
        print('%-20s %10.3f %8.3f' % (protocol, np.average(timelst)/60, np.std(timelst)/60))
    
main()



