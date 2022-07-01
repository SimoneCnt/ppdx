#!/usr/bin/env python3

import json
import numpy as np

# Read starting data
with open('descriptors-all.json', 'r') as fp:
    data = json.load(fp)

# Save averages
avgs = dict()
for cpx in data.keys():
    avgs[cpx] = dict()
    for protocol in data[cpx].keys():
        avgs[cpx][protocol] = dict()
        for desc in data[cpx][protocol].keys():
            avgs[cpx][protocol][desc] = np.average(list(data[cpx][protocol][desc].values()))
with open('descriptors-avg.json', 'w') as fp:
    json.dump(avgs, fp, indent=4, sort_keys=True)

# Save medians
avgs = dict()
for cpx in data.keys():
    avgs[cpx] = dict()
    for protocol in data[cpx].keys():
        avgs[cpx][protocol] = dict()
        for desc in data[cpx][protocol].keys():
            avgs[cpx][protocol][desc] = np.median(list(data[cpx][protocol][desc].values()))
with open('descriptors-med.json', 'w') as fp:
    json.dump(avgs, fp, indent=4, sort_keys=True)

