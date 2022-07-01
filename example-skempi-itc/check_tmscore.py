#!/usr/bin/env python3

import json
import numpy as np

# Check the TM-align score for all models
# TMscores less than 0.6/0.7 indicate the model is structurally
# different from the template. Could be a sign of a glitch in the
# modeling or an error somewhere.

def main():
    with open('descriptors-all.json', 'r') as fp:
        data = json.load(fp)
    for cpx, protocols in data.items():
        for protocol, descriptors in protocols.items():
            amin = np.amin(list(descriptors['TMscore'].values()))
            if amin<0.8:
                print(cpx, protocol, amin)

main()

