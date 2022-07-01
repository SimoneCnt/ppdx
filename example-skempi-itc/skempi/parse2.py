#!/usr/bin/env python3

import csv
import numpy as np

gdata = list()
with open('gdata_cleaned.csv') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=';', quotechar='"')
    for r in reader:

        if r['PDB'].startswith('1LFD'): continue

        # Mutant
        temperature = float(r['Temperature'])
        name = r['PDB'] + '__' + '_'.join(r['Mutations'].split(',')) + '__%.0f' % (temperature)
        dG = 1.987*temperature*np.log(float(r['Affinity_mut']))/1000
        dH = float(r['dH_mut'])
        dS = float(r['dS_mut'])
        TdS = temperature * dS / 1000
        check_mut = np.abs(dG - (dH - TdS))
        mut = dict(Name=name, Temperature=temperature, dG=dG, dH=dH, dS=dS, TdS=TdS)

        # WT
        temperature = float(r['Temperature'])
        if r['Notes'].startswith('wild type parameters measured at 303K'):
            temperature = 303.0
        name = r['PDB'] + '__%.0f' % (temperature)
        dG = 1.987*temperature*np.log(float(r['Affinity_wt']))/1000
        dH = float(r['dH_wt'])
        dS = float(r['dS_wt'])
        TdS = temperature * dS / 1000
        check_wt = np.abs(dG - (dH - TdS))
        wt = dict(Name=name, Temperature=temperature, dG=dG, dH=dH, dS=dS, TdS=TdS)

        addmut = True
        addwt = True
        for v in gdata:
            if v==mut: addmut = False
            if v==wt: addwt = False
        for v in gdata:
            if addmut and v['Name']==mut['Name']:
                print('Duplicate\n', v, '\n',  mut)
                addmut = False
        for v in gdata:
            if addwt and v['Name']==wt['Name']:
                print('Duplicate\n', v, '\n', wt)
                addwt = False

        if check_mut>0.2: addmut=False
        if check_wt>0.2: addwt=False

        if addwt: gdata.append(wt)        
        if addmut: gdata.append(mut)

print(len(gdata))
print(len(set( [ v['Name'] for v in gdata ])))

with open('gdata2.csv', 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=gdata[0].keys())
    writer.writeheader()
    for d in gdata:
        writer.writerow(d)


