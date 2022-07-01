#!/usr/bin/env python3

import csv

with open('skempi_v2.csv') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=';', quotechar='"')
    data = [ row for row in reader ]


methods = set([d['Method'] for d in data])
print(methods)

keys = ['Affinity_mut_parsed', 'Affinity_wt_parsed', 'dH_mut (kcal mol^(-1))', 'dH_wt (kcal mol^(-1))', 'dS_mut (cal mol^(-1) K^(-1))', 'dS_wt (cal mol^(-1) K^(-1))']

gdata = list()
for d in data:
    skip = False
    for k in keys:
        if d[k]=="": skip=True
    if skip:
        if d['Method']=='ITC':
            print('SKIP', d)
    else:
        if d['Method']!='ITC':
            print('KEEP', d)
        gdata.append(d)

print(len(gdata))

with open('gdata.csv', 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=gdata[0].keys())
    writer.writeheader()
    for d in gdata:
        writer.writerow(d)


