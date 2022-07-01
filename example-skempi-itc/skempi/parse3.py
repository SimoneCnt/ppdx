#!/usr/bin/env python3

import os
import csv
import numpy as np

import ppdx

ppdb = open('../ppdb/ppdb.txt', 'w')
ppdb.write("%-50s %6s %6s %15s %10s %10s %10s %10s %s\n" % ('#Name', 'RecN', 'LigN', 'Temperature', 'dG', 'dH', 'dS', 'TdS', 'Template'))
seqs = dict()

with open('gdata2.csv') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',', quotechar='"')
    for r in reader:
        splt = r['Name'].split('__')
        if len(splt)==3:
            name, mutations, temperature = splt
        else:
            name, temperature = splt
            mutations = None
        nameunq = name+'__'+mutations if mutations else name

        pdbf, rec, lig = name.split('_')

        pdb = ppdx.Pdb('PDBs/%s.pdb' % (pdbf))
        sequence = '/'.join([s for s in pdb.get_sequence().split('/') if s!=''])
        #pdb.keep_protein_only()
        pdb.remove_water_ions()
        chains = pdb.split_by_chain()
        chains_lst = list(chains.keys())
        if rec+lig != ''.join(chains_lst):
            raise ValueError('test1')
        if len(sequence.split('/'))!=len(rec+lig):
            raise ValueError('test2')
        for k in chains.keys():
            chains[k] = list(chains[k].get_sequence())
        if mutations:
            mutations = mutations.split('_')
            for m in mutations:
                fr = m[0]
                ch = m[1]
                to = m[-1]
                rr = m[2:-1]
                if m!=fr+ch+rr+to:
                    raise ValueError("%s, %s, %s %s %s" % (m, fr, ch, to, rr))
                rr = int(rr)-1
                frq = chains[ch][rr]
                if fr!=frq:
                    raise ValueError("%s %s" % (fr, frq))
                chains[ch][rr] = to

        newseq = '/'.join([''.join(chains[k]) for k in rec+lig])
        if nameunq in seqs.keys():
            if newseq!=seqs[nameunq]:
                raise ValueError('Sequences are different! %s' % (nameunq))
        else:
            seqs[nameunq] = newseq

        ppdb.write("%-50s %6d %6d %15d %10.3f %10.3f %10.3f %10.3f %s\n" % (nameunq, len(rec), len(lig), int(temperature), float(r['dG']), float(r['dH']), float(r['dS']), float(r['TdS']), pdbf+'.pdb'))

        pdb.write('../ppdb/%s.pdb' % (pdbf))

        if not mutations:
            if sequence!=newseq:
                print(r['Name'])
                print(sequence)
                print(newseq)
                quit()

with open('../ppdb/ppdb.seq', 'w') as fasta:
    for key, seq in seqs.items():
        fasta.write(">%s\n%s\n" % (key, seq))



