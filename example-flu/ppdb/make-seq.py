#!/usr/bin/env python3

import ppdg

import Bio.Seq
def seq_nt2aa(seq_nt):
    """
        Convert a nucleotide sequence to amino acids.
    """
    return str(Bio.Seq.Seq(seq_nt).translate())


def main():

    # Read sequences and database
    sequences_ag = ppdg.tools.read_multi_fasta('ag-cleaned.seq')
    sequences_ab = ppdg.tools.read_multi_fasta('kallewaard.seq')

    # Prepare list of what to compute
    inputs = list()
    with open('ppdb.txt') as fp:
        for line in fp:
            if line[0]=='#':
                continue
            name, recn, lign, dgexp, tpl = line.split()
            ab, ag = name.split('__')
            ab = ab.replace('_', '-')
            sequence = sequences_ag[ag] + '/'
            sequence = sequence*3
            sequence = sequence + seq_nt2aa(sequences_ab[ab+'-VH']) + '/' + seq_nt2aa(sequences_ab[ab+'-VK'])
            print(">%s\n%s\n" % (name, sequence))

main()

