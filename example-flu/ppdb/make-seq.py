#!/usr/bin/env python3

import ppdx

import Bio.Seq
def seq_nt2aa(seq_nt):
    """
        Convert a nucleotide sequence to amino acids.
    """
    return str(Bio.Seq.Seq(seq_nt).translate())


def main():

    # Read sequences and database
    sequences_ag = ppdx.tools.read_multi_fasta('ag-cleaned.seq')
    sequences_ab = ppdx.tools.read_multi_fasta('kallewaard.seq')
    sequences_ab.update(ppdx.tools.read_multi_fasta('corti.seq'))

    # Check that the length of the nucleotide sequence is a multiple of 3
    # Some VH sequences had an extra G at the end
    # Some VK sequences had an extra C at the end
    # Those have been modified manually in the sequence files
    for k,v in sequences_ab.items():
        if len(v)%3!=0:
            raise ValueError('Nucleotide sequence length for %s is not a multiple of 3 (%d%%3=%d)' % (k, len(v), len(v)%3))

    # Prepare list of what to compute
    inputs = list()
    with open('ppdb.txt', 'r') as fpin:
        with open('ppdb.seq', 'w') as fpout:
            for line in fpin:
                if line[0]=='#':
                    continue
                name, recn, lign, dgexp, tpl = line.split()
                ab, ag = name.split('__')
                ab = ab.replace('_', '-')
                sequence = sequences_ag[ag] + '/'
                sequence = sequence*3
                sequence = sequence + seq_nt2aa(sequences_ab[ab+'-VH']) + '/' + seq_nt2aa(sequences_ab[ab+'-VK'])
                fpout.write(">%s\n%s\n" % (name, sequence))

main()

