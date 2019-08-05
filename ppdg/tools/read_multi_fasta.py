#!/usr/bin/env python3

def read_multi_fasta(fname):
    '''
        Read a fasta file containing multiple sequences.
        Return a dictionary where the key is the ID of the sequence
    '''
    seqs = dict()
    seqid = None
    seq = ''
    with open(fname, 'r') as fp:
        for line in fp.readlines():
            line = line.strip()
            if len(line)==0:
                continue
            if line[0]=='>':
                if seq and seqid:
                    if seqid in seqs.keys():
                        log.error('Sequence ID %s is duplicate in file %s!' % (seqid, fname))
                        return None
                    seqs[seqid] = seq
                seqid = line.split()[0][1:]
                seq = ''
            else:
                seq += line
    if seq and seqid:
        if seqid in seqs.keys():
            log.error('Sequence ID %s is duplicate in file %s!' % (seqid, fname))
            return None
        seqs[seqid] = seq
    return seqs

