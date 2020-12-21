"""Make meta alignments from nucleotide alignments of genes."""

import os
from collections import namedtuple
from random import randrange

import Bio.AlignIO as AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


def is_redundant(col, cutoff):
    gapnum = 0
    for _, sym in col:
        if sym in ['-', 'N']:
            gapnum += 1
    return gapnum <= (1 - cutoff) * len(col)


def is_invariant(col):
    for i, (_, sym) in enumerate(col):
        if sym not in ['-', 'N']:
            break
    for _, sym in col[i+1:]:
        if sym not in ['-', 'N', col[i].sym]:
            return False
    return True


Column = namedtuple('Column', ['spid', 'sym'])

# Extract column pools
colpools = {'100red': (lambda col: is_redundant(col, 1), []),
            '100red_ni': (lambda col: is_redundant(col, 1) and not is_invariant(col), []),
            '50red': (lambda col: is_redundant(col, 0.5), []),
            '50red_ni': (lambda col: is_redundant(col, 0.5) and not is_invariant(col), []),
            '0red': (lambda col: is_redundant(col, 0), []),
            '0red_ni': (lambda col: is_redundant(col, 0) and not is_invariant(col), [])}
for file_id in filter(lambda x: x.endswith('.mfa'), os.listdir('../align_aa2nt1/out/')):
    align = AlignIO.read(f'../align_aa2nt1/out/{file_id}', 'fasta')
    for i in range(len(align[0])):
        col = [Column(seq.description[-4:], seq[i]) for seq in align]
        for condition, colpool in colpools.values():
            if condition(col):
                colpool.append(col)

# Make meta alignments
for label, (_, colpool) in colpools.items():
    if not os.path.exists(f'out/{label}/'):
        os.makedirs(f'out/{label}/')

    print(f'{label}:', len(colpool))
    for samplenum in range(100):
        sample = [colpool[randrange(len(colpool))] for _ in range(10000)]
        seqs = {}
        for col in sample:
            for spid, sym in col:
                try:
                    seqs[spid].append(sym)
                except KeyError:
                    seqs[spid] = [sym]

        align = MultipleSeqAlignment([SeqRecord(Seq(''.join(seq)), id=spid, description=f'{label}_{samplenum}')
                                      for spid, seq in seqs.items()])
        AlignIO.write(align, f'out/{label}/meta_{samplenum}.fasta', 'fasta')
        AlignIO.write(align, f'out/{label}/meta_{samplenum}.phy', 'phylip-relaxed')

"""
OUTPUT
100red: 2460073
100red_ni: 1273799
50red: 2646503
50red_ni: 1424250
0red: 3002523
0red_ni: 1515908

DEPENDENCIES
../align_aa2nt1/align_aa2nt1.py
    ../align_aa2nt1/out/*.mfa
"""