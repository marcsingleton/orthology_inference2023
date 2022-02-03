"""Make meta alignments from amino acid alignments of genes."""

import os
from collections import namedtuple
from random import randrange, seed

from src.utils import read_fasta


def is_redundant(col, cutoff):
    gapnum = 0
    for _, sym in col:
        if sym in ['-', 'X']:
            gapnum += 1
    return gapnum <= (1 - cutoff) * len(col)


def is_invariant(col):
    for i, (_, sym) in enumerate(col):
        if sym not in ['-', 'X']:
            break
    for _, sym in col[i+1:]:
        if sym not in ['-', 'X', col[i].sym]:
            return False
    return True


Column = namedtuple('Column', ['spid', 'sym'])
seed(930715)  # Set seed to make results consistent

# Extract column pools
colpools = [('100red', lambda col: is_redundant(col, 1), []),
            ('100red_ni', lambda col: is_redundant(col, 1) and not is_invariant(col), []),
            ('50red', lambda col: is_redundant(col, 0.5), []),
            ('50red_ni', lambda col: is_redundant(col, 0.5) and not is_invariant(col), []),
            ('0red', lambda col: is_redundant(col, 0), []),
            ('0red_ni', lambda col: is_redundant(col, 0) and not is_invariant(col), [])]
for file_id in filter(lambda x: x.endswith('.mfa'), os.listdir('../align_fastas/out/')):  # Because inputs are not sorted, results are not guaranteed to be consistent
    msa = read_fasta(f'../align_fastas/out/{file_id}')
    for i in range(len(msa[0][1])):
        col = [Column(header[-4:], seq[i]) for header, seq in msa]
        for _, condition, colpool in colpools:
            if condition(col):
                colpool.append(col)

# Make meta alignments
for label, _, colpool in colpools:
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

        with open(f'out/{label}/meta_{samplenum}.fasta', 'w') as file:
            for spid, seq in sorted(seqs.items()):
                seqstring = '\n'.join([''.join(seq[i:i+80]) for i in range(0, len(seq), 80)]) + '\n'
                file.write(f'>{spid} {label}_{samplenum}\n' + seqstring)

"""
OUTPUT
100red: 2240900
100red_ni: 853598
50red: 2526699
50red_ni: 1072823
0red: 2937882
0red_ni: 1219886

DEPENDENCIES
../align_fastas/align_fastas.py
    ../align_fastas/out/*.mfa
"""