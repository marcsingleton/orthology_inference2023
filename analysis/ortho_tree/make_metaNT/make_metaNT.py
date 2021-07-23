"""Make meta alignments from nucleotide alignments of genes."""

import os
from collections import namedtuple
from random import randrange, seed

from skbio import io, TabularMSA, Protein


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
seed(930715)  # Set seed to make results consistent

# Extract column pools
colpools = {'100red': (lambda col: is_redundant(col, 1), []),
            '100red_ni': (lambda col: is_redundant(col, 1) and not is_invariant(col), []),
            '50red': (lambda col: is_redundant(col, 0.5), []),
            '50red_ni': (lambda col: is_redundant(col, 0.5) and not is_invariant(col), []),
            '0red': (lambda col: is_redundant(col, 0), []),
            '0red_ni': (lambda col: is_redundant(col, 0) and not is_invariant(col), [])}
for file_id in filter(lambda x: x.endswith('.mfa'), os.listdir('../align_aa2nt/out/')):
    msa = io.read(f'../align_aa2nt/out/{file_id}', 'fasta', into=TabularMSA, constructor=Protein)
    for i in range(msa.shape[1]):
        col = [Column(seq.metadata['id'][-4:], str(seq)[i]) for seq in msa]
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

        msa = TabularMSA([Protein(''.join(seq), metadata={'id': spid, 'description': f'{label}_{samplenum}'})
                          for spid, seq in seqs.items()], minter='id')  # Sets index to id (for writing in Phylip)
        msa.write(f'out/{label}/meta_{samplenum}.fasta', 'fasta', max_width=80)

"""
OUTPUT
100red: 4434228
100red_ni: 2274285
50red: 4966029
50red_ni: 2684624
0red: 5759178
0red_ni: 2933188

DEPENDENCIES
../align_aa2nt/align_aa2nt.py
    ../align_aa2nt/out/*.mfa
"""