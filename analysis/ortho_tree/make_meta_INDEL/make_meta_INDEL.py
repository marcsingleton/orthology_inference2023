"""Make meta alignments from amino acid alignments of genes."""

import os
import re

import numpy as np
from src.utils import read_fasta


def is_invariant(column):
    for sym in column:
        if sym == '0':
            return False
    return True


rng = np.random.default_rng(seed=930715)
num_samples = 100
num_columns = 10000
spid_regex = r'spid=([a-z]+)'

# Load genomes
spids = []
with open('../config/genomes.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        spids.append(fields['spid'])
spid2idx = {spid: i for i, spid in enumerate(spids)}

# Extract column pools
column_pools = [('all', lambda x: True, []),
                ('NI', lambda x: not is_invariant(x), [])]
paths = sorted([path for path in os.listdir('../align_fastas/out/') if path.endswith('.afa')])  # Sort to ensure consistent order
for path in paths:
    msa = [(re.search(spid_regex, header).group(1), seq) for header, seq in read_fasta(f'../align_fastas/out/{path}')]
    msa = sorted(msa, key=lambda x: spid2idx[x[0]])
    for i in range(len(msa[0][1])):
        column = ['0' if seq[i] == '-' else '1' for _, seq in msa]
        for _, condition, column_pool in column_pools:
            if condition(column):
                column_pool.append(column)

if not os.path.exists('out/'):
    os.mkdir('out/')

# Write column pool sizes
with open('out/pool_sizes.tsv', 'w') as file:
    file.write('pool_name\tpool_size\n')
    for label, _, column_pool in column_pools:
        file.write(f'{label}\t{len(column_pool)}\n')

# Make meta alignments
for label, _, column_pool in column_pools:
    if not os.path.exists(f'out/{label}/'):
        os.mkdir(f'out/{label}/')

    for sample_id in range(num_samples):
        sample = rng.choice(column_pool, size=num_columns)
        seqs = {spid: [] for spid in spids}
        for column in sample:
            for i, sym in enumerate(column):
                seqs[spids[i]].append(sym)

        with open(f'out/{label}/{label}_{sample_id:03}.afa', 'w') as file:
            for spid, seq in sorted(seqs.items()):
                seqstring = '\n'.join([''.join(seq[i:i+80]) for i in range(0, len(seq), 80)])
                file.write(f'>{spid} {label}_{sample_id:03}\n{seqstring}\n')

"""
DEPENDENCIES
../config/genomes.tsv
../align_fastas/align_fastas.py
    ../align_fastas/out/*.afa
"""