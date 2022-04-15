"""Make meta alignments from nucleotide alignments of genes."""

import os
import re

import numpy.random
from src.utils import read_fasta


def is_redundant(column, cutoff):
    count = 0
    for sym in column:
        if sym in ['-', 'X']:
            count += 1
    return count <= (1 - cutoff) * len(column)


def is_invariant(column):
    sym0 = None
    for sym in column:
        if (sym0 is not None) and (sym not in ['-', 'X', sym0]):
            return False
        elif sym not in ['-', 'X']:
            sym0 = sym
    return True


rng = numpy.random.default_rng(seed=930715)
num_samples = 100
num_columns = 10000
spid_regex = r'spid=([a-z]+)'

# Load genomes
spids = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spids.append(line.rstrip('\n').split('\t')[0])
spid2idx = {spid: i for i, spid in enumerate(spids)}

# Extract column pools
column_pools = [('100R', lambda x: is_redundant(x, 1), []),
                ('100R_NI', lambda x: is_redundant(x, 1) and not is_invariant(x), []),
                ('50R', lambda x: is_redundant(x, 0.5), []),
                ('50R_NI', lambda x: is_redundant(x, 0.5) and not is_invariant(x), []),
                ('0R', lambda x: is_redundant(x, 0), []),
                ('0R_NI', lambda x: is_redundant(x, 0) and not is_invariant(x), [])]
paths = sorted([path for path in os.listdir('../align_AA2NT/out/') if path.endswith('.afa')])  # Sort to ensure consistent order
for path in paths:
    msa = [(re.search(spid_regex, header).group(1), seq) for header, seq in read_fasta(f'../align_fastas/out/{path}')]
    msa = sorted(msa, key=lambda x: spid2idx[x[0]])
    for i in range(len(msa[0][1])):
        column = [seq[i] for _, seq in msa]
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
../align_AA2NT/align_AA2NT.py
    ../align_AA2NT/out/*.afa
"""