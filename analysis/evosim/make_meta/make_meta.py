"""Make meta alignments of ordered and disordered regions."""

import os
import re

import numpy.random
from src.utils import read_fasta


def is_redundant(column, cutoff):
    count = 0
    for sym in column:
        if sym in ['-', '.', 'X']:
            count += 1
    return count <= (1 - cutoff) * len(column)


rng = numpy.random.default_rng(seed=930715)
num_columns = int(1E5)  # Scientific notation defaults to float
spid_regex = r'spid=([a-z]+)'

# Load genomes
spids = []
with open('../../ortho_MSA/config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spids.append(line.rstrip('\n').split('\t')[0])
spid2idx = {spid: i for i, spid in enumerate(spids)}

# Load regions
OGid2regions = {}
with open('../../brownian/aucpred_regions/out/regions.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder = line.rstrip('\n').split('\t')
        try:
            OGid2regions[OGid].append((int(start), int(stop), disorder))
        except KeyError:
            OGid2regions[OGid] = [(int(start), int(stop), disorder)]

# Extract column pools
column_pools = [('100R_disorder', True, lambda column: is_redundant(column, 1), []),
                ('100R_order', False, lambda column: is_redundant(column, 1), []),
                ('50R_disorder', True, lambda column: is_redundant(column, 0.5), []),
                ('50R_order', False, lambda column: is_redundant(column, 0.5), []),
                ('0R_disorder', True, lambda column: is_redundant(column, 0), []),
                ('0R_order', False, lambda column: is_redundant(column, 0), [])]
for OGid, regions in sorted(OGid2regions.items()):
    msa = [(re.search(spid_regex, header).group(1), seq) for header, seq in read_fasta(
        f'../../brownian/insertion_trim/out/{OGid}.afa')]
    msa = sorted(msa, key=lambda x: spid2idx[x[0]])
    if len(msa) < len(spids):  # Only use alignments with all species
        continue
    for start, stop, region_disorder in regions:
        for i in range(start, stop):
            column = [seq[i] for _, seq in msa]
            for _, pool_disorder, condition, column_pool in column_pools:
                if (region_disorder is pool_disorder) and condition(column):
                    column_pool.append(column)

if not os.path.exists('out/'):
    os.mkdir('out/')

# Write column pool sizes
with open('out/pool_sizes.tsv', 'w') as file:
    file.write('pool_name\tpool_size\n')
    for label, _, column_pool in column_pools:
        file.write(f'{label}\t{len(column_pool)}\n')

# Make meta alignments
for label, _, _, column_pool in column_pools:
    sample = rng.choice(column_pool, size=num_columns)
    seqs = {spid: [] for spid in spids}
    for column in sample:
        for i, sym in enumerate(column):
            seqs[spids[i]].append(sym)

    with open(f'out/{label}.afa', 'w') as file:
        for spid, seq in sorted(seqs.items()):
            seqstring = '\n'.join([''.join(seq[i:i+80]) for i in range(0, len(seq), 80)])
            file.write(f'>{spid} {label}\n{seqstring}\n')

"""
DEPENDENCIES
../../brownian/aucpred_regions/get_regions.py
    ../..brownian/aucpred_regions/out/regions.tsv
../../brownian/insertion_trim/extract.py
    ../../brownian/insertion_trim/out/*.afa
../../ortho_MSA/config/genomes.tsv
"""