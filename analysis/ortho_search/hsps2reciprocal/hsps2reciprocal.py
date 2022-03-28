"""Compute reciprocity for each HSP."""

import os
from itertools import groupby, permutations


def line2key(line):
    fields = line.rstrip('\n').split('\t')
    return fields[0], fields[3]


# Load genomes
spids = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, _, _ = line.split()
        spids.append(spid)

# Load graph
graph = {}
with open('../hsps2graph/out/hsp_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = set(adjs.split(','))

# Check reciprocity
for qspid, sspid in permutations(spids, 2):
    rows = []
    with open(f'../blast2hsps/out/hsps/{qspid}/{sspid}') as file:
        file.readline()  # Skip header
        for key, group in groupby(file, key=line2key):
            qppid, sppid = key
            try:
                reciprocal = qppid in graph[sppid]
            except KeyError:
                reciprocal = False
            rows.extend([(qppid, sppid, str(reciprocal))] * len(list(group)))

    # Write to file
    if not os.path.exists(f'out/{qspid}/'):
        os.makedirs(f'out/{qspid}/')

    with open(f'out/{qspid}/{sspid}.tsv', 'w') as file:
        file.write('qppid\tsppid\treciprocal\n')
        for row in rows:
            file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../config/genomes.tsv
../blast2hsps/blast2hsps.py
    ../blast2hsps/hsps/out/*/*.tsv
../hsps2graph/hsps2graph.py
    ../hsps2graph/out/hsp_graph.tsv
"""