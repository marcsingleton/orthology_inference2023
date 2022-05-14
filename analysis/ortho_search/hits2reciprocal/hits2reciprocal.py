"""Compute reciprocity for each hit."""

import os
from itertools import permutations


def is_reciprocal(qppid, sppid, graph):
    # Check both directions since not all hits are in graph
    try:
        reciprocal1 = qppid in graph[sppid]
        reciprocal2 = sppid in graph[qppid]
        reciprocal = reciprocal1 and reciprocal2
    except KeyError:
        reciprocal = False
    return reciprocal


columns = {'qppid': str, 'qgnid': str,
           'sppid': str, 'sgnid': str,
           'hspnum': int, 'chspnum': int,
           'qlen': int, 'nqa': int, 'cnqa': int,
           'slen': int, 'nsa': int, 'cnsa': int,
           'bitscore': float}

# Load genomes
spids = []
with open('../config/genomes.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        spids.append(fields['spid'])

# Load graph
graph = {}
with open('../hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = {adj.split(':')[0] for adj in adjs.split(',')}

# Check reciprocity
for qspid, sspid in permutations(spids, 2):
    rows = []
    with open(f'../hsps2hits/out/{qspid}/{sspid}.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            hit = {column: f(field) for (column, f), field in zip(columns.items(), line.rstrip('\n').split('\t'))}
            qppid, sppid = hit['qppid'], hit['sppid']
            reciprocal = is_reciprocal(qppid, sppid, graph)
            rows.append((qppid, sppid, str(reciprocal)))

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
../hsps2hits/hsps2hits.py
    ../hsps2hits/out/*/*.tsv
../hits2graph/hits2graph.py
    ../hits2graph/out/hit_graph.tsv
"""