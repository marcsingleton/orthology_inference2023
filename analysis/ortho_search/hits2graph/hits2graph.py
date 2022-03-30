"""Convert hits to a directed graph."""

import os
from itertools import groupby, permutations


def line2key(line):
    fields = line.rstrip('\n').split('\t')
    return fields[0]


columns = {'qppid': str, 'qgnid': str, 'qspid': str,
           'sppid': str, 'sgnid': str, 'sspid': str,
           'hspnum': int, 'chspnum': int,
           'qlen': int, 'nqa': int, 'cnqa': int,
           'slen': int, 'nsa': int, 'cnsa': int,
           'bitscore': float}

# Load genomes
spids = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spids.append(line.rstrip('\n').split('\t')[0])

# Make graph
graph = {}
for qspid, sspid in permutations(spids, 2):
    with open(f'../hsps2hits/out/{qspid}/{sspid}.tsv') as file:
        file.readline()  # Skip header
        for _, group in groupby(file, line2key):
            group = [line.rstrip('\n').split('\t') for line in group]
            max_bitscore = max([float(fields[14]) for fields in group])
            for fields in group:
                if max_bitscore == float(fields[14]):  # Only record hits with maximum bitscore
                    hit = {column: f(field) for (column, f), field in zip(columns.items(), fields)}
                    qppid, sppid = hit['qppid'], hit['sppid']
                    qlen, cnqa = hit['qlen'], hit['cnqa']
                    bitscore = hit['bitscore']

                    if cnqa / qlen >= 0.5:
                        try:
                            graph[qppid].add((sppid, bitscore))
                        except KeyError:
                            graph[qppid] = {(sppid, bitscore)}

# Write to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/hit_graph.tsv', 'w') as file:
    for qppid, edges in graph.items():
        file.write(qppid + '\t' + ','.join([sppid + ':' + str(bitscore) for sppid, bitscore in edges]) + '\n')

"""
DEPENDENCIES
../config/genomes.tsv
../hsps2hits/hsps2hits.py
    ../hsps2hits/out/*/*.tsv
"""