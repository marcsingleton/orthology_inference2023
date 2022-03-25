"""Convert HSPs to a directed graph."""

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

# Make graph
graph = {}
for qspid, sspid in permutations(spids, 2):
    with open(f'../blast2hsps/out/hsps/{qspid}/{sspid}.tsv') as file:
        file.readline()  # Skip header
        for key, _ in groupby(file, key=line2key):
            qppid, sppid = key
            try:
                graph[qppid].add(sppid)
            except KeyError:
                graph[qppid] = {sppid}

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write to file
with open('out/hsp_graph.tsv', 'w') as file:
    for qppid, sppids in graph.items():
        file.write(qppid + '\t' + ','.join(sppids) + '\n')

"""
DEPENDENCIES
../config/genomes.tsv
../blast2hsps/blast2hsps.py
    ../blast2hsps/out/hsps/*/*.tsv
"""