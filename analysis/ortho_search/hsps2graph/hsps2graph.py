"""Convert HSPs to a directed graph."""

import os
from itertools import groupby, permutations


def make_line2key(field_names):
    def line2key(line):
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        return fields['qppid'], fields['sppid']
    return line2key


# Load genomes
spids = []
with open('../config/genomes.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        spids.append(fields['spid'])

# Make graph
graph = {}
for qspid, sspid in permutations(spids, 2):
    with open(f'../blast2hsps/out/hsps/{qspid}/{sspid}.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        line2key = make_line2key(field_names)
        for key, _ in groupby(file, key=line2key):
            qppid, sppid = key
            try:
                graph[qppid].add(sppid)
            except KeyError:
                graph[qppid] = {sppid}

# Write to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/hsp_graph.tsv', 'w') as file:
    for qppid, sppids in graph.items():
        file.write(qppid + '\t' + ','.join(sppids) + '\n')

"""
DEPENDENCIES
../config/genomes.tsv
../blast2hsps/blast2hsps.py
    ../blast2hsps/out/hsps/*/*.tsv
"""