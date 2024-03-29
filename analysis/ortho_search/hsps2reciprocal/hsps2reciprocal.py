"""Compute reciprocity for each HSP."""

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

# Load graph
graph = {}
with open('../hsps2graph/out/hsp_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        graph[node] = set(adjs.split(','))

# Check reciprocity
for qspid, sspid in permutations(spids, 2):
    rows = []
    with open(f'../blast2hsps/out/hsps/{qspid}/{sspid}.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        line2key = make_line2key(field_names)
        for key, group in groupby(file, key=line2key):
            qppid, sppid = key
            try:
                reciprocal = qppid in graph[sppid]
            except KeyError:
                reciprocal = False
            rows.extend([(qppid, sppid, str(reciprocal))] * len(list(group)))  # List multiply by scalar makes shallow copy of contents

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