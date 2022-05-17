"""Convert hits to a directed graph."""

import os
from itertools import groupby, permutations


def make_line2key(field_names):
    def line2key(line):
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        return fields['qppid']
    return line2key


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

# Make graph
graph = {}
for qspid, sspid in permutations(spids, 2):
    with open(f'../hsps2hits/out/{qspid}/{sspid}.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        line2key = make_line2key(field_names)
        for _, group in groupby(file, line2key):
            hits = []
            for line in group:
                hit = {key: columns[key](value) for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
                hits.append(hit)
            max_bitscore = max([hit['bitscore'] for hit in hits])

            for hit in hits:
                if max_bitscore == hit['bitscore']:
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