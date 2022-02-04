"""Convert hits to a directed pgraph."""

import os
from itertools import groupby


def add_edge(qppid, sppid, bitscore, pgraph):
    try:
        pgraph[qppid].add((sppid, bitscore))
    except KeyError:
        pgraph[qppid] = {(sppid, bitscore)}


columns = {'qppid': str, 'qgnid': str, 'qspid': str,
           'sppid': str, 'sgnid': str, 'sspid': str,
           'hspnum': int, 'chspnum': int,
           'qlen': int, 'nqa': int, 'cnqa': int,
           'slen': int, 'nsa': int, 'cnsa': int,
           'bitscore': float}

# Make pgraphs
graph1 = {}
graph2 = {}
for qspid in os.listdir('../hsps2hits/out/'):
    for sspid in os.listdir(f'../hsps2hits/out/{qspid}/'):
        with open(f'../hsps2hits/out/{qspid}/{sspid}') as file:
            file.readline()  # Skip header
            for _, group in groupby(file, lambda x: x.split()[0]):
                group = list([line.split() for line in group])
                bs_max = max([float(fields[14]) for fields in group])
                for fields in group:
                    if bs_max == float(fields[14]):
                        d = {column: f(field) for (column, f), field in zip(columns.items(), fields)}
                        qppid, sppid = d['qppid'], d['sppid']
                        qlen, cnqa = d['qlen'], d['cnqa']
                        bitscore = d['bitscore']

                        add_edge(qppid, sppid, bitscore, graph1)
                        if cnqa / qlen >= 0.5:
                            add_edge(qppid, sppid, bitscore, graph2)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write to file
with open('out/pgraph1.tsv', 'w') as file:
    for qppid, edges in graph1.items():
        file.write(qppid + '\t' + ','.join([sppid + ':' + str(bitscore) for sppid, bitscore in edges]) + '\n')
with open('out/pgraph2.tsv', 'w') as file:
    for qppid, edges in graph2.items():
        file.write(qppid + '\t' + ','.join([sppid + ':' + str(bitscore) for sppid, bitscore in edges]) + '\n')

"""
DEPENDENCIES
../hsps2hits/hsps2hits.py
    ../hsps2hits/out/*/*.tsv
"""