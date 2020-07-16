"""Remove non-reciprocal hits between polypeptides in polypeptide graph."""

import os

# Parse best hits as graph
ipgraph = {}
with open('../blast2pgraph/out/pgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        ipgraph[node] = adjs.split(',')

# Copy graph, including only reciprocal hits
opgraph = {}
for qppid, sppids in ipgraph.items():
    for sppid in sppids:
        try:
            reciprocal = qppid in ipgraph[sppid]
        except KeyError:
            reciprocal = False

        if reciprocal:
            try:
                opgraph[qppid].append(sppid)
            except KeyError:
                opgraph[qppid] = [sppid]

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write graph as adjacency list to file
with open('out/pgraph.tsv', 'w') as outfile:
    for qppid, sppids in opgraph.items():
        outfile.write(qppid + '\t' + ','.join(sppids) + '\n')

"""
DEPENDENCIES
../blast2pgraph/blast2pgraph.py
    ../blast2pgraph/out/pgraph.tsv
"""