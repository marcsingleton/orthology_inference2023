"""Find connected components of pgraph."""

import os
from DFS import connect

# Parse best hits as graph
pgraph = {}
with open('../make_xreciprocal/out/pgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        pgraph[node] = adjs.split(',')

# Cluster by triangle criterion
CCs = connect(pgraph)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write clusters to file
with open('out/pconnect.txt', 'w') as outfile:
    for i, CC in enumerate(filter(lambda x: len(x) > 2, CCs)):
        CCid = hex(i)[2:].zfill(4)
        outfile.write(CCid + ':' + ','.join(CC) + '\n')

"""
DEPENDENCIES
../make_xreciprocal/make_preciprocal.py
    ../make_xreciprocal/out/pgraph.tsv
./DFS.py
"""