"""Find connected components of ggraph."""

import os
from DFS import connect

# Load ggraph
ggraph = {}
with open('../hsps2ggraph/out/ggraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        ggraph[node] = adjs.split(',')

# Cluster by triangle criterion
CCs = connect(ggraph)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write clusters to file
with open('out/gconnect.txt', 'w') as outfile:
    for i, CC in enumerate(filter(lambda x: len(x) > 2, CCs)):
        CCid = hex(i)[2:].zfill(4)
        outfile.write(CCid + ':' + ','.join(CC) + '\n')

"""
DEPENDENCIES
../hsps2ggraph/hsps2ggraph.py
    ../hsps2ggraph/out/ggraph.tsv
./DFS.py
"""