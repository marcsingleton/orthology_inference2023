"""Compute reciprocity for each HSP."""

import os

# Load pgraph
pgraph = {}
with open('../hsps2pgraph/out/pgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        pgraph[node] = adjs.split(',')

# Check reciprocity
rows = []
with open('../blast2hsps/out/hsps.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        fields = line.split()
        qppid, sppid = fields[0], fields[3]
        try:
            r = qppid in pgraph[sppid]
        except KeyError:
            r = False
        rows.append((qppid, sppid, str(r)))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write to file
with open('out/hsps.tsv', 'w') as file:
    file.write('qppid\tsppid\treciprocal\n')
    for row in rows:
        file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../blast2hsps/blast2hsps.py
    ../blast2hsps/hsps.tsv
../hsps2pgraph/hsps2pgraph.py
    ../hsps2pgraph/out/pgraph.tsv
"""