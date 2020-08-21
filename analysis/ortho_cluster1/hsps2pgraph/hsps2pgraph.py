"""Convert HSPs to a directed pgraph."""

import os

# Make pgraph
pgraph = {}
with open('../blast2hsps/out/hsps.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        fields = line.split()
        qppid, sppid = fields[0], fields[3]
        try:
            pgraph[qppid].add(sppid)
        except KeyError:
            pgraph[qppid] = set([sppid])

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write to file
with open('out/pgraph.tsv', 'w') as file:
    for qppid, sppids in pgraph.items():
        file.write(qppid + '\t' + ','.join(sppids) + '\n')

"""
DEPENDENCIES
../blast2hsps/blast2hsps.py
    ../blast2hsps/hsps.tsv
"""