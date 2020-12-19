"""Compute reciprocity for each HSP."""

import os
from itertools import groupby


def line2key(line):
    fields = line.split()
    return fields[0], fields[3]


# Load pgraph
pgraph = {}
with open('../hsps2pgraph/out/pgraph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        pgraph[node] = adjs.split(',')

# Check reciprocity
for qspid in os.listdir('../blast2hsps/out/hsps/'):
    for sspid in os.listdir(f'../blast2hsps/out/hsps/{qspid}/'):
        rows = []
        with open(f'../blast2hsps/out/hsps/{qspid}/{sspid}') as file:
            file.readline()  # Skip header
            for key, group in groupby(file, key=line2key):
                qppid, sppid = key[0], key[1]
                try:
                    r = qppid in pgraph[sppid]
                except KeyError:
                    r = False
                rows.extend([(qppid, sppid, str(r))] * len(list(group)))

        # Make output directory
        if not os.path.exists(f'out/{qspid}/'):
            os.makedirs(f'out/{qspid}/')  # Recursive folder creation

        # Write to file
        with open(f'out/{qspid}/{sspid}', 'w') as file:
            file.write('qppid\tsppid\treciprocal\n')
            for row in rows:
                file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../blast2hsps/blast2hsps.py
    ../blast2hsps/hsps/out/*/*.tsv
../hsps2pgraph/hsps2pgraph.py
    ../hsps2pgraph/out/pgraph.tsv
"""