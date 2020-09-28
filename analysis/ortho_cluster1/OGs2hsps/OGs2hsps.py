"""Extract HSPs which are contained in OGs."""

import os
import pandas as pd


def get_ids(qgnid, sgnid, edge2ids):
    try:
        return edge2ids[(qgnid, sgnid)]
    except KeyError:
        try:
            return edge2ids[(sgnid, qgnid)]
        except KeyError:
            return 'NA', 'NA'


dtypes = {'qgnid': 'string', 'sgnid': 'string'}

# Load OG edges
edge2ids = {}
with open('../subcluster_xgraph/out/ggraph/gclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        edge2ids.update({tuple(edge.split(',')): (CCid, OGid) for edge in edges.split('\t')})

# Iterate through hits
df = pd.read_csv('../blast2hsps/out/hsps.tsv', sep='\t',
                 usecols=dtypes.keys(), dtype=dtypes, memory_map=True)
r = pd.read_csv('../hsps2reciprocal/out/hsps.tsv', sep='\t',
                usecols=['reciprocal'], memory_map=True)
dfr = df.join(r)

orows = []
gnids = set()
for irow in dfr.itertuples():
    qgnid, sgnid = irow.qgnid, irow.sgnid
    if (qgnid, sgnid) in gnids:
        continue

    CCid, OGid = get_ids(qgnid, sgnid, edge2ids)

    orows.append((qgnid, sgnid, CCid, OGid))
    gnids.add((qgnid, sgnid))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write to file
with open('out/hsps.tsv', 'w') as file:
    file.write('qgnid\tsgnid\tCCid\tOGid\n')
    for orow in orows:
        file.write('\t'.join(orow) + '\n')

"""
DEPENDENCIES
../blast2hsps/blast2hsps.py
    ../blast2hsps/hsps.tsv
../hsps2reciprocal/hsps2reciprocal.py
    ../hsps2reciprocal/hsps.tsv
../subcluster_xgraph/subcluster_ggraph.py
    ../subcluster_xgraph/out/ggraph/gclusters.txt
"""