"""Convert HSPs to an undirected ggraph."""

import os
import pandas as pd

dtypes = {'qppid': 'string', 'qgnid': 'string',
          'sppid': 'string', 'sgnid': 'string'}
df = pd.read_csv(f'../blast2hsps/out/hsps.tsv', sep='\t',
                 usecols=dtypes.keys(), dtype=dtypes, memory_map=True)
r = pd.read_csv(f'../hsps2reciprocal/out/hsps.tsv', sep='\t',
                usecols=['reciprocal'], memory_map=True)
hsps = df.join(r)

ggraph = {}
for row in hsps.itertuples():
    qgnid, sgnid = row.qgnid, row.sgnid
    r = row.reciprocal

    if r:
        try:
            ggraph[qgnid].add(sgnid)
        except KeyError:
            ggraph[qgnid] = {sgnid}

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write to file
with open('out/ggraph.tsv', 'w') as file:
    for qgnid, sgnids in ggraph.items():
        file.write(qgnid + '\t' + ','.join(sgnids) + '\n')

"""
DEPENDENCIES
../hsps2hits/hsps2hits.py
    ../hsps2hits/out/*/*.tsv
../hits2reciprocal/hits2reciprocal.py
    ../hits2reciprocal/out/*/*.tsv
"""