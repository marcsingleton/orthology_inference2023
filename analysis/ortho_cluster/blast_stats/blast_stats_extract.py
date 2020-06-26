"""Convert BLAST best hit graph into tabular format."""

import json
import os
import pandas as pd

# Load gn metadata
gnid2spid = {}
with open('../ppid2meta/out/ppid2meta.tsv') as infile:
    for line in infile:
        _, gnid, spid = line.split()
        gnid2spid[gnid] = spid

# Load best hits graph
with open('../blast2ggraph/out/ggraph.json') as file:
    ggraph = json.load(file)

# Extract data from graph
rows = []
for qgnid, sgnids in ggraph.items():
    # Remove None first to not loop over
    if 'null' in sgnids:
        del sgnids['null']
    qspid = gnid2spid[qgnid]

    for sgnid, qppids in sgnids.items():
        try:
            reciprocal = qgnid in ggraph[sgnid]
        except KeyError:
            reciprocal = False
        sspid = gnid2spid[sgnid]

        for qppid, sppids in qppids.items():
            for sppid, params in sppids.items():
                row = {'qppid': qppid, 'qgnid': qgnid, 'qspid': qspid,
                       'sppid': sppid, 'sgnid': sgnid, 'sspid': sspid,
                       'reciprocal': reciprocal, **params}
                rows.append(row)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Save as tsv
df = pd.DataFrame(rows)  # Use df to ensure columns are aligned
df.to_csv('out/ggraph.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../blast2ggraph/blast2ggraph.py
    ../blast2ggraph/out/ggraph.json
../ppid2meta/ppid2meta.py
    ../ppid2meta/out/ppid2meta.tsv
"""