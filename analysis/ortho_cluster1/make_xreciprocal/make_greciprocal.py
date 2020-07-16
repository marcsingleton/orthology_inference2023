"""Remove non-reciprocal hits between polypeptides in gene graph."""

import json
import os

# Load best hits graph
with open('../blast2ggraph/out/ggraph.json') as file:
    iggraph = json.load(file)

# Copy graph, including only reciprocal hits
oggraph = {}
for qgnid, sgnids in iggraph.items():
    if 'null' in sgnids:
        del sgnids['null']  # Remove None first to not loop over

    for sgnid, qppids in sgnids.items():
        for qppid, sppids in qppids.items():
            for sppid, params in sppids.items():
                try:
                    reciprocal = qppid in iggraph[sgnid][qgnid][sppid]
                except KeyError:
                    reciprocal = False

                if reciprocal:
                    try:
                        oggraph[qgnid][sgnid][qppid][sppid] = params
                    except KeyError as err:
                        if err.args[0] == qgnid:
                            oggraph[qgnid] = {sgnid: {qppid: {sppid: params}}}
                        elif err.args[0] == sgnid:
                            oggraph[qgnid][sgnid] = {qppid: {sppid: params}}
                        elif err.args[0] == qppid:
                            oggraph[qgnid][sgnid][qppid] = {sppid: params}

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write graph to file
with open('out/ggraph.json', 'w') as outfile:
    json.dump(oggraph, outfile, indent=1)

# Write graph as adjacency list to file
with open('out/ggraph.tsv', 'w') as outfile:
    for qgnid, sgnids in oggraph.items():
        outfile.write(qgnid + '\t' + ','.join(sgnids.keys()) + '\n')

"""
DEPENDENCIES
../blast2ggraph/blast2ggraph.py
    ../blast2ggraph/out/ggraph.json
"""