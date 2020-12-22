"""Convert hits to an undirected ggraph."""

import multiprocessing as mp
import os
import pandas as pd
from itertools import permutations


def load_hit(qspid, sspid):
    df = pd.read_csv(f'../../ortho_search/hsps2hits/out/{qspid}/{sspid}.tsv', sep='\t',
                     usecols=dtypes.keys(), dtype=dtypes, memory_map=True)
    r = pd.read_csv(f'../../ortho_search/hits2reciprocal/out/{qspid}/{sspid}.tsv', sep='\t',
                    usecols=['reciprocal2'], memory_map=True)
    dfr = df.join(r)

    return dfr[dfr['reciprocal2']].groupby(['qgnid', 'sgnid'])['bitscore'].mean()


dtypes = {'qppid': 'string', 'qgnid': 'string',
          'sppid': 'string', 'sgnid': 'string',
          'bitscore': float}
num_processes = 2

if __name__ == '__main__':
    # Parse parameters
    spids = []
    with open('params.tsv') as file:
        fields = file.readline().split()  # Skip header
        for line in file:
            spids.append(line.split()[0])

    # Load data
    with mp.Pool(processes=num_processes) as pool:
        hits = pd.concat(pool.starmap(load_hit, permutations(spids, 2)))

    ggraph = {}
    for (qgnid, sgnid), bitscore in hits.iteritems():
        try:
            ggraph[qgnid].append((sgnid, round(bitscore)))
        except KeyError:
            ggraph[qgnid] = [(sgnid, round(bitscore))]

    # Make output directory
    if not os.path.exists('out/'):
        os.mkdir('out/')

    # Write to file
    with open('out/ggraph.tsv', 'w') as file:
        for qgnid, edges in ggraph.items():
            file.write(qgnid + '\t' + ','.join([sgnid + ':' + str(bitscore) for sgnid, bitscore in edges]) + '\n')

"""
DEPENDENCIES
../../ortho_search/hsps2hits/hsps2hits.py
    ../../ortho_search/hsps2hits/out/*/*.tsv
../../ortho_search/hits2reciprocal/hits2reciprocal.py
    ../../ortho_search/hits2reciprocal/out/*/*.tsv
./params.tsv
"""