"""Convert hits to an undirected ggraph."""

import multiprocessing as mp
import os
import pandas as pd


def load_hit(qspid, sspid):
    df = pd.read_csv(f'../hsps2hits/out/{qspid}/{sspid}', sep='\t',
                     usecols=dtypes.keys(), dtype=dtypes, memory_map=True)
    r = pd.read_csv(f'../hits2reciprocal/out/{qspid}/{sspid}', sep='\t',
                    usecols=['reciprocal2'], memory_map=True)
    dfr = df.join(r)

    return dfr[dfr['reciprocal2']].groupby(['qgnid', 'sgnid'])['bitscore'].mean()


dtypes = {'qppid': 'string', 'qgnid': 'string',
          'sppid': 'string', 'sgnid': 'string',
          'bitscore': float}
num_processes = 2

if __name__ == '__main__':
    # Load data
    with mp.Pool(processes=num_processes) as pool:
        tsvs = [(qspid, sspid) for qspid in os.listdir('../hsps2hits/out/')
                for sspid in os.listdir(f'../hsps2hits/out/{qspid}/')]
        hits = pd.concat(pool.starmap(load_hit, tsvs))

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
    with open('out/ggraph2.tsv', 'w') as file:
        for qgnid, edges in ggraph.items():
            file.write(qgnid + '\t' + ','.join([sgnid + ':' + str(bitscore) for sgnid, bitscore in edges]) + '\n')

"""
DEPENDENCIES
../hsps2hits/hsps2hits.py
    ../hsps2hits/out/*/*.tsv
../hits2reciprocal/hits2reciprocal.py
    ../hits2reciprocal/out/*/*.tsv
"""