"""Calculate features for different length cutoffs."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ
from sys import argv

# Input variables
path = argv[1]  # Path to segmented sequences
type_name = argv[2]  # Name of column denoting segment type
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    segs = pd.read_csv(path, sep='\t', keep_default_na=False)
    segs['seq'] = segs['seq'].map(lambda x: x.translate({ord('-'): None}))
    subs_lt = segs[segs['seq'].map(lambda x: len(x) >= 1)]  # Select entries where sequence is non-empty

    # Extract indices
    seg_ids = subs_lt['seg_id']
    types = subs_lt[type_name]
    lengths = subs_lt['seq'].map(len).rename('length')

    # Compute features
    with mp.Pool(processes=num_processes) as pool:
        results = pd.DataFrame(pool.imap(seqfeat.feat_all, subs_lt['seq'], chunksize=50))

    # Merge subsets and save
    results.set_index([seg_ids, types, lengths], inplace=True)
    results.to_csv('features.tsv', sep='\t')
