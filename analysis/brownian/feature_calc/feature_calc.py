"""Calculate features for all non-empty sequences."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ

# Input variables
path = '../segment_avg/segment_avg.tsv'  # Path to segmented sequences
type_name = 'ordered'  # Name of column denoting block class
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    # Clean data by removing gaps and empty sequences
    segs = pd.read_csv(path, sep='\t', keep_default_na=False)
    segs['seq'] = segs['seq'].map(lambda x: x.translate({ord('-'): None}))
    segs_lt = segs[segs['seq'].map(lambda x: len(x) >= 1)]  # Select entries where sequence is non-empty

    # Extract indices
    block_ids = segs_lt['block_id']
    species = segs_lt['seq_id'].map(lambda x: x[0:4]).rename('species_id')
    seg_ids = segs_lt['seg_id']
    types = segs_lt[type_name]
    lengths = segs_lt['seq'].map(len).rename('length')

    # Compute features
    with mp.Pool(processes=num_processes) as pool:
        results = pd.DataFrame(pool.imap(seqfeat.feat_all, segs_lt['seq'], chunksize=50))

    # Merge subsets and save
    results.set_index([block_ids, species, seg_ids, types, lengths], inplace=True)
    results.to_csv('features.tsv', sep='\t')

"""
DEPENDENCIES
../segment_avg/segment_avg.py
    ../segment_avg/segment_avg.tsv
"""