"""Calculate features for all non-empty sequences."""

import multiprocessing as mp
import pandas as pd
import seqfeat
from os import environ

# Input variables
path = '../segment_avg/segment_avg.tsv'  # Path to segmented sequences
key = 'ordered'  # Key of column denoting block class
num_processes = int(environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    # Clean data by removing gaps and empty sequences
    blocks = pd.read_csv(path, sep='\t', keep_default_na=False)
    blocks['seq'] = blocks['seq'].map(lambda x: x.translate({ord('-'): None}))
    blocks_lt = blocks[blocks['seq'].map(lambda x: len(x) >= 1)]  # Select entries where sequence is non-empty

    # Extract indices
    block_ids = blocks_lt['block_id']
    classes = blocks_lt[key]
    species = blocks_lt['seq_id'].map(lambda x: x[0:4])

    # Compute features
    with mp.Pool(processes=num_processes) as pool:
        results = pd.DataFrame(pool.imap(seqfeat.feat_all, blocks_lt['seq'], chunksize=50))

    # Merge subsets and save
    results.set_index([block_ids, classes, species], inplace=True)
    results.to_csv(f'features.tsv', sep='\t')

"""
DEPENDENCIES
../segment_avg/segment_avg.py
    ../segment_avg/segment_avg.tsv
"""