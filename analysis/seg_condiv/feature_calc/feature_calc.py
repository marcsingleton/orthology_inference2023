"""Calculate features for all non-empty sequences"""

import multiprocessing as mp
import os
import pandas as pd
import src.seg_scripts.seqfeat as seqfeat

# Input variables
path = '../segment_aliscore/out/segment_aliscore.tsv'  # Path to segmented sequences
type_name = 'conserved'  # Name of column denoting segment type
num_processes = int(os.environ['SLURM_NTASKS'])

if __name__ == '__main__':  # Multiprocessing can only occur in top-level script (likely to prevent recursion)
    # Clean data by removing gaps and empty sequences
    segs = pd.read_csv(path, sep='\t', keep_default_na=False)
    segs['seq'] = segs['seq'].map(lambda x: x.translate({ord('-'): None}))
    segs_lt = segs[segs['seq'].map(lambda x: len(x) >= 1)]  # Select entries where sequence is non-empty

    # Extract indices
    seg_ids = segs_lt['seg_id']
    types = segs_lt[type_name]
    lengths = segs_lt['seq'].map(len).rename('length')

    # Compute features
    with mp.Pool(processes=num_processes) as pool:
        features = pd.DataFrame(pool.imap(seqfeat.feat_all, segs_lt['seq'], chunksize=50))

    # Make output directory
    if not os.path.exists('out/'):
        os.mkdir('out/')

    # Set index and save
    features.set_index([seg_ids, types, lengths], inplace=True)
    features.to_csv('out/features.tsv', sep='\t')

"""
DEPENDENCIES
../../../src/seg_scripts/seqfeat.py
../segment_aliscore/segment_aliscore.py
    ../segment_aliscore/out/segment_aliscore.tsv
"""