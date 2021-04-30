"""Extract type information from segmented blocks."""

import os
import pandas as pd

# Read data
path = '../segment_avg/out/segment_avg.tsv'
segs = pd.read_csv(path, sep='\t', keep_default_na=False)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Extract types and save
reduced = segs[['block_id', 'ordered']].drop_duplicates()
reduced.to_csv('out/ordered.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../segment_avg/segment_avg.py
    ../segment_avg/out/segment_avg.tsv
"""