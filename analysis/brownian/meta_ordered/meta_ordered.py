"""Extract type information from segmented blocks."""

import pandas as pd

# Read data
path = '../segment_avg/segment_avg.tsv'
segs = pd.read_csv(path, sep='\t', keep_default_na=False)

# Extract types and save
reduced = segs[['block_id', 'ordered']].drop_duplicates()
reduced.to_csv('ordered.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../segment_avg/segment_avg.py
    ../segment_avg/segment_avg.tsv
"""