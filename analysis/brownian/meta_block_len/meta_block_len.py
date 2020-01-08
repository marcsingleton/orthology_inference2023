"""Calculate the minimum block length for each block."""

import pandas as pd

# Read data
path = '../segment_avg/segment_avg.tsv'
segs = pd.read_csv(path, sep='\t', keep_default_na=False)

# Create lengths dataframe
lengths = segs['seq'].map(lambda x: len(x.translate({ord('-'): None})))
lengths.name = 'min_length'
lengths.index = segs['block_id']

# Extract lengths and save
min_lengths = lengths.groupby('block_id').min()
min_lengths.to_csv('block_len.tsv', sep='\t', header=True)

"""
DEPENDENCIES
../segment_avg/segment_avg.py
    ../segment_avg/segment_avg.tsv
"""