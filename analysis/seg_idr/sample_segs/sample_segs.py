"""Create samples sequences at different cutoffs."""

import os
import pandas as pd

# Input variables
path = '../segment_iupred2a/out/segment_iupred2a.tsv'
type_name = 'ordered'  # Name of column denoting segment type

# Constants
n = 7000

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

segs = pd.read_csv(path, sep='\t', keep_default_na=False)
segs['seq'] = segs['seq'].map(lambda x: x.translate({ord('-'): None}))  # Remove gaps if present
for i in [2 ** x for x in range(6)]:
    # Read data and split segments
    segs_lt = segs[segs['seq'].map(lambda x: len(x) >= i)]
    T_segs = segs_lt[segs_lt[type_name]]
    F_segs = segs_lt[~segs_lt[type_name]]

    # Sample segments
    T_sample = T_segs.sample(n).sort_values('seg_id')
    F_sample = F_segs.sample(n).sort_values('seg_id')

    # Merge samples and save
    sample = pd.concat([T_sample, F_sample])
    sample.to_csv(f'out/segments_{i}.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../segment_iupred2a/segment_iupred2a.py
    ../segment_iupred2a/out/segment_iupred2a.tsv
"""