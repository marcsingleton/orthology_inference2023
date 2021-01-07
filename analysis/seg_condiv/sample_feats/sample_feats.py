"""Extract features corresponding to sampled segments."""

import os
import pandas as pd
import re

# Input variables
num_idx = 3  # Number of index columns
type_name = 'conserved'  # Name of column denoting segment type

# Constants
segment_dir = '../sample_segs/out/'
features_path = '../feature_calc/out/features.tsv'

features = pd.read_csv(features_path, sep='\t', index_col=list(range(num_idx)))

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

paths = filter(lambda x: re.match('segments_[0-9]+\.tsv', x), os.listdir(segment_dir))
for path in paths:
    # Read data
    segs = pd.read_csv(segment_dir + path, sep='\t', keep_default_na=False)

    # Get file index
    j0 = path.find('_')
    j1 = path.find('.tsv')
    i = path[j0 + 1:j1]

    # Extract features and save
    num_id = features.index.names.index('seg_id')
    idx = tuple([slice(None) if i != num_id else segs['seg_id'] for i in range(num_idx)])

    sample = features.loc[idx, :].sort_index(level=[type_name, 'seg_id'], ascending=[False, True])
    sample.to_csv(f'out/features_{i}.tsv', sep='\t')

"""
DEPENDENCIES
../feature_calc/feature_calc.py
    ../feature_calc/out/features.tsv
../sample_segs/sample_segs.py
    ../sample_segs/out/segments_*.tsv
"""