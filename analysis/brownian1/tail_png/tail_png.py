"""Create alignment png of blocks at tails of feature distributions."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from src.draw import draw_msa

# Input variables
path_pics = '../pic_calc/out/pics.tsv'
path_segs = '../segment_avg/out/segment_avg.tsv'
lt = 32  # Length threshold
idx = 25  # Number of contrasts at each tail

# Read data and filter
pics = pd.read_csv(path_pics, sep='\t', index_col=list(range(3)))
pics_lt = pics[pics.index.get_level_values('min_length') >= lt]
segs = pd.read_csv(path_segs, sep='\t', keep_default_na=False)

for feature in pics:
    # Make output directories for features
    cur_dir = f'out/{feature}/'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)  # Recursive folder creation

    sort = pics_lt[feature].sort_values().index.get_level_values('block_id')
    extremes = set(sort[:idx].unique()) | set(sort[-idx:].unique())

    for extreme in extremes:
        # Extract block MSA from segs
        MSA = segs.loc[segs['block_id'] == extreme, 'seq'].array
        im = draw_msa(MSA)
        plt.imsave(cur_dir + extreme + '.png', im)

"""
DEPENDENCIES
../../../src/draw.py
../pic_calc/pic_calc.py
    ../pic_calc/out/pics.tsv
../segment_avg/segment_avg.py
    ../segment_avg/out/segment_avg.tsv
"""