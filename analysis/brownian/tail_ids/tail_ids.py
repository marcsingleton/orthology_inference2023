"""Plot source of contrasts at feature distribution tails."""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from math import ceil

# Input variables
path = '../pic_calc/out/pics.tsv'
num_contrasts = 9
max_percentile = 0.20
lt = 32  # Length threshold

# Read data and filter
pics = pd.read_csv(path, sep='\t', index_col=list(range(3)))
pics_lt = pics[pics.index.get_level_values('min_length') >= lt]
pics_lt.insert(0, 'pic_id', list(range(num_contrasts)) * (len(pics_lt) // num_contrasts))
pics_lt.set_index('pic_id', append=True, inplace=True)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Sort and count representation in tails
size = ceil(0.005 * len(pics_lt))
max_idx = int(max_percentile * len(pics_lt)) - 1  # Subtract 1 due to 0 indexing
x = [i / len(pics_lt) for i in range(size + 1, max_idx - size + 1)]
for feature in pics_lt:
    sort = pics_lt[feature].sort_values()
    for pic_id in range(num_contrasts):
        lcount = sort[:max_idx].index.get_level_values('pic_id') == pic_id
        rcount = sort[:-max_idx - 1: -1].index.get_level_values('pic_id') == pic_id
        tail = np.convolve(lcount, np.ones(2 * size + 1), 'valid') + np.convolve(rcount, np.ones(2 * size + 1), 'valid')
        plt.plot(x, tail, label=pic_id, linewidth=1)

    plt.title(f'{feature}:\nContrast Counts in 1% Sliding Windows')
    plt.xlabel('Left and Right Tail Percentiles')
    plt.ylabel('Count')
    plt.legend(title='Contrast ID', bbox_to_anchor=(1.025, 0.5), loc='center left')
    plt.subplots_adjust(right=0.8)
    plt.savefig(f'out/{feature}.png')
    plt.close()

"""
DEPENDENCIES
../pic_calc/pic_calc.py
    ../pic_calc/out/pics.tsv

NOTES
The counts are taken in non-overlapping windows rather than cumulatively to accurately represent each fraction of the tail.
"""
