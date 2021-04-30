"""Plot correlations of conscore and the maximum contrast value in a block."""

import matplotlib.pyplot as plt
import os
import pandas as pd


def to_dict(path, dtype):
    dict = {}
    with open(path) as file:
        file.readline()  # Skip first line
        for line in file:
            key, val = line.rstrip().split('\t')
            dict[key] = dtype(val)
    return dict


# Input variables
path_pics = '../pic_calc/out/pics.tsv'
path_scores = '../meta_conscore/out/scores_match.tsv'

# Read data
pics = pd.read_csv(path_pics, sep='\t', index_col=list(range(3)))
scores = to_dict(path_scores, float)

for lt in [2 ** x for x in range(6)]:
    # Filter data and create x, y vectors
    pics_lt = pics[pics.index.get_level_values('min_length') >= lt]
    pics_max = pics.abs().groupby('block_id').max()
    scores_max = [scores[block_id] for block_id in pics_max.index.get_level_values('block_id')]

    # Make output directories for cutoffs
    cur_dir = f'out/{lt}/'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)  # Recursive folder creation

    # Plot feature scatters
    for feature in pics_max:
        plt.scatter(scores_max, pics_max[feature], s=5, edgecolors='none')
        plt.xlabel('Match Score')
        plt.ylabel(feature + ' Absolute Maximum')
        plt.savefig(cur_dir + feature + '.png')
        plt.close()

"""
DEPENDENCIES
../meta_conscore/meta_conscore.py
    ../meta_conscore/out/scores_match.tsv
../pic_calc/pic_calc.py
    ../pic_calc/out/pics.tsv
"""