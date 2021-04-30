"""Distribution of pics where highly conserved blocks are removed."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from math import ceil


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
type_name = 'ordered'
name_T = 'Ordered'
name_F = 'Disordered'
cutoff_ext = 0.01
cutoff_conscore = 0.8

# Read data and filter
scores = to_dict(path_scores, float)
pics = pd.read_csv(path_pics, sep='\t', index_col=list(range(3)))
for lt in [2 ** x for x in range(6)]:
    idx1 = pics.index.get_level_values('min_length') >= lt
    idx2 = pics.index.get_level_values('block_id').map(lambda x: scores[x] <= cutoff_conscore).array.astype(bool)
    pics_lt = pics[idx1 & idx2]

    # Make output directories for length thresholds
    cur_dir = f'out/{lt}/'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)  # Recursive folder creation

    # Plot distributions
    idx = ceil(cutoff_ext * len(pics_lt))  # Ceil ensures the cutoff fraction is removed
    for feature in pics:
        trimmed = pics_lt[feature].sort_values()[idx:-idx]

        # Get indices for plotting
        idx_T = trimmed.index.get_level_values(type_name).array.astype(bool)
        idx_F = ~idx_T

        # Plot classes as one series
        plt.figure()
        _, bins, _ = plt.hist(trimmed, bins=100)
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.savefig(cur_dir + f'{feature}_joint.png')
        plt.close()

        # Plot classes as two series
        plt.figure()
        plt.hist(trimmed[idx_T], label=name_T, bins=bins, alpha=0.5)
        plt.hist(trimmed[idx_F], label=name_F, bins=bins, alpha=0.5)
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.legend()
        plt.savefig(cur_dir + f'{feature}_separate.png')
        plt.close()

        # Plot classes as separate graphs
        plt.figure()
        plt.hist(trimmed[idx_T], label=name_T, bins=100, color='C0')
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.legend()
        plt.savefig(cur_dir + f'{feature}_{name_T}.png')
        plt.close()

        plt.figure()
        plt.hist(trimmed[idx_F], label=name_F, bins=100, color='C1')
        plt.xlabel(feature)
        plt.ylabel('Count')
        plt.legend()
        plt.savefig(cur_dir + f'{feature}_{name_F}.png')
        plt.close()

"""
DEPENDENCIES
../meta_conscore/meta_conscore.py
    ../meta_conscore/out/scores_match.tsv
../pic_calc/pic_calc.py
    ../pic_calc/out/pics.tsv
"""