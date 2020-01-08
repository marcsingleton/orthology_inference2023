"""Plot scatters of mean, iqr, and std for block feature distributions."""

import matplotlib.pyplot as plt
import os
import pandas as pd

path = '../pic_calc/pics.tsv'
lt = 32
params = {'s': 2.5, 'color': 'black', 'edgecolors': 'none', 'alpha': '0.25'}

pics = pd.read_csv(path, sep='\t', index_col=list(range(3)))
pics_lt = pics[(pics.index.get_level_values('min_length') >= lt) &
               (~pics.index.get_level_values('ordered').array.astype(bool))]

for feature in pics_lt:
    # Make output directories for length thresholds
    cur_dir = f'out/{feature}/'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)  # Recursive folder creation

    data = pics_lt[feature].sort_values()
    means = data.groupby('block_id').mean()
    stds = data.groupby('block_id').std()
    q75s = data.groupby('block_id').quantile(0.75)
    q25s = data.groupby('block_id').quantile(0.25)
    iqrs = q75s - q25s

    plt.figure()
    plt.scatter(means, stds, **params)
    plt.title(f'{feature}:\nMean vs Standard Deviation of Contrasts Within Blocks')
    plt.xlabel('Mean')
    plt.ylabel('Standard Deviation')
    plt.savefig(cur_dir + 'corr_mean_std.png')

    plt.figure()
    plt.scatter(means, iqrs, **params)
    plt.title(f'{feature}:\nMean vs IQR of of Contrasts Within Blocks')
    plt.xlabel('Mean')
    plt.ylabel('IQR')
    plt.savefig(cur_dir + 'corr_mean_iqr.png')

    plt.figure()
    plt.scatter(stds, iqrs, **params)
    plt.title(f'{feature}:\nStandard Deviation vs IQR of Contrasts Within Blocks')
    plt.xlabel('Standard Deviation')
    plt.ylabel('IQR')
    plt.savefig(cur_dir + 'corr_std_iqr.png')

    plt.close('all')

"""
DEPENDENCIES
../pic_calc/pic_calc.py
    ../pic_calc/pics.tsv
"""