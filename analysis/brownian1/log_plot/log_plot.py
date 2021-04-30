"""Plot log-linear and log-log plots of tails."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from numpy import log

# Input variables
path = '../pic_calc/out/pics.tsv'
lt = 32

# Read data and filter
pics = pd.read_csv(path, sep='\t', index_col=list(range(3)))
pics_lt = pics[(pics.index.get_level_values('min_length') >= lt) &
               (~pics.index.get_level_values('ordered').array.astype(bool))]
rates = (pics_lt ** 2).groupby('block_id').mean()

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

for feature in pics_lt:
    data = rates.loc[rates[feature] != 0, feature].sort_values()
    q = [1 - (i - 0.5) / len(data) for i in range(1, len(data) + 1)]

    idx = int(0.5 * len(data))
    plt.scatter(data[idx:], log(q)[idx:], s=2)
    plt.title(f'{feature}: Tail Probability Plot')
    plt.xlabel(feature)
    plt.ylabel('log(tail probability)')
    plt.savefig(f'out/{feature}_loglin.png')
    plt.close()

    plt.scatter(log(data)[idx:], log(q)[idx:], s=2)
    plt.title(f'{feature}: Tail Probability Plot')
    plt.xlabel(f'log({feature})')
    plt.ylabel('log(tail probability)')
    plt.savefig(f'out/{feature}_loglog.png')
    plt.close()

"""
DEPENDENCIES
../pic_calc/pic_calc.py
    ../pic_calc/out/pics.tsv
"""