"""Plot number of blocks removed by trimming the tails of feature distributions."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from numpy import linspace

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Read data and filter
pics = pd.read_csv('../pic_calc/out/pics.tsv', sep='\t', index_col=list(range(3)))
for lt in [2 ** x for x in range(6)]:
    pics_lt = pics[(pics.index.get_level_values('min_length') >= lt) & (~pics.index.get_level_values('ordered').array.astype(bool))]

    # Initialize counts and and cutoff range
    idxs = []
    num_blocks = []
    cutoffs = linspace(0.00005, 0.01, 150)
    for cutoff in cutoffs:
        idx = int(cutoff * len(pics_lt))  # Adding 1 unnecessary due to slice indexing
        idxs.append(idx)

        # Extract extremes and counts
        extremes = set()
        for feature in pics_lt:
            feature_sort = pics_lt[feature].sort_values()
            lower = set(feature_sort[:idx].index.get_level_values('block_id'))
            upper = set(feature_sort[-idx:].index.get_level_values('block_id'))
            extremes |= lower | upper
        num_blocks.append(len(extremes))

    # Plot cutoff idx
    plt.plot(cutoffs, idxs)
    plt.xlabel('Tail Fraction Removed')
    plt.ylabel('Number of Contrasts Removed at Tail')
    plt.savefig(f'out/{lt}_contrasts.png')
    plt.close()

    # Plot blocks removed
    fig, ax1 = plt.subplots()
    ax1.plot(cutoffs, num_blocks)
    ax1.set_xlabel('Tail Fraction Removed')
    ax1.set_ylabel('Number of Blocks Removed')

    ax2 = ax1.twinx()
    mn, mx = ax1.get_ylim()
    n = len(pics_lt.index.get_level_values('block_id').unique())
    ax2.set_ylim(mn / n, mx / n)
    ax2.set_ylabel('Fraction of Blocks Removed', labelpad=10)
    fig.subplots_adjust(left=0.125, right=0.875)
    plt.savefig(f'out/{lt}_blocks.png')
    plt.close()

"""
NOTES
The number of blocks removed increases surprisingly rapidly with the cutoff threshold.
    A cutoff of 0.0005 always removes over 5% of the blocks.
My expectation was that the most extreme contrasts across all features generally correspond to the same contrasts or blocks.
There are two factors that can account for the rapid increase.
    The features are generally independent, so an extreme value in one feature does not broadly correlate with extreme values in other features.
    The number of contrasts is so large, that removing a small fraction will correspond to a large number of contrasts.
        Given that the extreme values are generally uncorrelated, the number of unique blocks will increase rapidly.
            This has some similarities to the birthday problem.
A cutoff of 0.0002 balances removing the extremes without inflating to block removal fraction to unacceptable levels.
    8 contrasts are removed from the left and right tails. This is consistent with the conscore_corr plots showing the number of extreme values is on the order of 10.
    208 blocks are removed from the data set. This corresponds to 4.5% of the blocks.
        Removing the most conserved blocks will further inflate this rate.

DEPENDENCIES
../pic_calc/pic_calc.py
    ../pic_calc/out/pics.tsv
"""