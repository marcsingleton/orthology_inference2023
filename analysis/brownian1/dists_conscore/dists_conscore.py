"""Plot the conservation score distributions by length."""

import matplotlib.pyplot as plt
import os


def to_dict(path, dtype):
    dict = {}
    with open(path) as file:
        file.readline()  # Skip header
        for line in file:
            key, val = line.rstrip().split('\t')
            dict[key] = dtype(val)
    return dict


# Input variables
path_blosum_score = '../meta_conscore/out/scores_blosum.tsv'
path_match_score = '../meta_conscore/out/scores_match.tsv'
path_block_len = '../meta_block_len/out/block_len.tsv'
path_ordered = '../meta_ordered/out/ordered.tsv'

# Load dictionaries
scores_blosum = to_dict(path_blosum_score, float)
scores_match = to_dict(path_match_score, float)
block_len = to_dict(path_block_len, int)
ordered = to_dict(path_ordered, lambda x: True if x == 'True' else False)

for lt in [2 ** x for x in range(6)]:
    for scores, out_dir, name in [(scores_blosum, 'blosum', 'BLOSUM'), (scores_match, 'match', 'Match')]:
        # Make output directory
        cur_dir = f'out/{out_dir}/'
        if not os.path.exists(cur_dir):
            os.makedirs(cur_dir)  # Recursive folder creation

        # Plot scores as single series
        _, data = zip(*filter(lambda x: block_len[x[0]] >= lt, scores.items()))
        plt.hist(data, bins=100)
        plt.title(f'Average Pairwise {name} Score of Blocks \nwith Minimum Length Greater than {lt}')
        plt.xlabel(f'{name} Score')
        plt.ylabel('Count')
        plt.savefig(cur_dir + f'{lt}.png')
        plt.close()

        # Plot scores as ordered and disordered series
        _, ord = zip(*filter(lambda x: block_len[x[0]] >= lt and ordered[x[0]], scores.items()))
        plt.hist(ord, bins=100, color='C0')
        plt.title(f'Average Pairwise {name} Score of Ordered Blocks \nwith Minimum Length Greater than {lt}')
        plt.xlabel(f'{name} Score')
        plt.ylabel('Count')
        plt.savefig(cur_dir + f'ord{lt}.png')
        plt.close()

        _, dis = zip(*filter(lambda x: block_len[x[0]] >= lt and not ordered[x[0]], scores.items()))
        plt.hist(dis, bins=100, color='C1')
        plt.title(f'Average Pairwise {name} Score of Disordered Blocks \nwith Minimum Length Greater than {lt}')
        plt.xlabel(f'{name} Score')
        plt.ylabel('Count')
        plt.savefig(cur_dir + f'dis{lt}.png')
        plt.close()

"""
DEPENDENCIES
../meta_block_len/meta_block_len.py
    ../meta_block_len/out/block_len.tsv
../meta_conscore/meta_conscore.py
    ../meta_conscore/out/scores_blosum.tsv
    ../meta_conscore/out/scores_match.tsv
../meta_ordered/meta_ordered.py
    ../meta_ordered/out/ordered.tsv
"""