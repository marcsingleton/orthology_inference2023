"""Create alignment png of blocks corresponding to different conservation score quantiles."""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from math import ceil
from random import sample


def to_dict(path, dtype):
    dict = {}
    with open(path) as file:
        file.readline()  # Skip header
        for line in file:
            key, val = line.rstrip().split('\t')
            dict[key] = dtype(val)
    return dict


# Path variables
path_blosum_scores = '../meta_conscore/out/scores_blosum.tsv'
path_match_scores = '../meta_conscore/out/scores_match.tsv'
path_block_len = '../meta_block_len/out/block_len.tsv'
path_ordered = '../meta_ordered/out/ordered.tsv'
path_segs = '../segment_avg/out/segment_avg.tsv'

# Filter and sampling variables
lt = 32
num_partitions = 10
num_samples = 10

# PNG parameters
colors_lite = {'A': 'bfbfbf', 'C': 'f2f27f', 'D': 'f28484', 'E': 'f28484', 'F': '9898d4',
               'G': 'dfdfdf', 'H': 'c0c0e8', 'I': '86c086', 'K': '89acff', 'L': '86c086',
               'M': 'f2f27f', 'N': '7feded', 'P': 'edcac0', 'Q': '7feded', 'R': '89acff',
               'S': 'fcca7f', 'T': 'fcca7f', 'V': '86c086', 'W': 'd9acd9', 'Y': '9898d4',
               '-': 'f6f6f6', 'X': '7f7f7f'}
cols_im = 140  # Number of columns in output
spacing = 25  # Spacing between blocks
sym_length = 7  # Length of a symbol rectangle
sym_height = 7  # Height of a symbol rectangle

# Read data
scores_blosum = to_dict(path_blosum_scores, float)
scores_match = to_dict(path_match_scores, float)
block_len = to_dict(path_block_len, int)
ordered = to_dict(path_ordered, lambda x: True if x == 'True' else False)
segs = pd.read_csv(path_segs, sep='\t', keep_default_na=False)

for scores, out_dir in ((scores_blosum, 'blosum'), (scores_match, 'match')):
    # Make output directory
    cur_dir = f'out/{out_dir}/'
    if not os.path.exists(cur_dir):
        os.makedirs(cur_dir)  # Recursive folder creation

    # Filter scores by length and calculate partition length
    scores_lt = [(key, val) for key, val in scores.items() if block_len[key] >= lt and not ordered[key]]
    min_val = min(scores_lt, key=lambda x: x[1])[1]
    len_partitions = (max(scores_lt, key=lambda x: x[1])[1] - min_val) / num_partitions

    for i in range(num_partitions):
        # Sample scores within partitions
        scores_pt = [(key, val) for key, val in scores_lt if (i + 1) * len_partitions >= val - min_val >= i * len_partitions]
        scores_sample = sample(scores_pt, num_samples)

        for key, val in scores_sample:
            # Extract block MSA from segs and calculate image parameters
            MSA = segs.loc[segs['block_id'] == key, 'seq'].array
            rows_MSA, cols_MSA = len(MSA), len(MSA[0])
            length_im = sym_length * min(cols_MSA, cols_im)  # Length of final image
            blocks_im = cols_MSA // cols_im  # Number of blocks in addition to the first
            height_im = (sym_height * rows_MSA + spacing) * blocks_im + sym_height * rows_MSA  # Height of final image

            # Instantiate array and fill with values
            im = np.full((height_im, length_im, 3), 255, dtype='uint8')
            for row in range(rows_MSA):
                for col in range(cols_MSA):
                    # Position of symbol rectangle
                    block = col // cols_im
                    y = (sym_height * rows_MSA + spacing) * block + sym_height * row
                    x = col % cols_im * sym_length

                    # Create color tuple
                    hex = colors_lite[MSA[row][col]]
                    color = [int(hex[i:i+2], 16) for i in (0, 2, 4)]

                    # Fill slice with color
                    im[slice(y, y + sym_height), slice(x, x + sym_length), :] = color
                    if MSA[row][col] == '-':
                        color = [int('3f3f3f'[i:i + 2], 16) for i in (0, 2, 4)]
                        y1, y2 = y + ceil((sym_height - 2) / 2), y + ceil((sym_height + 1) / 2)
                        x1, x2 = x + ceil((sym_length - 2) / 2), x + ceil((sym_length + 1) / 2)
                        im[slice(y1, y2), slice(x1, x2), :] = color

            val_str = str(val)[:5].translate({ord('.'): '_'})  # Convert decimal point to underscore
            plt.imsave(cur_dir + f'pt{i}_{val_str}_{key}.png', im)

    # Output partition intervals
    print(out_dir.upper())
    for i in range(num_partitions):
        print(min_val + i * len_partitions, min_val + (i + 1) * len_partitions, sep=', ')

"""
OUTPUT
BLOSUM
-1.1465319865319865, -0.2751944845278178
-0.2751944845278178, 0.5961430174763509
0.5961430174763509, 1.4674805194805198
1.4674805194805198, 2.3388180214846885
2.3388180214846885, 3.210155523488857
3.210155523488857, 4.081493025493026
4.081493025493026, 4.952830527497195
4.952830527497195, 5.824168029501363
5.824168029501363, 6.695505531505532
6.695505531505532, 7.5668430335097
MATCH
0.08418079096045197, 0.17576271186440678
0.17576271186440678, 0.2673446327683616
0.2673446327683616, 0.3589265536723164
0.3589265536723164, 0.4505084745762712
0.4505084745762712, 0.542090395480226
0.542090395480226, 0.6336723163841809
0.6336723163841809, 0.7252542372881357
0.7252542372881357, 0.8168361581920904
0.8168361581920904, 0.9084180790960452
0.9084180790960452, 1.0

DEPENDENCIES
../meta_block_len/meta_block_len.py
    ../meta_block_len/out/block_len.tsv
../meta_conscore/meta_conscore.py
    ../meta_conscore/out/scores_blosum.tsv
    ../meta_conscore/out/scores_match.tsv
../meta_ordered/meta_ordered.py
    ../meta_ordered/out/ordered.tsv
../segment_avg/segment_avg.py
    ../segment_avg/out/segment_avg.tsv
"""