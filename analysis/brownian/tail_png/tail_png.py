"""Create alignment png of blocks at tails of feature distributions."""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from math import ceil

# Input variables
path_pics = '../pic_calc/pics.tsv'
path_segs = '../segment_avg/segment_avg.tsv'
lt = 32  # Length threshold
idx = 25  # Number of contrasts at each tail

# PNG parameters
colors_lite = {'A': 'bfbfbf', 'C': 'f2f27f', 'D': 'f28484', 'E': 'f28484', 'F': '9898d4',
               'G': 'dfdfdf', 'H': 'c0c0e8', 'I': '86c086', 'K': '89acff', 'L': '86c086',
               'M': 'f2f27f', 'N': '7feded', 'P': 'edcac0', 'Q': '7feded', 'R': '89acff',
               'S': 'fcca7f', 'T': 'fcca7f', 'V': '86c086', 'W': 'd9acd9', 'Y': '9898d4',
               '-': 'ffffff', 'X': '7f7f7f'}
cols_im = 140  # Number of columns in output
spacing = 25  # Spacing between blocks
sym_length = 7  # Length of a symbol rectangle
sym_height = 7  # Height of a symbol rectangle

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
        # Extract block MSA from segs and calculate image parameters
        MSA = segs.loc[segs['block_id'] == extreme, 'seq'].array
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
                color = [int(hex[i:i + 2], 16) for i in (0, 2, 4)]

                # Fill slice with color
                im[slice(y, y + sym_height), slice(x, x + sym_length), :] = color
                if MSA[row][col] == '-':
                    color = [int('3f3f3f'[i:i + 2], 16) for i in (0, 2, 4)]
                    y1, y2 = y + ceil((sym_height - 2) / 2), y + ceil((sym_height + 1) / 2)
                    x1, x2 = x + ceil((sym_length - 2) / 2), x + ceil((sym_length + 1) / 2)
                    im[slice(y1, y2), slice(x1, x2), :] = color

        plt.imsave(cur_dir + extreme + '.png', im)

"""
DEPENDENCIES
../pic_calc/pic_calc.py
    ../pic_calc/pics.tsv
../segment_avg/segment_avg.py
    ../segment_avg/segment_avg.tsv
"""