"""Draw alignments with largest gap contrasts."""

import os
import re
from math import ceil

import Bio.Phylo as Phylo
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def draw_alignment(path, MSA):
    # Parameters
    ratio = 2.5  # Aspect ratio (length:height)
    spacing = 25  # Spacing between blocks
    sym_length = 7  # Length of a symbol rectangle
    sym_height = 7  # Height of a symbol rectangle
    aa2color = {'A': '6dd7a1', 'I': '55c08c', 'L': '55c08c', 'V': '55c08c', 'M': '55c08c',
                'F': 'b897ec', 'Y': 'b897ec', 'W': 'a180d2',
                'S': 'ffbe74', 'T': 'ffbe74',
                'N': '77eaf4', 'Q': '77eaf4',
                'D': 'ee8485', 'E': 'ee8485',
                'H': '96c4ff', 'K': '7fadea', 'R': '7fadea',
                'C': 'faed70', 'G': 'e2dedd', 'P': 'ffb1f1',
                'X': '93908f', '-': 'ffffff'}

    # Calculated parameters
    rows_MSA, cols_MSA = len(MSA), len(MSA[0][1])

    def get_dims(cols_im):
        length_im = sym_length * cols_im  # Length of final image
        blocks_im = cols_MSA // cols_im - (1 if cols_MSA % cols_im == 0 else 0)  # Number of blocks in addition to the first
        height_im = (sym_height * rows_MSA + spacing) * blocks_im + sym_height * rows_MSA  # Height of final image
        return length_im, height_im

    def get_aspect(cols_im):
        length_im, height_im = get_dims(cols_im)
        return length_im / height_im

    # Get interval
    interval = (1, cols_MSA)
    while interval[1] - interval[0] > 1:
        i1 = (interval[0], (interval[0] + interval[1]) // 2)
        i2 = ((interval[0] + interval[1]) // 2, interval[1])
        if (get_aspect(i1[0]) - ratio) * (get_aspect(i1[1]) - ratio) < 0:
            interval = i1
        elif (get_aspect(i2[0]) - ratio) * (get_aspect(i2[1]) - ratio) < 0:
            interval = i2
        else:
            break
    cols_im = min(interval, key=lambda x: abs(get_aspect(x) - ratio))
    if cols_MSA % cols_im < 0.5 * cols_im:  # Guarantees at least two blocks
        blocks_im = cols_MSA // cols_im  # Total blocks minus 1
        cols_im += ceil((cols_MSA % cols_im) / blocks_im)  # Distribute excess to other blocks

    # Instantiate array and fill with values
    length_im, height_im = get_dims(cols_im)
    im = np.full((height_im, length_im, 3), 255, dtype='uint8')
    for i, record in enumerate(MSA):
        for j, sym in enumerate(record[1]):
            # Position of symbol rectangle
            block = j // cols_im
            y = (sym_height * rows_MSA + spacing) * block + sym_height * i
            x = j % cols_im * sym_length

            # Create color tuple
            hex = aa2color[sym]
            color = [int(hex[i:i + 2], 16) for i in (0, 2, 4)]

            # Fill slice with color
            im[slice(y, y + sym_height), slice(x, x + sym_length), :] = color
            if sym == '-':
                color = [int('3f3f3f'[i:i + 2], 16) for i in (0, 2, 4)]
                y1, y2 = y + ceil((sym_height - 2) / 2), y + ceil((sym_height + 1) / 2)
                x1, x2 = x + ceil((sym_length - 2) / 2), x + ceil((sym_length + 1) / 2)
                im[slice(y1, y2), slice(x1, x2), :] = color
    plt.imsave(path, im)


def load_alignment(path):
    MSA = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                spid = re.search(r'spid=([a-z]+)', line).group(1)
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            MSA.append((spid, seq))
    return MSA


tree = Phylo.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick')
tree.prune('sleb')

OGid2meta = pd.read_table('../OGid2meta/out/OGid2meta.tsv').drop(['CCid', 'edgenum'], axis=1)
df = pd.read_table('out/row_sums.tsv').merge(OGid2meta, on='OGid', how='left')
df['c_sum'] = df[[f'c{i}' for i in range(25)]].sum(axis=1)
df['c_avg'] = df['c_sum'] / df['len2']

# 1 PLOT STATISTICS
# 1.1 Tail fraction
num_contrasts = 25
max_percentile = 0.30

size = ceil(0.015 * len(df))
max_idx = int(max_percentile * len(df)) - 1  # Subtract 1 due to 0 indexing
x = [i / len(df) for i in range(size + 1, max_idx - size + 1)]

contrasts = pd.concat([df[f'c{i}'] for i in range(num_contrasts)], keys=[f'c{i}' for i in range(num_contrasts)],
                      names=['contrast_id', 'OGid']).sort_values(ascending=False)
counts = [(i, contrasts[:max_idx].index.get_level_values('contrast_id') == f'c{i}') for i in range(num_contrasts)]
tails = sorted([(i, np.convolve(count, np.ones(2 * size + 1) / (2 * size + 1), 'valid')) for i, count in counts],
               key=lambda y: sum(y[1]), reverse=True)
for i, tail in tails[:9]:
    plt.plot(x, tail, label=f'c{i}', linewidth=1)
plt.plot(x, sum([tail for _, tail in tails[9:]]), label='others', linewidth=1)
plt.legend()
plt.title('Contrast Fraction in 3% Sliding Windows')
plt.xlabel('Right Tail Percentile')
plt.ylabel('Fraction')
plt.legend(title='Contrast ID', bbox_to_anchor=(1.025, 0.5), loc='center left')
plt.subplots_adjust(right=0.8)
plt.savefig('out/line_contrast_window.png')
plt.close()

# 1.2 Overall contrast distribution
plt.hist(contrasts, bins=200)
plt.xlabel('Contrast Value')
plt.ylabel('Count')
plt.savefig('out/hist_contrast.png')
plt.yscale('log')
plt.savefig('out/hist_contrast_log.png')
plt.close()

# 1.3 Distribution of averages within OGs
plt.hist(contrasts.groupby('OGid').mean(), bins=200)
plt.xlabel('Contrast mean in OG')
plt.ylabel('Count')
plt.savefig('out/hist_contrast_OGmean.png')
plt.yscale('log')
plt.savefig('out/hist_contrast_OGmean_log.png')
plt.close()

# 1.4 Contrast averages across all OGs
plt.bar(list(range(num_contrasts)), [df[f'c{i}'].mean() for i in range(num_contrasts)],
        yerr=[df[f'c{i}'].std()/50 for i in range(num_contrasts)])
plt.xlabel('Contrast ID')
plt.ylabel('Mean ± STD/50')
plt.savefig('out/bar_contrast_mean.png')
plt.close()

# 2 DRAW ALIGNMENTS
if not os.path.exists('out/sum/'):
    os.mkdir('out/sum/')

# 2.1 Ranked by sum
head1 = df.sort_values(by='c_sum', ascending=False).head(100)
for i, row in enumerate((head1.itertuples())):
    if row.sqidnum == 26:
        MSA = load_alignment(f'../align_fastas1/out/{row.OGid}.mfa')
    else:
        MSA = load_alignment(f'../align_fastas2-2/out/{row.OGid}.mfa')

    order = {terminal.name: i for i, terminal in enumerate(tree.get_terminals())}
    MSA = sorted(MSA, key=lambda x: order[x[0]])  # Re-order sequences
    draw_alignment(f'out/sum/{i}_{row.OGid}.png', MSA)

# 2.2 Ranked by avg
if not os.path.exists('out/avg/'):
    os.mkdir('out/avg/')

head1 = df.sort_values(by='c_avg', ascending=False).head(100)
for i, row in enumerate((head1.itertuples())):
    if row.sqidnum == 26:
        MSA = load_alignment(f'../align_fastas1/out/{row.OGid}.mfa')
    else:
        MSA = load_alignment(f'../align_fastas2-2/out/{row.OGid}.mfa')

    order = {terminal.name: i for i, terminal in enumerate(tree.get_terminals())}
    MSA = sorted(MSA, key=lambda x: order[x[0]])  # Re-order sequences
    draw_alignment(f'out/avg/{i}_{row.OGid}.png', MSA)

"""
../../ortho_tree/consensus_tree/consensus_tree.py
    ../../ortho_tree/consensus_tree/out/100red_ni.txt
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-2/align_fastas2-2.py
    ../align_fastas2-2/out/*.mfa
../OGid2meta/OGid2meta.py
    ../OGid2meta/out/OGid2meta.tsv
./align_contrasts_calc.py
    ./out/row_sums.tsv
"""