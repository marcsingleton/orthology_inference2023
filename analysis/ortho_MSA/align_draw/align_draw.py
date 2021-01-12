"""Draw test alignments."""

import os
import re
from math import ceil

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skbio


def draw_alignment(OGid, MSA):
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
    plt.imsave(f'out/{OGid}.png', im)


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


rows = []
with open('../../ortho_cluster3/clique5+_community/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        for gnid in gnids:
            rows.append({'CCid': CCid, 'OGid': OGid, 'gnid': gnid})
OGs = pd.DataFrame(rows)

test_genes = pd.read_table('test_genes.tsv')
OGid2meta = pd.read_table('../OGid2meta/out/OGid2meta.tsv').drop(['CCid', 'edgenum'], axis=1)
tree = skbio.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick', skbio.TreeNode)
tree = tree.shear([tip.name for tip in tree.tips() if tip.name != 'sleb'])

if not os.path.exists('out/'):
    os.mkdir('out/')

df = OGs.merge(test_genes, on='gnid', how='right').merge(OGid2meta, on='OGid', how='left')
df.to_csv('out/OGids.tsv', sep='\t', index=False)

for row in df.dropna().itertuples():
    OGid = row.OGid
    gnidnum, spidnum, sqidnum = row.gnidnum, row.spidnum, row.sqidnum
    if gnidnum == spidnum == sqidnum == 26:
        MSA = load_alignment(f'../align_fastas1/out/{OGid}.mfa')
    else:
        MSA = load_alignment(f'../align_fastas2-1/out/{OGid}.mfa')

    order = {tip.name: i for i, tip in enumerate(tree.tips())}
    MSA = sorted(MSA, key=lambda x: order[x[0]])  # Re-order sequences
    draw_alignment(OGid, MSA)

"""
DEPENDENCIES
../../ortho_cluster3/clique5+_community/clique5+_community2.py
    ../../ortho_cluster3/clique5+_community/out/ggraph2/5clique/gclusters.txt
../ortho_tree/consensus_tree/consensus_tree.py
    ../ortho_tree/consensus_tree/out/100red_ni.txt
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.mfa
../align_fastas2-1/align_fastas2-2.py
    ../align_fastas2-1/out/*.mfa
../OGid2meta/OGid2meta.py
    ../OGid2meta/out/OGid2meta.tsv
./test_genes.tsv
"""