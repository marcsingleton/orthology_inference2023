"""Make consensus trees from bootstrapped IQ-TREE runs."""

import os

import matplotlib.pyplot as plt
import skbio
from src.draw import plot_tree
from src.ortho_tree.tree import majority_consensus


if not os.path.exists('out/'):
    os.mkdir('out/')

for label in os.listdir('../iqtree_WAG/out/'):
    trees = []
    for file in filter(lambda x: x.endswith('.treefile'), os.listdir(f'../iqtree_WAG/out/{label}/')):
        tree = skbio.read(f'../iqtree_WAG/out/{label}/{file}', 'newick', skbio.TreeNode)
        outgroup = tree.find('sleb').ancestors()[0]
        tree = tree.root_at(outgroup)
        trees.append(tree)
    ctree = majority_consensus(trees)
    for node in ctree.traverse():
        node.children = sorted(node.children, key=lambda x: len(list(x.tips())))
    skbio.write(ctree, 'newick', f'out/{label}.txt')

    # Save image as PNG
    fig, ax = plot_tree(ctree, tip_fontsize=8.5)
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('')
    plt.savefig(f'out/{label}.png')
    plt.close()

    # Save image as with supports
    for node in ctree.traverse():
        if node.support == 1:
            node.support = None
    fig, ax = plot_tree(ctree, tip_fontsize=8.5, support_labels=True, support_fontsize=8.5,
                        support_ha='right', support_hoffset=-0.005, support_voffset=0.25)
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('')
    plt.savefig(f'out/{label}_supports.png')
    plt.savefig(f'out/{label}_supports.svg')
    plt.close()

"""
DEPENDENCIES
../../../src/draw.py
../../../src/ortho_tree/tree.py
../iqtree_WAG/iqtree_WAG.py
    ../iqtree_WAG/out/*/meta_*.treefile
"""