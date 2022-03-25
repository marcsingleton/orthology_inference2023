"""Make consensus trees from bootstrapped IQ-TREE runs."""

import os

import matplotlib.pyplot as plt
import skbio
from src.draw import plot_tree
from src.ortho_tree.tree import majority_consensus


if not os.path.exists('out/'):
    os.mkdir('out/')

for label in os.listdir('../iqtree_LG/out/'):
    trees = []
    for file in filter(lambda x: x.endswith('.treefile'), os.listdir(f'../iqtree_LG/out/{label}/')):
        tree = skbio.read(f'../iqtree_LG/out/{label}/{file}', 'newick', skbio.TreeNode)
        outgroup = tree.find('sleb').ancestors()[0]
        tree = tree.root_at(outgroup)
        trees.append(tree)
    consensus_tree = majority_consensus(trees)
    for node in consensus_tree.traverse():
        node.children = sorted(node.children, key=lambda x: len(list(x.tips())))
    skbio.write(consensus_tree, 'newick', f'out/{label}.nwk')

    # Save image as PNG
    fig, ax = plot_tree(consensus_tree, tip_fontsize=8.5)
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('')
    plt.savefig(f'out/{label}.png')
    plt.close()

    # Save image as with supports
    for node in consensus_tree.traverse():
        if node.support == 1:
            node.support = None
    fig, ax = plot_tree(consensus_tree, tip_fontsize=8.5, support_labels=True, support_fontsize=8.5,
                        support_ha='right', support_hoffset=-0.005, support_voffset=0.25)
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('')
    plt.savefig(f'out/{label}_supports.png')
    plt.savefig(f'out/{label}_supports.svg')
    plt.close()

"""
DEPENDENCIES
../iqtree_LG/iqtree_LG.sh
    ../iqtree_WAG/out/*/meta_*.treefile
"""