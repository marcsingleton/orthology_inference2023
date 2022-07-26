"""Make consensus trees from bootstrapped IQ-TREE runs."""

import os

import matplotlib.pyplot as plt
import skbio
from src.draw import plot_tree
from src.ortho_tree.tree import majority_consensus

if not os.path.exists('out/'):
    os.mkdir('out/')

for label in [path for path in os.listdir('../iqtree_GTR/out/') if os.path.isdir(f'../iqtree_GTR/out/{path}')]:
    trees = []
    for file in [path for path in os.listdir(f'../iqtree_GTR/out/{label}/') if path.endswith('.treefile')]:
        tree = skbio.read(f'../iqtree_GTR/out/{label}/{file}', 'newick', skbio.TreeNode)
        outgroup = tree.find('sleb').ancestors()[0]
        tree = tree.root_at(outgroup)
        trees.append(tree)
    consensus_tree = majority_consensus(trees)
    for node in consensus_tree.postorder():
        if node.is_tip():
            node.sort_name = node.name
        else:
            node.children = sorted(node.children, key=lambda x: (len(list(x.tips())), x.sort_name))
            node.sort_name = ''.join([child.sort_name for child in node.children])
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
            node.support = None  # Clear supports of 1 to simplify plot
    fig, ax = plot_tree(consensus_tree, tip_fontsize=8.5, support_labels=True,
                        support_format_spec='.2f', support_fontsize=8.5,
                        support_ha='right', support_hoffset=-0.005, support_voffset=-0.005)
    ax.yaxis.set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('')
    plt.savefig(f'out/{label}_supports.png')
    plt.savefig(f'out/{label}_supports.svg')
    plt.close()

"""
DEPENDENCIES
../iqtree_GTR/iqtree_GTR.sh
    ../iqtree_GTR/out/*/meta_*.treefile
"""