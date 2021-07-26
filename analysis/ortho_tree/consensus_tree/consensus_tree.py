"""Make consensus trees from bootstrapped IQ-TREE runs."""

import os

import matplotlib.pyplot as plt
import skbio
from src.draw import plot_tree
from src.ortho_tree.tree import majority_consensus


if not os.path.exists('out/'):
    os.mkdir('out/')

for label in os.listdir('../phyml_GTR/out/'):
    trees = []
    for file in filter(lambda x: x.endswith('.phy_phyml_tree.txt'), os.listdir(f'../phyml_GTR/out/{label}/')):
        tree = skbio.read(f'../phyml_GTR/out/{label}/{file}', 'newick', skbio.TreeNode)
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
NOTES
When a taxon is designated as an outgroup, this algorithm is equivalent to tree popping described in the Mathematics of
Phylogenetics. Since a clade is the subset of a split that does not contain the root, all the clades obtained from the
find_clades method are equivalent to a split where the clade is one subset and the other subset are the remaining taxa
including the root. For concreteness, color the nodes in the clade blue and those not in the clade (necessarily
including the root) as red.

Now we show that searching through the consensus taxa is equivalent to finding the minimum spanning tree of the nodes
colored by each clade. The tree induced by a node whose terminals are a superset of the clade clearly spans that clade.
Thus, by iterating through all the nodes in the consensus tree sorted by decreasing size, the last node that meets this
criterion induces the minimal spanning tree (MST) for the blue nodes (those within the clade). The MST for the red nodes
(those not within the clade and containing the root) intersects with this MST at a single node via the tree popping
theorem. By maintaining the directionality of the parent-child relationships away from the root, we know this is the
node we have already identified. To see why, assume another node in the blue MST is the node of intersection. The parent
of this node must also be blue because otherwise the MST is not minimal. The parent must also lie on the red MST since
by construction following the chain of parents creates a path to the root and this path to the root must be on the red
MST as trees have no cycles. Thus, we have a contradiction of the blue and red MSTs only intersecting at a single node,
so the intersection node must be the previously identified node.

From here we move the labels in the parent associated with the blue MST to the new node. While conceptually simple, the
execution is slightly involved since we first have to copy over the labels and then delete them.

DEPENDENCIES
../../../src/draw.py
../../../src/ortho_tree/tree.py
../phyml_GTR/phyml_GTR.py
    ../phyml_GTR/out/*/meta_*.phy_phyml_tree.txt
"""