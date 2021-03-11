"""Make consensus trees from bootstrapped PhyML runs."""

import os

import matplotlib.pyplot as plt
import skbio
from src.draw import plot_tree


def majority_consensus(trees, cutoff=0.5):
    # Count nodes
    counts = {}
    for tree in trees:
        for node in tree.traverse():
            if not node.is_tip():
                tip_set = frozenset([tip.name for tip in node.tips()])
                counts[tip_set] = counts.get(tip_set, 0) + 1

    # Make count list and node dictionary
    # This uses the top-level node that includes all taxon labels and is equivalent to initiating the algorithm with a
    # single node containing said labels
    counts = sorted(counts.items(), key=lambda x: (x[1], len(x[0])), reverse=True)
    consensus_nodes = {frozenset([name]): skbio.TreeNode(name=name) for name in counts[0][0]}

    # Make root
    root = skbio.TreeNode(children=list(consensus_nodes.values()))
    consensus_nodes[frozenset([node.name for node in root.children])] = root

    # Add nodes
    for tip_set1, count in counts[1:]:
        if count / len(trees) < cutoff:
            break

        # Find parent
        for tip_set2 in sorted(consensus_nodes, key=lambda x: len(x), reverse=True):
            compatible = (tip_set1 <= tip_set2 or
                          tip_set2 <= tip_set1 or
                          not (tip_set1 & tip_set2))
            if not compatible:
                break
            if tip_set2 >= tip_set1:
                parent_set = tip_set2  # Smallest superset since sorted by size
        if not compatible:
            continue
        parent_node = consensus_nodes[parent_set]

        # Construct current node
        children = []
        for node in parent_node.children:
            for tip in node.tips(include_self=True):
                if tip.name in tip_set1:
                    children.append(node)
                    break
        node = skbio.TreeNode(children=children, support=count/len(trees))  # Instantiation automatically updates parents
        consensus_nodes[tip_set1] = node
        parent_node.append(node)

    # Get branch lengths
    bl_sums = {}
    for tree in trees:
        for node in tree.traverse():
            if node.is_tip():
                continue
            tip_set = frozenset([tip.name for tip in node.tips()])

            # Add clade
            if tip_set in consensus_nodes:
                count, bl_sum = bl_sums.get(tip_set, (0, 0))
                count, bl_sum = count + 1, bl_sum + node.length if node.length else None
                bl_sums[tip_set] = (count, bl_sum)

                # Add terminal children of clade
                for child_node in filter(lambda x: not x.is_tip(), node.children):
                    if not frozenset([tip.name for tip in child_node.tips()]) in consensus_nodes:
                        break
                else:
                    for child_tip in filter(lambda x: x.is_tip(), node.children):
                        tip_set = frozenset([child_tip.name])
                        count, bl_sum = bl_sums.get(tip_set, (0, 0))
                        count, bl_sum = count + 1, bl_sum + child_tip.length if child_tip.length else None
                        bl_sums[tip_set] = (count, bl_sum)

    # Set branch lengths
    for tip_set, (count, bl_sum) in bl_sums.items():
        node = consensus_nodes[tip_set]
        node.length = bl_sum / count if bl_sum else None

    return root


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
../phyml_GTR/phyml_GTR.py
    ../phyml_GTR/out/*/meta_*.phy_phyml_tree.txt
"""