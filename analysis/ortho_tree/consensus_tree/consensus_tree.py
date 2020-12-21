"""Make consensus trees from bootstrapped PhyML runs."""

import os

import Bio.Phylo as Phylo
import matplotlib.pyplot as plt


def majority_consensus(trees, cutoff=0.5):
    # Count clades
    counts = {}
    for tree in trees:
        for clade in tree.find_clades(terminal=False):
            terminal_set = frozenset([terminal.name for terminal in clade.get_terminals()])
            counts[terminal_set] = counts.get(terminal_set, 0) + 1

    # Make count list and clade dictionary
    # This uses the top-level clade that includes all taxon labels and is equivalent to initiating the algorithm with a
    # single node containing said labels
    counts = sorted(counts.items(), key=lambda x: (x[1], len(x[0])), reverse=True)
    consensus_clades = {frozenset([name]): Phylo.BaseTree.Clade(name=name) for name in counts[0][0]}

    # Make root
    root = Phylo.BaseTree.Clade()
    root.clades = list(consensus_clades.values())
    consensus_clades[frozenset([clade.name for clade in root.clades])] = root

    # Add clades
    for terminal_set1, count in counts[1:]:
        if count / len(trees) < cutoff:
            break

        # Find parent
        for terminal_set2 in sorted(consensus_clades, key=lambda x: len(x), reverse=True):
            compatible = (terminal_set1 <= terminal_set2 or
                          terminal_set2 <= terminal_set1 or
                          not (terminal_set1 & terminal_set2))
            if not compatible:
                break
            if terminal_set2 >= terminal_set1:
                parent_set = terminal_set2  # Smallest superset since sorted by size
        if not compatible:
            continue
        parent_clade = consensus_clades[parent_set]

        # Construct current clade
        child_clades = []
        for clade in parent_clade.clades:
            for terminal in clade.get_terminals():
                if terminal.name in terminal_set1:
                    child_clades.append(clade)
                    break
        clade = Phylo.BaseTree.Clade(clades=child_clades, confidence=count/len(trees))
        consensus_clades[terminal_set1] = clade

        # Remove clades from parent
        parent_clade.clades.append(clade)  # Add current clade
        parent_clade.clades = [clade for clade in parent_clade.clades if clade not in child_clades]

    # Get branch lengths
    bl_sums = {}
    for tree in trees:
        for clade in tree.find_clades(terminal=False):
            terminal_set = frozenset([terminal.name for terminal in clade.get_terminals()])

            # Add clade
            if terminal_set in consensus_clades:
                count, bl_sum = bl_sums.get(terminal_set, (0, 0))
                count, bl_sum = count + 1, bl_sum + clade.branch_length if clade.branch_length else None
                bl_sums[terminal_set] = (count, bl_sum)

                # Add terminal children of clade
                for child_clade in filter(lambda x: not x.is_terminal(), clade.clades):
                    if not frozenset([terminal.name for terminal in child_clade.get_terminals()]) in consensus_clades:
                        break
                else:
                    for child_terminal in filter(lambda x: x.is_terminal(), clade.clades):
                        terminal_set = frozenset([child_terminal.name])
                        count, bl_sum = bl_sums.get(terminal_set, (0, 0))
                        count, bl_sum = count + 1, bl_sum + child_terminal.branch_length if child_terminal.branch_length else None
                        bl_sums[terminal_set] = (count, bl_sum)

    # Set branch lengths
    for terminal_set, (count, bl_sum) in bl_sums.items():
        clade = consensus_clades[terminal_set]
        clade.branch_length = bl_sum / count if bl_sum else None

    return Phylo.BaseTree.Tree(root=root)


if not os.path.exists('out/'):
    os.mkdir('out/')

for label in os.listdir('../phyml_GTR/out/'):
    trees = []
    for file in filter(lambda x: x.endswith('.phy_phyml_tree.txt'), os.listdir(f'../phyml_GTR/out/{label}/')):
        tree = Phylo.read(f'../phyml_GTR/out/{label}/{file}', 'newick')
        for clade in tree.find_clades():
            clade.confidence = None
        path = tree.get_path('dwil')
        if len(path) == 1:
            outgroup = tree.clade
        else:
            outgroup = path[-2]
        tree.root_with_outgroup(outgroup)
        trees.append(tree)
    ctree = majority_consensus(trees)
    ctree.ladderize()
    Phylo.write(ctree, f'out/{label}.txt', 'newick')

    # Save image as PNG
    Phylo.draw(ctree, show_confidence=False, do_show=False)
    plt.gca().yaxis.set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.xlabel('')
    plt.savefig(f'out/{label}.png')
    plt.close()

    # Save image as SVG (with confidences)
    Phylo.draw(ctree, do_show=False)
    plt.gca().yaxis.set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.xlabel('')
    plt.savefig(f'out/{label}.svg')
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
../phyml_GTR/phyml_GTR.py
    ../phyml_GTR/out/*/meta_*.phy_phyml_tree.txt
"""