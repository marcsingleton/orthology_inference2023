"""Algorithms for trees not found in skbio."""

import skbio


def majority_consensus(trees, cutoff=0.5):
    """Compute the majority consensus tree.

    When a taxon is designated as an outgroup, this algorithm is equivalent to
    tree popping described in The Mathematics of Phylogenetics. Since a clade
    is the subset of a split that does not contain the root, all sets of
    terminal labels obtained via the tips method when traversing the tree
    correspond to the half of the split without the root. Thus operating on the
    clades is implicitly operating on the equivalent split.

    The overall approach of the algorithm is to first tally the clades and
    then to progressively build the consensus tree by applying the clades in
    decreasing order of frequency (as long as they are compatible with all
    clades already displayed on the tree). The proof presented in The
    Mathematics of Phylogenetics is considerably more complex than the more
    intuitive algorithm used here. Whereas that proof required finding the
    minimum spanning tree, this approach instead uses the idea that as long as
    clades are applied in order of size (which is guaranteed by the principle
    that a parent clade will appear at least as many times in counts as its
    subclades and subsequent sorting on counts then size), then the right node
    onto which to apply a clade is always the smallest one that is compatible.

    This actually equivalent to finding the minimum spanning tree of the nodes
    colored by each clade (blue for those in the clade and red for those not in
    the clade). The tree induced by a node whose terminals are a superset of
    the clade clearly spans that clade. Thus, by iterating through all the
    nodes in the consensus tree sorted by decreasing size, the last node that
    meets this criterion induces the minimal spanning tree (MST) for the blue
    nodes (those within the clade). The MST for the red nodes (those not within
    the clade and containing the root) intersects with this MST at a single
    node via the tree popping theorem. By maintaining the directionality of the
    parent-child relationships away from the root, we know this is the node we
    have already identified. To see why, assume another node in the blue MST is
    the node of intersection. The parent of this node must also be blue because
    otherwise the MST is not minimal.The parent must also lie on the red MST
    since by construction following the chain of parents creates a path to the
    root and this path to the root must be on the red MST as trees have no
    cycles. Thus, we have a contradiction of the blue and red MSTs only
    intersecting at a single node, so the intersection node must be the
    previously identified node.

    Parameters
    ----------
        trees: list of skbio TreeNodes
            Input trees from which to compute the consensus.
        cutoff: float
            The minimum fraction of trees a clade must be displayed on in order
            (not including the value itself) to contribute to the consensus.
            The default value of 0.5 ensures all incorporated clades are
            compatible. If the cutoff is set lower, a clade is only added to
            the consensus if it compatible with all clades already displayed.
            Clades are added in order of decreasing number of counts in the
            input trees.

    Returns
    -------
        root: TreeNode
    """
    # Count nodes and record tip names
    counts, tip_names = {}, set()
    for tree in trees:
        tip_names.update([tip.name for tip in tree.tips()])
        for node in tree.traverse():
            if not node.is_tip():
                tip_set = frozenset([tip.name for tip in node.tips()])
                counts[tip_set] = counts.get(tip_set, 0) + 1
    counts = sorted(counts.items(), key=lambda x: (x[1], len(x[0])), reverse=True)  # Sort by count then size

    # Make count list and node dictionary
    # Initialize consensus nodes with a tip for each unique taxon label and a root node containing all these tips
    # (The tips are instantiated first, and then the root node containing them is created.)
    consensus_nodes = {frozenset([name]): skbio.TreeNode(name=name) for name in tip_names}
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
                break  # Break search over parents
            if tip_set2 >= tip_set1:
                parent_set = tip_set2  # Smallest superset since sorted by size
        if not compatible:
            continue  # Each clade should be compatible with all others in consensus tree
        parent_node = consensus_nodes[parent_set]

        # Construct current node
        children = []
        for node in parent_node.children:
            for tip in node.tips(include_self=True):
                if tip.name in tip_set1:
                    children.append(node)
                    break  # Compatibility requires all tips of a subnode are included if one is
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
                # The lengths of terminal children are only included if all the non-terminal children are clades
                # displayed in the consensus tree. This is expressed by checking this condition in a for loop and
                # proceeding to the else clause only if the loop isn't broken
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
