"""Common functions for ASR."""

import numpy as np
import scipy.linalg as linalg
import skbio


def get_conditional(node, matrix):
    """Return conditional probabilities of tree given tips and node state."""
    child1, child2 = node.children

    if child1.is_tip():
        s1, conditional1 = 0, child1.conditional
    else:
        s1, conditional1 = get_conditional(child1, matrix)
    if child2.is_tip():
        s2, conditional2 = 0, child2.conditional
    else:
        s2, conditional2 = get_conditional(child2, matrix)

    p1 = linalg.expm(matrix * child1.length)
    p2 = linalg.expm(matrix * child2.length)
    conditional = np.matmul(p1, conditional1) * np.matmul(p2, conditional2)
    s = conditional.sum(axis=0)
    conditional = conditional / s  # Normalize to 1 to prevent underflow
    s = np.log(s) + s1 + s2  # Pass forward scaling constant in log space

    return s, conditional


def get_tree(tree1, tree2):
    """Root fitted tree such that root creates splits compatible with original tree.

    Since an arbitrarily rooted tree can split only one of the two true root clades, we check for both when traversing
    the tree. The root is temporarily set as the parent node of the first of the root clades that is found. As the
    re-rooted tree has three children at the root, this has the effect of splitting that root clade while preserving the
    other root clade. For simplicity, the unsplit root clade is designated as nodeA and the other as nodeB. Once these
    labels are defined, the new root is introduced between these nodes with a small amount of "tree surgery."
    """
    # Find split and calculate ratio
    tree1 = tree1.shear([tip.name for tip in tree2.tips()])
    nodeA, nodeB = tree1.children
    tipsA, tipsB = set([tip.name for tip in nodeA.tips()]), set([tip.name for tip in nodeB.tips()])
    for node in tree2.traverse():
        tips = set([tip.name for tip in node.tips()])
        if tips == tipsA:
            nodeA, nodeB = nodeB, nodeA  # Swap labels so nodeA is the preserved root clade
            ratio = nodeA.length / (nodeA.length + nodeB.length)
            break
        elif tips == tipsB:
            ratio = nodeA.length / (nodeA.length + nodeB.length)
            break
    tree2 = tree2.root_at(node)

    # Insert root node
    nodeB1, nodeB2 = None, None
    for child in tree2.children:
        if child.count(tips=True) == nodeA.count(tips=True):
            nodeA = child
        elif nodeB1 is None:
            nodeB1 = child
        else:
            nodeB2 = child
    lengthA = ratio * nodeA.length
    lengthB = (1 - ratio) * nodeA.length

    nodeA = skbio.TreeNode(children=nodeA.children, length=lengthA)
    nodeB = skbio.TreeNode(children=[nodeB1, nodeB2], length=lengthB)
    tree = skbio.TreeNode(children=[nodeA, nodeB])

    return tree
