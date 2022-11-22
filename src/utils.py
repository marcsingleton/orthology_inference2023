"""Functions for common operations in this project."""

from itertools import product

import numpy as np


def get_brownian_weights(tree):
    """Get weights of tip characters using Brownian motion process.

    This model assumes the tip characters are continuous traits that evolve
    according to a Brownian motion process, i.e. they follow a multivariate
    normal distribution where the covariance between two tips is proportional
    to the length of their shared path from the root. The MLE for the root
    value is then a weighted average of the observed tip values with the
    returned weights.

    The formula is given in the appendix of:
    Weights for Data Related by a Tree. SF Altschul et al.
    J. Mol. Biol. (1989) 207, 647-653. 10.1016/0022-2836(89)90234-9

    Parameters
    ----------
    tree: TreeNode (skbio)

    Returns
    -------
    weights: list of tuples of (tipm, weight)
    """
    tree = tree.copy()  # Make copy so computations do not change original tree
    tips = list(tree.tips())
    tip2idx = {tip: i for i, tip in enumerate(tips)}
    idx2tip = {i: tip for i, tip in enumerate(tips)}

    # Accumulate tip names up to root
    for node in tree.postorder():
        if node.is_tip():
            node.tip_nodes = {node}
        else:
            node.tip_nodes = set.union(*[child.tip_nodes for child in node.children])

    # Fill in covariance matrix from root to tips
    tree.root_length = 0
    cov = np.zeros((len(tips), len(tips)))
    for node in tree.traverse(include_self=True):
        for child in node.children:
            child.root_length = node.root_length + child.length
        if not node.is_tip():
            child1, child2 = node.children
            idxs1, idxs2 = [tip2idx[tip] for tip in child1.tip_nodes], [tip2idx[tip] for tip in child2.tip_nodes]
            for idx1, idx2 in product(idxs1, idxs2):
                cov[idx1, idx2] = node.root_length
                cov[idx2, idx1] = node.root_length
        else:
            idx = tip2idx[node]
            cov[idx, idx] = node.root_length

    # Compute weights
    # The formula below is from the appendix of the referenced work
    inv = np.linalg.inv(cov)
    row_sum = inv.sum(axis=1)
    total_sum = inv.sum()
    weights = row_sum / total_sum
    return [(idx2tip[idx], weight) for idx, weight in enumerate(weights)]


def read_fasta(path):
    """Read FASTA file at path and return list of headers and sequences.

    Parameters
    ----------
    path: str
        Path to FASTA file

    Returns
    -------
    fasta: list of tuples of (header, seq)
        The first element is the header line with the >, and the second
        element is the corresponding sequence.
    """
    fasta = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line.rstrip()
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            fasta.append((header, seq))
    return fasta
