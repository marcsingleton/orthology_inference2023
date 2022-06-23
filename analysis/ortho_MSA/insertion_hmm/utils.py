"""Classes and functions shared between fitting and plotting."""

import numpy as np
from numpy import exp, log


class BinomialArrayRV:
    """Class that allows HMM class to interface with pre-computed probability array.

    The rvs method is only defined so the HMM recognizes it as a proper
    random variable. Since random variates are not needed in this script, the
    body of the method is left blank.
    """
    def __init__(self, p, array):
        self.p = p
        self.array = array

    def pmf(self, x):
        p1 = self.p if x[0] else 1 - self.p
        p2 = self.array[x[1]]
        return p1 * p2

    def rvs(self, size=None):
        pass


def get_tip_posterior(tree, ppid, pi, q0, q1):
    """Return probability of tip given tree."""
    tree1 = tree
    tree2 = tree.deepcopy()
    for tip in tree2.tips():
        if tip.ppid == ppid:
            tip.conditional = 1 - tip.conditional  # Flip state of given tip

    p1 = get_tree_probability(tree1, pi, q0, q1)
    p2 = get_tree_probability(tree2, pi, q0, q1)
    return p1 / (p1 + p2)


def get_tree_probability(tree, pi, q0, q1):
    """Return probability of tree given tips."""
    s, conditional = get_conditional(tree, q0, q1)
    l = ((exp(s) * conditional) * [[1-pi], [pi]]).sum(axis=0)
    return l


def get_conditional(node, q0, q1):
    """Return conditional probabilities of tree given tips and node state."""
    child1, child2 = node.children

    if child1.is_tip():
        s1, conditional1 = 0, child1.conditional
    else:
        s1, conditional1 = get_conditional(child1, q0, q1)
    if child2.is_tip():
        s2, conditional2 = 0, child2.conditional
    else:
        s2, conditional2 = get_conditional(child2, q0, q1)

    p1 = get_transition_matrix(q0, q1, child1.length)
    p2 = get_transition_matrix(q0, q1, child2.length)
    conditional = np.matmul(p1, conditional1) * np.matmul(p2, conditional2)
    s = conditional.sum(axis=0)
    conditional = conditional / s  # Normalize to 1 to prevent underflow
    s = log(s) + s1 + s2  # Pass forward scaling constant in log space

    return s, conditional


def get_transition_matrix(q0, q1, t):
    """Return transition matrix for two-state CTMC."""
    q = q0 + q1
    p00 = (q1 + q0 * exp(-q*t)) / q
    p01 = 1 - p00
    p11 = (q0 + q1 * exp(-q*t)) / q
    p10 = 1 - p11
    return np.array([[p00, p01], [p10, p11]])
