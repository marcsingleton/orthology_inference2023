"""Classes and functions shared between fitting and plotting."""

import numpy as np
from numpy import exp, log
from scipy.special import beta, comb


class ArrayRV:
    """Class that allows HMM class to interface with pre-computed probability array.

    The rvs method is only defined so the HMM recognizes it as a proper
    random variable. Since random variates are not needed in this script, the
    body of the method is left blank.
    """
    def __init__(self, array):
        self.array = array

    def pmf(self, x):
        return self.array[x]

    def rvs(self, random_state=None):
        pass


def get_betabinom_pmf(x, n, a, b):
    """Return pmf of beta-binomial distribution evaluated at x."""
    return comb(n, x) * beta(x + a, n - x + b) / beta(a, b)


def get_tree_pmf(tree, pi, q0, q1):
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
