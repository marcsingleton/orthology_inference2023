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


def get_tree_pmf(tree, pi, q0, q1, r):
    """Return probability of tree given tips."""
    s, conditional = get_conditional(tree, q0, q1, r)
    l = ((exp(s) * conditional) * [[1-pi], [pi]]).sum(axis=0)
    return l


def get_conditional(node, q0, q1, r):
    """Return conditional probabilities of tree given tips and node state."""
    ss, ps = [], []
    for child in node.children:
        if child.is_tip():
            s, conditional = np.zeros(child.conditional.shape[1]), child.conditional
            r_child = r
        else:
            s, conditional = get_conditional(child, q0, q1, r)
            r_child = 0
        m = get_transition_matrix(q0, q1, r_child, child.length)
        p = np.matmul(m, conditional)

        ss.append(s)
        ps.append(p)

    conditional = np.product(np.stack(ps), axis=0)
    s = conditional.sum(axis=0)
    conditional = conditional / s  # Normalize to 1 to prevent underflow
    s = log(s) + np.sum(np.stack(ss), axis=0)  # Pass forward scaling constant in log space

    return s, conditional


def get_transition_matrix(q0, q1, r, t):
    """Return transition matrix for two-state CTMC."""
    q = q0 + q1
    t += r
    p00 = (q1 + q0 * exp(-q*t)) / q
    p01 = 1 - p00
    p11 = (q0 + q1 * exp(-q*t)) / q
    p10 = 1 - p11
    return np.array([[p00, p01], [p10, p11]])
