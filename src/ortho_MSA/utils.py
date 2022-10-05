"""Classes and functions shared between fitting and plotting."""

import numpy as np
from numpy import exp, log
from scipy.special import beta, comb, digamma


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


# Probability functions
def get_bernoulli_pmf(x, p):
    """Return pmf of Bernoulli distribution evaluated at x."""
    return p**x * (1 - p)**(1 - x)


def get_betabinom_pmf(x, n, a, b):
    """Return pmf of beta-binomial distribution evaluated at x."""
    return comb(n, x) * beta(x + a, n - x + b) / beta(a, b)


def get_tree_pmf(tree, pi, q0, q1, p0, p1):
    """Return probability of tree given tips."""
    s, conditional = get_conditional(tree, q0, q1, p0, p1)
    l = ((exp(s) * conditional) * [[1-pi], [pi]]).sum(axis=0)
    return l


def get_tip_pmf(tree, tip_name, pi, q0, q1, p0, p1):
    """Return pmf of tip with tip_name given other tips."""
    pmf0 = get_tree_pmf(tree, pi, q0, q1, p0, p1)

    tip = tree.tip_dict[tip_name]
    tip.conditional = 1 - tip.conditional  # Flip state of given tip
    pmf1 = get_tree_pmf(tree, pi, q0, q1, p0, p1)
    tip.conditional = 1 - tip.conditional  # Flip state of given tip back

    return pmf0 / (pmf0 + pmf1)


def get_conditional(node, q0, q1, p0, p1):
    """Return conditional probabilities of tree given tips and node state."""
    ss, ps = [], []
    for child in node.children:
        if child.is_tip():
            s, conditional = np.zeros(child.conditional.shape[1]), child.conditional
            conditional = np.matmul([[1-p0, p0], [p1, 1-p1]], conditional)
        else:
            s, conditional = get_conditional(child, q0, q1, p0, p1)
        m = get_transition_matrix(q0, q1, child.length)
        p = np.matmul(m, conditional)

        ss.append(s)
        ps.append(p)

    conditional = np.product(np.stack(ps), axis=0)
    s = conditional.sum(axis=0)
    conditional = conditional / s  # Normalize to 1 to prevent underflow
    s = log(s) + np.sum(np.stack(ss), axis=0)  # Pass forward scaling constant in log space

    return s, conditional


def get_transition_matrix(q0, q1, t):
    """Return transition matrix for two-state CTMC."""
    q = q0 + q1
    p00 = (q1 + q0 * exp(-q*t)) / q
    p01 = 1 - p00
    p11 = (q0 + q1 * exp(-q*t)) / q
    p10 = 1 - p11
    return np.array([[p00, p01], [p10, p11]])


# Gradient functions
def get_bernoulli_prime(x, p):
    """Return derivative of Bernoulli pmf relative to p."""
    return get_bernoulli_pmf(x, p) * (x/p - (1-x)/(1-p))


def get_betabinom_prime(x, n, a, b, param):
    """Return derivative of beta-binomial pmf relative to a given parameter."""
    if param == 'a':
        return comb(n, x) * (beta_prime(x + a, n - x + b) * beta(a, b) - beta(x + a, n - x + b) * beta_prime(a, b)) / (beta(a, b)) ** 2
    elif param == 'b':
        return comb(n, x) * (beta_prime(n - x + b, x + a) * beta(a, b) - beta(x + a, n - x + b) * beta_prime(b, a)) / (beta(a, b)) ** 2
    else:
        raise ValueError('"param" is not "a" or "b"')


def beta_prime(a, b):
    """Return derivative of beta function relative to its first parameter, a."""
    return beta(a, b) * (digamma(a) - digamma(a + b))


def get_tree_prime(tree, pi, q0, q1, p0, p1, param):
    """Return derivative of probability of tree relative to a given parameter."""
    if param == 'pi':
        s, conditional = get_conditional(tree, q0, q1, p0, p1)
        d = ((exp(s) * conditional) * [[-1], [1]]).sum(axis=0)  # Broadcasting magic
        return d
    elif param in ['q0', 'q1']:
        derivative = get_conditional_prime_q(tree, q0, q1, p0, p1, param)
        d = (derivative * [[1 - pi], [pi]]).sum(axis=0)  # Broadcasting magic
        return d
    elif param in ['p0', 'p1']:
        derivative = get_conditional_prime_p(tree, q0, q1, p0, p1, param)
        d = (derivative * [[1 - pi], [pi]]).sum(axis=0)  # Broadcasting magic
        return d
    else:
        raise ValueError('"param" is not "pi", "q0", "q1", "p0", or "p1"')


def get_tip_prime(tree, tip_name, pi, q0, q1, p0, p1, param):
    """Return derivative of probability of tip relative to a given parameter."""
    pmf0 = get_tree_pmf(tree, pi, q0, q1, p0, p1)
    pmf0_prime = get_tree_prime(tree, pi, q0, q1, p0, p1, param)

    tip = tree.tip_dict[tip_name]
    tip.conditional = 1 - tip.conditional  # Flip state of given tip
    pmf1 = get_tree_pmf(tree, pi, q0, q1, p0, p1)
    pmf1_prime = get_tree_prime(tree, pi, q0, q1, p0, p1, param)
    tip.conditional = 1 - tip.conditional  # Flip state of given tip back

    d = pmf0 + pmf1
    d_prime = pmf0_prime + pmf1_prime
    return (pmf0_prime * d - pmf0 * d_prime) / d ** 2  # Quotient rule applied to get_tip_pmf


def get_conditional_prime_q(node, q0, q1, p0, p1, param):
    """Return derivative of conditional probabilities relative to q."""
    # Collect product and derivative of product for each branch
    ps = []
    dps = []
    for child in node.children:
        if child.is_tip():
            conditional = child.conditional
            conditional = np.matmul([[1-p0, p0], [p1, 1-p1]], conditional)
            derivative = np.zeros((2, conditional.shape[1]))
        else:
            s, conditional = get_conditional(child, q0, q1, p0, p1)
            conditional = exp(s) * conditional  # Un-normalize
            derivative = get_conditional_prime_q(child, q0, q1, p0, p1, param)

        m = get_transition_matrix(q0, q1, child.length)
        p = np.matmul(m, conditional)

        dm = get_transition_matrix_prime_q(q0, q1, child.length, param)
        dp = np.matmul(dm, conditional) + np.matmul(m, derivative)

        ps.append(p)
        dps.append(dp)

    # Assemble products and derivatives into terms and then sum (product rule)
    terms = []
    for i in range(len(node.children)):
        term = [dps[j] if i == j else ps[j] for j in range(len(node.children))]
        term = np.product(np.stack(term), axis=0)
        terms.append(term)
    derivative = np.sum(np.stack(terms), axis=0)

    return derivative


def get_conditional_prime_p(node, q0, q1, p0, p1, param):
    """Return derivative of conditional probabilities relative to p."""
    # Collect product and derivative of product for each branch
    ps = []
    dps = []
    for child in node.children:
        if child.is_tip():
            conditional = child.conditional
            conditional = np.matmul([[1 - p0, p0], [p1, 1 - p1]], conditional)
            if param == 'p0':
                dm = [[-1, 1], [0, 0]]
            elif param == 'p1':
                dm = [[0, 0], [1, -1]]
            else:
                raise ValueError('"param" is not "p0" or "p1"')
            derivative = np.matmul(dm, child.conditional)  # Use original conditional
        else:
            s, conditional = get_conditional(child, q0, q1, p0, p1)
            conditional = exp(s) * conditional  # Un-normalize
            derivative = get_conditional_prime_p(child, q0, q1, p0, p1, param)

        m = get_transition_matrix(q0, q1, child.length)
        p = np.matmul(m, conditional)

        dp = np.matmul(m, derivative)

        ps.append(p)
        dps.append(dp)

    # Assemble products and derivatives into terms and then sum (product rule)
    terms = []
    for i in range(len(node.children)):
        term = [dps[j] if i == j else ps[j] for j in range(len(node.children))]
        term = np.product(np.stack(term), axis=0)
        terms.append(term)
    derivative = np.sum(np.stack(terms), axis=0)

    return derivative


def get_transition_matrix_prime_q(q0, q1, t, param):
    """Return derivative transition matrix for two-state CTMC relative to q."""
    q = q0 + q1
    if param == 'q0':
        d00 = -(q1 + (t * q0 ** 2 + t * q0 * q1 - q1) * exp(-q * t)) / q ** 2
        d01 = -d00
        d11 = q1*(1 - (q0 * t + q1 * t + 1) * exp(-q * t)) / q ** 2
        d10 = -d11
        return np.array([[d00, d01], [d10, d11]])
    elif param == 'q1':
        d00 = q0 * (1 - (t * q0 + t * q1 + 1) * exp(-q * t)) / q ** 2
        d01 = -d00
        d11 = -(q0 + (t * q1 ** 2 + t * q0 * q1 - q0) * exp(-q * t)) / q ** 2
        d10 = -d11
        return np.array([[d00, d01], [d10, d11]])
    else:
        raise ValueError('"param" is not "q0" or "q1"')
