"""Classes and functions for fitting phylogenetic models."""

import numpy as np
from numpy import exp, log
from scipy.special import beta, comb, digamma, gammainc
from scipy.stats import gamma as gamma_rv


class ArrayRV:
    """Class that allows HMM class to interface with a pre-computed probability array.

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


def get_tree_pmf(tree, pinv, k, alpha, pi, q0, q1, p0, p1):
    """Return probability of tree given tips."""
    likelihoods = []
    rates = get_rates(pinv, k, alpha)
    for rate, prior in rates:
        s, conditional = get_conditional(tree, rate * q0, rate * q1, p0, p1)
        likelihood = ([[1-pi], [pi]] * exp(s) * conditional).sum(axis=0)
        likelihoods.append(likelihood * prior)
    return np.stack(likelihoods).sum(axis=0)


def get_tip_pmf(tree, tip, pinv, k, alpha, pi, q0, q1, p0, p1):
    """Return pmf of tip given other tips."""
    pmf0 = get_tree_pmf(tree, pinv, k, alpha, pi, q0, q1, p0, p1)

    tip.value = 1 - tip.value  # Flip state of given tip
    pmf1 = get_tree_pmf(tree, pinv, k, alpha, pi, q0, q1, p0, p1)
    tip.value = 1 - tip.value  # Flip state of given tip back

    return pmf0 / (pmf0 + pmf1)


def get_conditional(tree, q0, q1, p0, p1, inplace=False):
    """Return conditional probabilities of tree given tips and node state."""
    if not inplace:
        tree = tree.copy()  # Make copy so computations do not change original tree

    for node in tree.postorder():
        if node.is_tip():
            node.s = np.zeros(node.value.shape[1])
            node.conditional = np.matmul([[1-p0, p0], [p1, 1-p1]], node.value)
        else:
            ss, ps = [], []
            for child in node.children:
                s, conditional = child.s, child.conditional
                m = get_transition_matrix(q0, q1, child.length)
                p = np.matmul(m, conditional)

                ss.append(s)
                ps.append(p)

            conditional = np.product(np.stack(ps), axis=0)
            s = conditional.sum(axis=0)
            node.conditional = conditional / s  # Normalize to 1 to prevent underflow
            node.s = log(s) + np.sum(np.stack(ss), axis=0)  # Pass forward scaling constant in log space

    return tree.s, tree.conditional


def get_rates(pinv, k, alpha):
    """Return rates and priors for models with site-specific rate variation."""
    igfs = []  # Incomplete gamma function evaluations
    for i in range(k+1):
        x = gamma_rv.ppf(i/k, a=alpha, scale=1/alpha)
        igfs.append(gammainc(alpha+1, alpha*x))
    rates = [(0, pinv)]
    for i in range(k):
        rate = k/(1-pinv) * (igfs[i+1] - igfs[i])
        rates.append((rate, (1-pinv)/k))
    return rates


def get_transition_matrix(q0, q1, t):
    """Return transition matrix for two-state CTMC."""
    q = q0 + q1
    if q == 0 or t == 0:
        p00 = 1
        p11 = 1
    else:
        p00 = (q1 + q0 * exp(-q*t)) / q
        p11 = (q0 + q1 * exp(-q*t)) / q
    p01 = 1 - p00
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


def get_tree_prime(tree, pinv, k, alpha, pi, q0, q1, p0, p1, param):
    """Return derivative of probability of tree relative to a given parameter."""
    derivatives = []
    rates = get_rates(pinv, k, alpha)
    if param == 'pinv':
        s, conditional = get_conditional(tree, 0, 0, p0, p1)
        likelihood = ([[1-pi], [pi]] * exp(s) * conditional).sum(axis=0)
        derivatives.append(likelihood)
        for rate, prior in rates[1:]:  # Skip invariant because calculated above
            s, conditional = get_conditional(tree, rate * q0, rate * q1, p0, p1)
            likelihood = ([[1 - pi], [pi]] * exp(s) * conditional).sum(axis=0)
            term1 = -likelihood / k

            derivative0 = get_conditional_prime_q(tree, rate * q0, rate * q1, p0, p1, 'q0')
            derivative0 = ([[1 - pi], [pi]] * derivative0).sum(axis=0)  # Broadcasting magic
            derivative1 = get_conditional_prime_q(tree, rate * q0, rate * q1, p0, p1, 'q1')
            derivative1 = ([[1 - pi], [pi]] * derivative1).sum(axis=0)  # Broadcasting magic
            term2 = prior * (derivative0 * q0 + derivative1 * q1) * rate / (1 - pinv)

            derivatives.append(term1 + term2)
    elif param == 'alpha':
        rates_prime = get_rates_prime_alpha(pinv, k, alpha)
        for (rate, prior), rate_prime in zip(rates[1:], rates_prime[1:]):  # Skip invariant because alpha is not parameter of invariant category
            derivative0 = get_conditional_prime_q(tree, rate * q0, rate * q1, p0, p1, 'q0')
            derivative0 = ([[1 - pi], [pi]] * derivative0).sum(axis=0)  # Broadcasting magic
            derivative1 = get_conditional_prime_q(tree, rate * q0, rate * q1, p0, p1, 'q1')
            derivative1 = ([[1 - pi], [pi]] * derivative1).sum(axis=0)  # Broadcasting magic
            derivative = derivative0 * q0 + derivative1 * q1

            derivatives.append(rate_prime * derivative * prior)
    elif param == 'pi':
        for rate, prior in rates:
            s, conditional = get_conditional(tree, rate * q0, rate * q1, p0, p1)
            derivative = ([[-1], [1]] * exp(s) * conditional).sum(axis=0)  # Broadcasting magic
            derivatives.append(derivative * prior)
    elif param in ['q0', 'q1']:
        for rate, prior in rates[1:]:  # Skip invariant because qs are not parameters of invariant category
            derivative = get_conditional_prime_q(tree, rate * q0, rate * q1, p0, p1, param)
            derivative = ([[1 - pi], [pi]] * derivative).sum(axis=0)  # Broadcasting magic
            derivatives.append(rate * derivative * prior)
    elif param in ['p0', 'p1']:
        for rate, prior in rates:
            derivative = get_conditional_prime_p(tree, rate * q0, rate * q1, p0, p1, param)
            derivative = ([[1 - pi], [pi]] * derivative).sum(axis=0)  # Broadcasting magic
            derivatives.append(derivative * prior)
    else:
        raise ValueError('"param" is not "pinv", "alpha", "pi", "q0", "q1", "p0", or "p1"')

    return np.stack(derivatives).sum(axis=0)


def get_tip_prime(tree, tip, pinv, k, alpha, pi, q0, q1, p0, p1, param):
    """Return derivative of probability of tip relative to a given parameter."""
    pmf0 = get_tree_pmf(tree, pinv, k, alpha, pi, q0, q1, p0, p1)
    pmf0_prime = get_tree_prime(tree, pinv, k, alpha, pi, q0, q1, p0, p1, param)

    tip.value = 1 - tip.value  # Flip state of given tip
    pmf1 = get_tree_pmf(tree, pinv, k, alpha, pi, q0, q1, p0, p1)
    pmf1_prime = get_tree_prime(tree, pinv, k, alpha, pi, q0, q1, p0, p1, param)
    tip.value = 1 - tip.value  # Flip state of given tip back

    d = pmf0 + pmf1
    d_prime = pmf0_prime + pmf1_prime
    return (pmf0_prime * d - pmf0 * d_prime) / d ** 2  # Quotient rule applied to get_tip_pmf


def get_conditional_prime_q(tree, q0, q1, p0, p1, param, inplace=False):
    """Return derivative of conditional probabilities relative to q."""
    if not inplace:
        tree = tree.copy()  # Make copy so computations do not change original tree

    get_conditional(tree, q0, q1, p0, p1, inplace=True)  # Calculate conditionals for use during traversal
    for node in tree.postorder():
        if node.is_tip():
            node.conditional = np.matmul([[1-p0, p0], [p1, 1-p1]], node.value)
            node.derivative = np.zeros((2, node.value.shape[1]))
        else:
            # Collect product and derivative of product for each branch
            ps = []
            dps = []
            for child in node.children:
                s, conditional, derivative = child.s, child.conditional, child.derivative
                conditional = exp(s) * conditional  # Un-normalize

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
            node.derivative = np.sum(np.stack(terms), axis=0)

    return tree.derivative


def get_conditional_prime_p(tree, q0, q1, p0, p1, param, inplace=False):
    """Return derivative of conditional probabilities relative to p."""
    if param == 'p0':
        dm = [[-1, 1], [0, 0]]
    elif param == 'p1':
        dm = [[0, 0], [1, -1]]
    else:
        raise ValueError('"param" is not "p0" or "p1"')
    if not inplace:
        tree = tree.copy()  # Make copy so computations do not change original tree

    get_conditional(tree, q0, q1, p0, p1, inplace=True)  # Calculate conditionals for use during traversal
    for node in tree.postorder():
        if node.is_tip():
            node.conditional = np.matmul([[1 - p0, p0], [p1, 1 - p1]], node.value)
            node.derivative = np.matmul(dm, node.value)  # Use original value
        else:
            # Collect product and derivative of product for each branch
            ps = []
            dps = []
            for child in node.children:
                s, conditional, derivative = child.s, child.conditional, child.derivative
                conditional = exp(s) * conditional  # Un-normalize

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
            node.derivative = np.sum(np.stack(terms), axis=0)

    return tree.derivative


def get_transition_matrix_prime_q(q0, q1, t, param):
    """Return derivative transition matrix for two-state CTMC relative to q."""
    q = q0 + q1
    if param == 'q0':
        d00 = -(q1 + (t * q0 ** 2 + t * q0 * q1 - q1) * exp(-q * t)) / q ** 2
        d11 = q1*(1 - (q0 * t + q1 * t + 1) * exp(-q * t)) / q ** 2
        d01 = -d00
        d10 = -d11
        return np.array([[d00, d01], [d10, d11]])
    elif param == 'q1':
        d00 = q0 * (1 - (t * q0 + t * q1 + 1) * exp(-q * t)) / q ** 2
        d11 = -(q0 + (t * q1 ** 2 + t * q0 * q1 - q0) * exp(-q * t)) / q ** 2
        d01 = -d00
        d10 = -d11
        return np.array([[d00, d01], [d10, d11]])
    else:
        raise ValueError('"param" is not "q0" or "q1"')


def get_rates_prime_alpha(pinv, k, alpha, eps=1E-8):
    """Return derivatives of rates relative to alpha for models with site-specific rate variation.

    Rate derivatives are calculated using the symmetric finite difference
    formula for the gamma function evaluations because the expression for the
    exact derivative was extremely complex and involved a non-standard integral
    derived from the gamma function.
    """
    igfs_prime = []  # Incomplete gamma function evaluations
    for i in range(k+1):
        alpha0 = alpha - eps/2
        alpha1 = alpha + eps/2
        x0 = gamma_rv.ppf(i/k, a=alpha0, scale=1/alpha0)
        x1 = gamma_rv.ppf(i/k, a=alpha1, scale=1/alpha1)
        igfs_prime.append((gammainc(alpha1+1, alpha1*x1) - gammainc(alpha0+1, alpha0*x0)) / eps)
    rates_prime = [0]
    for i in range(k):
        rate_prime = k/(1-pinv) * (igfs_prime[i+1] - igfs_prime[i])
        rates_prime.append(rate_prime)
    return rates_prime
