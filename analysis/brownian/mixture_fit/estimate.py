"""Functions for estimators of distributions with weighted data."""

import numpy as np
import scipy.optimize as opt
from math import exp, log, pi, sqrt
from scipy.special import digamma


# Functions for MLEs with no closed-form solutions
def create_fisk_scale(data, expt=None):
    expt = np.full(len(data), 1) if expt is None else expt

    def fisk_scale(scale):
        # Compute sums
        e = expt.sum()
        q = ((expt * data) / (scale + data)).sum()

        return 2 * q - e

    return fisk_scale


def create_fisk_shape(data, expt=None, scale=1):
    expt = np.full(len(data), 1) if expt is None else expt

    def fisk_shape(c):
        # Compute summands
        r = data / scale
        s = 1 / c + np.log(r) - 2 * np.log(r) * r ** c / (1 + r ** c)

        return (expt * s).sum()

    return fisk_shape


def create_gamma_shape(data, expt=None):
    expt = np.full(len(data), 1) if expt is None else expt

    def gamma_shape(a):
        # Compute sums
        e = expt.sum()
        ed = (expt * data).sum()
        elogd = (expt * np.log(data)).sum()

        return elogd - e * log(ed / e) + e * (log(a) - digamma(a))

    return gamma_shape


# Maximum likelihood estimators
def mle_expon(data, param_fix=None, expt=None, **kwargs):
    param_fix, expt = param_check(data, param_fix, expt)

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        e = expt.sum()
        ed = (expt * data).sum()
        scale = ed / e

    return {'scale': scale}


def mle_fisk(data, param_fix=None, param=None, expt=None):
    param_fix, expt = param_check(data, param_fix, expt)
    param = mm_fisk(data) if param is None else param

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        fisk_scale = create_fisk_scale(data, expt)
        scale = opt.newton(fisk_scale, param['scale'])

    # Shape
    if 'c' in param_fix:
        c = param_fix['c']
    else:
        fisk_shape = create_fisk_shape(data, expt, scale)
        c = opt.newton(fisk_shape, param['c'])

    return {'c': c, 'scale': scale}


def mle_gamma(data, param_fix=None, param=None, expt=None):
    param_fix, expt = param_check(data, param_fix, expt)
    param = mm_gamma(data) if param is None else param

    # Shape
    if 'a' in param_fix:
        a = param_fix['a']
    else:
        gamma_shape = create_gamma_shape(data, expt)
        try:
            a = opt.newton(gamma_shape, param['a'])
        except ValueError:  # Catch an error raised by a trial value below zero
            a = opt.brentq(gamma_shape, 1E-10, 2 * param['a'])  # Attempt Brent's method on (1E-9, 2a)

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        scale = (expt * data).sum() / (a * expt.sum())

    return {'a': a, 'scale': scale}


def mle_levy(data, param_fix=None, expt=None, **kwargs):
    param_fix, expt = param_check(data, param_fix, expt)

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        e = expt.sum()
        edivd = (expt / data).sum()
        scale = e / edivd

    return {'scale': scale}


def mle_lognorm(data, param_fix=None, expt=None, **kwargs):
    param_fix, expt = param_check(data, param_fix, expt)

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        e = expt.sum()
        elogd = (expt * np.log(data)).sum()
        scale = exp(elogd / e)

    # Shape
    if 's' in param_fix:
        s = param_fix['s']
    else:
        logd_sqerr = (np.log(data) - log(scale)) ** 2
        s = sqrt((expt * logd_sqerr).sum() / e)

    return {'s': s, 'scale': scale}


def mle_pareto(data, param_fix, expt=None, **kwargs):
    param_fix, expt = param_check(data, param_fix, expt)

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        scale = min(data)

    # Shape
    if 'b' in param_fix:
        b = param_fix['b']
    else:
        e = expt.sum()
        elogd = (expt * np.log(data)).sum()
        b = e / (elogd - e * log(scale))

    return {'b': b, 'scale': scale}


def mle_uniform(data, param_fix=None, **kwargs):
    param_fix = {} if param_fix is None else param_fix

    # Loc
    if 'loc' in param_fix:
        loc = param_fix['loc']
    else:
        loc = min(data)

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        scale = max(data) - loc

    return {'loc': loc, 'scale': scale}


# Method of moments estimators (for providing initial values for MLEs without closed forms)
def mm_fisk(data, expt=None, **kwargs):
    expt = np.full(len(data), 1) if expt is None else expt

    # Moments
    logdata = np.log(data)
    m1 = (logdata * expt).sum() / expt.sum()
    m2 = (logdata ** 2 * expt).sum() / expt.sum()

    return {'c': pi / sqrt(3 * (m2 - m1 ** 2)), 'scale': exp(m1)}


def mm_gamma(data, expt=None, **kwargs):
    expt = np.full(len(data), 1) if expt is None else expt

    # Moments
    m1 = (data * expt).sum() / expt.sum()
    m2 = (data ** 2 * expt).sum() / expt.sum()

    return {'a': m1 ** 2 / (m2 - m1 ** 2), 'scale': (m2 - m1 ** 2) / m1}


def mm_lognorm(data, expt=None, **kwargs):
    expt = np.full(len(data), 1) if expt is None else expt

    # Moments
    m1 = (data * expt).sum() / expt.sum()
    m2 = (data ** 2 * expt).sum() / expt.sum()

    # Parameters of transformed lognorm
    mu = 2 * log(m1) - 0.5 * log(m2)
    var = log(m2) - 2 * log(m1)

    return {'s': sqrt(var), 'scale': exp(mu)}


# Miscellaneous functions
def param_check(data, param_fix, expt):
    if param_fix is None:
        param_fix = {}
    if expt is None:
        expt = np.full(len(data), 1)

    return param_fix, expt


# MLEs and MMEs for access by distribution name
mles = {'expon': mle_expon,
        'fisk': mle_fisk,
        'gamma': mle_gamma,
        'levy': mle_levy,
        'lognorm': mle_lognorm,
        'pareto': mle_pareto,
        'uniform': mle_uniform}

mmes = {'fisk': mm_fisk,
        'gamma': mm_gamma,
        'lognorm': mm_lognorm}

# Closed-form estimators
cfes = {**mles,
        'fisk': mm_fisk,
        'gamma': mm_gamma}
