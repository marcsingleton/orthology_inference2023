"""Fit mixture models to rate distributions."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.optimize as opt
import scipy.stats as stats
from math import exp, log, pi, sqrt
from random import random
from scipy.special import digamma


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


def mixture_model(data, dists, params, params_fix=None, weights_init=None, max_iter=250, tol=1E-3, verbose=False):
    # Initialize variables
    weights = np.full(len(dists), 1 / len(dists)) if weights_init is None else weights_init.copy()
    params_fix = [{} for _ in range(len(dists))] if params_fix is None else params_fix
    params_opt = []
    for dist, param, param_fix in zip(dists, params, params_fix):
        if type(param) is dict:
            params_opt.append(param)
        elif type(param) is list:
            sample = data.sample(frac=random())
            mm = mms[dist.name](sample)
            params_opt.append({**mm, **param_fix})

    for n in range(1, max_iter + 1):

        # Expectation
        p = np.stack([dist.pdf(data, **param_fix, **param_opt) for dist, param_fix, param_opt in zip(dists, params_fix, params_opt)], axis=0)  # Apply pdfs for each component and stack results
        expts = p / p.sum(axis=0)  # Normalize result
        weights = expts.sum(axis=1) / expts.sum()

        # Maximization
        ll0 = log_likelihood(data, dists, params_fix, params_opt, weights)
        for dist, param_fix, param_opt, expt in zip(dists, params_fix, params_opt, expts):
            mle = mles[dist.name]  # Get MLE function
            sig = sigs[dist.name](data, param_fix, params_opt, expt)  # Get signature for MLE function
            est = mle(**sig)  # Get updated parameters
            param_opt.update(est)
        ll = log_likelihood(data, dists, params_fix, params_opt, weights)

        # Print output
        if verbose:
            print(n, ll, sep=': ')

        # Test convergence
        if ll - ll0 < tol:
            return param_opt, weights, ll, n, True

    return params_opt, weights, ll, n, False


def mle_expon(data, param_fix=None, expt=None):
    param_fix, expt = param_check(data, param_fix, expt)

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        e = expt.sum()
        ed = (expt * data).sum()
        scale = ed / e

    return {'scale': scale}


def mle_fisk(data, param_fix=None, expt=None, param_opt=None):
    param_fix, expt = param_check(data, param_fix, expt)

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        fisk_scale = create_fisk_scale(data, expt)
        scale = opt.newton(fisk_scale, param_opt['scale'])

    # Shape
    if 'c' in param_fix:
        c = param_fix['c']
    else:
        fisk_shape = create_fisk_shape(data, expt, scale)
        c = opt.newton(fisk_shape, param_opt['c'])

    return {'c': c, 'scale': scale}


def mle_gamma(data, param_fix=None, expt=None, param_opt=None):
    param_fix, expt = param_check(data, param_fix, expt)

    # Shape
    if 'a' in param_fix:
        a = param_fix['a']
    else:
        gamma_shape = create_gamma_shape(data, expt)
        try:
            a = opt.newton(gamma_shape, param_opt['a'])
        except ValueError:  # Catch an error raised by a trial value below zero
            a = opt.brentq(gamma_shape, 1E-9, 2 * param_opt['a'])  # Attempt Brent's method on (1E-6, 2a)

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        scale = (expt * data).sum() / (a * expt.sum())

    return {'a': a, 'scale': scale}


def mle_levy(data, param_fix=None, expt=None):
    param_fix, expt = param_check(data, param_fix, expt)

    # Scale
    if 'scale' in param_fix:
        scale = param_fix['scale']
    else:
        e = expt.sum()
        edivd = (expt / data).sum()
        scale = e / edivd

    return {'scale': scale}


def mle_lognorm(data, param_fix=None, expt=None):
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


def mle_pareto(data, param_fix, expt=None):
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


def mle_uniform(data, param_fix=None):
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


def mm_fisk(data, expt=None):
    expt = np.full(len(data), 1) if expt is None else expt

    # Moments
    logdata = np.log(data)
    m1 = (logdata * expt).sum() / expt.sum()
    m2 = (logdata ** 2 * expt).sum() / expt.sum()

    return {'c': pi / sqrt(3 * (m2 - m1 ** 2)), 'scale': exp(m1)}


def mm_gamma(data, expt=None):
    expt = np.full(len(data), 1) if expt is None else expt

    # Moments
    m1 = (data * expt).sum() / expt.sum()
    m2 = (data ** 2 * expt).sum() / expt.sum()

    return {'a': m1 ** 2 / (m2 - m1 ** 2), 'scale': (m2 - m1 ** 2) / m1}


def mm_lognorm(data, expt=None):
    expt = np.full(len(data), 1) if expt is None else expt

    # Moments
    m1 = (data * expt).sum() / expt.sum()
    m2 = (data ** 2 * expt).sum() / expt.sum()

    # Parameters of transformed lognorm
    mu = 2 * log(m1) - 0.5 * log(m2)
    var = log(m2) - 2 * log(m1)

    return {'s': sqrt(var), 'scale': exp(mu)}


def log_likelihood(data, dists, params_fix, params_opt, weights):
    p = 0
    for dist, param_fix, param_opt, weight in zip(dists, params_fix, params_opt, weights):
        p += weight * dist.pdf(data, **param_fix, **param_opt)  # This step is not log b/c we are summing the contribution of each component
    return np.log(p).sum()


def param_check(data, param_fix, expt):
    if param_fix is None:
        param_fix = {}
    if expt is None:
        expt = np.full(len(data), 1)

    return param_fix, expt


def sig1(data, param_fix, param_opt, expt):
    return {'data': data, 'param_fix': param_fix}


def sig2(data, param_fix, param_opt, expt):
    return {'data': data, 'param_fix': param_fix, 'expt': expt}


def sig3(data, param_fix, param_opt, expt):
    return {'data': data, 'param_fix': param_fix, 'param_opt': param_opt, 'expt': expt}


# Maximum likelihood estimators for EM algorithm
mles = {'expon': mle_expon,
        'fisk': mle_fisk,
        'gamma': mle_gamma,
        'levy': mle_levy,
        'lognorm': mle_lognorm,
        'pareto': mle_pareto,
        'uniform': mle_uniform}

# Method of moments estimators for initialization
mms = {**mles,
       'fisk': mm_fisk,
       'gamma': mm_gamma}

# Signatures for MLEs
sigs = {'expon': sig2,
        'fisk': sig3,
        'gamma': sig3,
        'levy': sig2,
        'lognorm': sig2,
        'pareto': sig2,
        'uniform': sig1}

# Input variables
path = r'C:\Users\marcs\Documents\Berkeley\Projects\IDREvoDevo\analysis\brownian\pic_calc\pics.tsv'
lt = 32

# Read data and filter
pics = pd.read_csv(path, sep='\t', index_col=list(range(3)))
pics_lt = pics[(pics.index.get_level_values('min_length') >= lt) & (~pics.index.get_level_values('ordered').array.astype(bool))]

# Calculate rates
rates = (pics_lt ** 2).groupby('block_id').mean()
feature = 'SCD'
data = rates.loc[rates[feature] != 0, feature].sort_values()

# Initialization
dists = [stats.lognorm, stats.lognorm, stats.lognorm]
params = [['s', 'scale'], ['s', 'scale'], ['s', 'scale']]

# Plot results
plot = data[:-50]
plt.hist(plot, bins=500, density=True)
ylim = plt.ylim()
x = np.linspace(plot.min(), plot.max(), 2000)
total = 0
for dist, param, weight in zip(dists, params_opt, weights):
    y = weight * dist.pdf(x, **param)
    total += y
    plt.plot(x, weight * dist.pdf(x, **param), label=dist.name)
plt.plot(x, total, label='total')
plt.legend()
plt.ylim(ylim)
