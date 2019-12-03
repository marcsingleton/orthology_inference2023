"""Fit mixture models with arbitrary components to data."""

import numpy as np
from random import random
from estimate import cfes, mles, mmes


def log_likelihood(data, dists, params_fix, params, weights):
    p = 0
    for dist, param_fix, param, weight in zip(dists, params_fix, params, weights):
        p += weight * dist.pdf(data, **param_fix, **param)  # This step is not log b/c we are summing the contribution of each component
    return np.log(p).sum()


class MixtureModel:

    def __init__(self, dists, params, params_fix=None, weights=None):
        # Model parameters
        self.dists = dists.copy()
        self.params = params.copy()
        self.params_fix = [{} for _ in range(len(dists))] if params_fix is None else params_fix.copy()
        self.weights = np.full(len(dists), 1 / len(dists)) if weights is None else weights.copy()

        # Fit metadata
        self.ll = 0
        self.converged = False
        self.n_iter = 0

    def fit(self, data, max_iter=250, tol=1E-3, verbose=False):
        # Use temporary values to preserve originals in case of error
        params_opt = []
        for dist, param, param_fix in zip(self.dists, self.params, self.params_fix):
            if type(param) is dict:  # Copy values if initial values provided
                params_opt.append(param.copy())
            elif type(param) is list:  # Estimate from random sample otherwise
                sample = data.sample(frac=random())
                cfe = cfes[dist.name]  # Get closed-form estimator
                param_init = cfe(sample, param_fix=param_fix)  # Get initial estimates
                params_opt.append({**param_init, **param_fix})
        weights_opt = self.weights.copy()

        for i in range(1, max_iter + 1):
            ll0 = log_likelihood(data, self.dists, self.params_fix, params_opt, weights_opt)

            # Expectation
            model_params = zip(self.dists, self.params_fix, params_opt, weights_opt)
            p = np.stack([weight_opt * dist.pdf(data, **param_fix, **param_opt) for dist, param_fix, param_opt, weight_opt in model_params], axis=0)  # Apply pdfs for each component and stack results
            expts = p / p.sum(axis=0)  # Normalize result
            weights_opt = expts.sum(axis=1) / expts.sum()

            # Maximization
            for dist, param_fix, param_opt, expt in zip(self.dists, self.params_fix, params_opt, expts):
                mle = mles[dist.name]  # Get MLE function
                opt = mle(data, param_fix=param_fix, param=param_opt, expt=expt)  # Get updated parameters
                param_opt.update(opt)
            ll = log_likelihood(data, self.dists, self.params_fix, params_opt, weights_opt)

            # Print output
            if verbose:
                print(i, ll, sep=': ')

            # Test convergence
            if ll - ll0 < tol:
                self.converged = True
                break

        self.params = params_opt
        self.weights = weights_opt.tolist()
        self.ll = ll
        self.n_iter = i

    def pdf(self, x, comp=None):
        model_params = zip(self.dists, self.params_fix, self.params, self.weights)
        if comp is None:
            p = np.stack([weight * dist.pdf(x, **param_fix, **param) for dist, param_fix, param, weight in model_params], axis=0)  # Apply pdfs for each component and stack results
            return p.sum(axis=0)
        else:
            dist, param_fix, param, weight = list(model_params)[comp]
            p = weight * dist.pdf(x, **param_fix, **param)
            return p

    def post_prob(self, x):
        model_params = zip(self.dists, self.params_fix, self.params, self.weights)
        p = np.stack([weight * dist.pdf(x, **param_fix, **param) for dist, param_fix, param, weight in model_params], axis=0)  # Apply pdfs for each component and stack results
        return p / p.sum(axis=0)  # Normalize result
