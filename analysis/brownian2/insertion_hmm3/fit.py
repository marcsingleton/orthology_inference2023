"""Fit HMM to segmented alignments using conditional maximum likelihood via gradient descent."""

import json
import multiprocessing as mp
import os
from math import comb
from functools import reduce
from itertools import accumulate, product

import scipy.stats as stats
import src.hmm as hmm
from numpy import exp, log
from scipy.special import beta, digamma
from src.utils import read_fasta


# Probability classes and functions
class HMM(hmm.HMM):
    """A class that modifies the default HMM to accept pre-calculated intermediates.

    Discriminative training requires decoding of both states and transitions.
    Both of these calculations rely on the forward and backward variables.
    Since these are the most computationally intensive steps in decoding and
    both are used in state and transition decoding, this modified class accepts
    the forward and backward variables as arguments rather than calculating
    them from scratch each time.
    """
    def forward_backward1(self, emits, fs, ss_f, bs, ss_b):
        """Posterior decoding of states."""
        p = reduce(lambda x, y: x+y, map(log, ss_f))
        ss_f = list(accumulate(map(log, ss_f)))
        ss_b = list(accumulate(map(log, ss_b[::-1])))[::-1]

        fbs = {state: [] for state in self.states}
        for i in range(len(emits)):
            for state, fb in fbs.items():
                fbs[state].append(fs[state][i]*bs[state][i]*exp(ss_f[i]+ss_b[i]-p))

        return fbs

    def forward_backward2(self, emits, fs, ss_f, bs, ss_b):
        """Posterior decoding of transitions."""
        p = reduce(lambda x, y: x+y, map(log, ss_f))
        ss_f = list(accumulate(map(log, ss_f)))
        ss_b = list(accumulate(map(log, ss_b[::-1])))[::-1]

        fbs = {p: [] for p in product(self.states, self.states)}
        for i in range(len(emits)-1):
            for s1, s2 in product(self.states, self.states):
                fbs[(s1, s2)].append(fs[s1][i]*self.t_dists[s1][s2]*self.e_dists[s2].pmf(emits[i+1])*bs[s2][i+1]*exp(ss_f[i]+ss_b[i+1]-p))

        return {p: sum(fb) for p, fb in fbs.items()}

    def joint_likelihood(self, emits, states):
        """Log-likelihood of emission and state sequence."""
        states = [self._state2idx[state] for state in states]
        p, state0 = log(self._e_dists_rv[states[0]].pmf(emits[0]) * self._start_dist_rv.pmf(states[0])), states[0]
        for emit, state in zip(emits[1:], states[1:]):
            p += log(self._e_dists_rv[state].pmf(emit) * self._t_dists_rv[state0].pmf(state))
            state0 = state
        return p


class bernoulli_betabinom_gen:
    def pmf(self, x, p, n, a, b):
        pmf0 = stats.bernoulli.pmf(x[0], p)
        pmf1 = stats.betabinom.pmf(x[1], n, a, b)
        return pmf0 * pmf1

    def rvs(self, p, n, a, b, size=None, random_state=None):
        rvs0 = stats.bernoulli.rvs(p, size=size, random_state=random_state)
        rvs1 = stats.betabinom.rvs(n, a, b, size=size, random_state=random_state)
        if size is None:
            return rvs0, rvs1
        else:
            return list(zip(rvs0, rvs1))

bernoulli_betabinom = bernoulli_betabinom_gen()


class bernoulli_betabinom_frozen:
    def __init__(self, p, n, a, b):
        self._dist = bernoulli_betabinom_gen()
        self.p = p
        self.n = n
        self.a = a
        self.b = b

    def pmf(self, x):
        return self._dist.pmf(x, self.p, self.n, self.a, self.b)

    def rvs(self, size=None):
        return self._dist.rvs(self.p, self.n, self.a, self.b, size=size)


# Gradient functions
def beta_prime(a, b):
    """Return derivative of beta function relative to its first parameter, a."""
    return beta(a, b) * (digamma(a) - digamma(a + b))


def beta_binom_pmf(x, n, a, b):
    """Return pmf of beta-binomial distribution evaluated at x."""
    return comb(n, x) * beta(x + a, n - x + b) / beta(a, b)


def beta_binom_pmf_prime1(x, n, a, b):
    """Return derivative of beta-binomial pmf relative to its first shape parameter, a."""
    return comb(n, x) * (beta_prime(x + a, n - x + b) * beta(a, b) - beta(x + a, n - x + b) * beta_prime(a, b)) / (beta(a, b))**2


def beta_binom_pmf_prime2(x, n, a, b):
    """Return derivative of beta-binomial pmf relative to its second shape parameter, b."""
    return comb(n, x) * (beta_prime(n - x + b, x + a) * beta(a, b) - beta(x + a, n - x + b) * beta_prime(b, a)) / (beta(a, b))**2


def bernoulli_pmf(x, p):
    """Return pmf of Bernoulli distribution evaluated at x."""
    return p**x * (1 - p)**(1 - x)


def bernoulli_pmf_prime(x, p):
    """Return derivative of Bernoulli pmf evaluated at x."""
    return bernoulli_pmf(x, p) * (x/p - (1-x)/(1-p))


# Utility functions
def count_transitions(state_seq, state_set):
    """Return counts of transitions between states."""
    mijs, state0 = {p: 0 for p in product(state_set, state_set)}, state_seq[0]
    for state in state_seq[1:]:
        try:
            mijs[(state0, state)] += 1
        except KeyError:
            mijs[(state0, state)] = 1
        state0 = state
    return mijs


def count_states(state_seq, state_set):
    """Return counts of states."""
    mis = {state: [] for state in state_set}
    for state in state_seq:
        for s in state_set:
            mis[s].append(1 if state == s else 0)
    return mis


def norm_params(t_dists, e_dists):
    """Return parameters as their normalized values."""
    t_dists_norm = {}
    for s1, t_dist in t_dists.items():
        z_sum = sum([exp(z) for z in t_dist.values()])
        t_dists_norm[s1] = {s2: exp(z)/z_sum for s2, z in t_dist.items()}
    e_dists_norm = {}
    for s, (zp, za, zb) in e_dists.items():
        e_dists_norm[s] = 1 / (1 + exp(-zp)), exp(za), exp(zb)
    return t_dists_norm, e_dists_norm


def get_expectations(t_dists_norm, e_dists_norm, start_dist, record):
    """Return record updated with expected values of states and transitions given model parameters."""
    # Instantiate model
    n, state_seq, emit_seq = record['n'], record['state_seq'], record['emit_seq']
    model = HMM(t_dists_norm,
                {state: bernoulli_betabinom_frozen(p, n-1, a, b) for state, (p, a, b) in e_dists_norm.items()},
                start_dist)

    # Get expectations
    fs, ss_f = model.forward(emit_seq)
    bs, ss_b = model.backward(emit_seq)
    record['nis'] = model.forward_backward1(emit_seq, fs, ss_f, bs, ss_b)
    record['nijs'] = model.forward_backward2(emit_seq, fs, ss_f, bs, ss_b)

    # Calculate likelihood
    px = reduce(lambda x, y: x + y, map(log, ss_f))
    pxy = model.joint_likelihood(emit_seq, state_seq)
    record['ll'] = pxy - px

    return record


eta = 1E-4  # Learning rate
epsilon = 1E-2  # Convergence criterion
iter_num = 200  # Max number of iterations
num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])

if __name__ == '__main__':
    # Load regions
    OGid2regions = {}
    state_set = set()
    with open('../config/segments.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            OGid, start, stop, state = line.split()
            if state != '0':  # Skip terminal insertions as actual state
                state_set.add(state)
            try:
                OGid2regions[OGid].append((int(start), int(stop), state))
            except KeyError:
                OGid2regions[OGid] = [(int(start), int(stop), state)]

    # Convert MSAs to records containing state-emissions sequences and other data
    records = []
    for OGid, regions in OGid2regions.items():
        # Load MSA and trim terminal insertions
        msa = read_fasta(f'../../ortho_MSA/realign_hmmer1/out/{OGid}.mfa')
        if regions[-1][2] == '0':
            start, _, _ = regions[-1]
            regions = regions[:-1]
            trim = []
            for header, seq in msa:
                trim.append((header, seq[:start]))
            msa = trim
        if regions[0][2] == '0':
            _, stop, _ = regions[0]
            trim = []
            for header, seq in msa:
                trim.append((header, seq[stop:]))
            msa = trim

            offset = regions[0][1]
            rs = []
            for start, stop, state in regions[1:]:
                rs.append((start-offset, stop-offset, state))
            regions = rs

        # Create emission sequence
        col0 = []
        emit_seq = []
        for j in range(len(msa[0][1])):
            col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
            emit0 = all([c0 == c for c0, c in zip(col0, col)])
            emit1 = sum(col)
            emit_seq.append((emit0, emit1))
            col0 = col

        # Create state sequence
        state_seq = []
        for (start, stop, state) in regions:
            state_seq.extend((stop - start) * [state])

        # Create count dictionaries
        mis = count_states(state_seq, state_set)
        mijs = count_transitions(state_seq, state_set)

        records.append({'OGid': OGid, 'n': len(msa), 'state_seq': state_seq, 'emit_seq': emit_seq,
                        'mis': mis, 'mijs': mijs})

    # Initialize parameters
    with open('../insertion_hmm2/out/model.json') as file:
        params = json.load(file)
    t_dists = params['t_dists']
    e_dists = params['e_dists']
    start_dist = params['start_dist']

    t_dists = {s1: {s2: log(v) for s2, v in t_dist.items()} for s1, t_dist in t_dists.items()}
    e_dists = {s: (-log(1/p-1), log(a), log(b)) for s, (p, a, b) in e_dists.items()}

    # Gradient descent
    j, ll0 = 0, None
    models = []
    while j < iter_num:
        # Calculate expectations and likelihoods
        t_dists_norm, e_dists_norm = norm_params(t_dists, e_dists)
        with mp.Pool(processes=num_processes) as pool:
            records = pool.starmap(get_expectations, [(t_dists_norm, e_dists_norm, start_dist, record) for record in records])

        # Save and report updated parameters
        ll = sum([record['ll'] for record in records])
        models.append((ll, t_dists_norm, e_dists_norm))

        print('ITERATION', j)
        print('ll:', ll)
        print('t_dists_norm:', t_dists_norm)
        print('e_dists_norm:', e_dists_norm)
        print()

        # Check convergence
        if ll0 is not None and abs(ll - ll0) < epsilon:
            break
        ll0 = ll

        # Update parameters
        for record in records:
            # Unpack variables
            n, emit_seq = record['n'], record['emit_seq']
            mis, mijs = record['mis'], record['mijs']
            nis, nijs = record['nis'], record['nijs']

            # Update t_dists
            for s1, t_dist in t_dists.items():
                mn_sum = sum([mijs[(s1, s2)] - nijs[(s1, s2)] for s2 in state_set])
                z_sum = sum([exp(z) for z in t_dist.values()])
                for s2, z in t_dist.items():
                    d = -(mijs[(s1, s2)] - nijs[(s1, s2)] - exp(z)/z_sum * mn_sum)  # Equation 2.17
                    t_dist[s2] -= eta * d

            # Update e_dists
            for s, (zp, za, zb) in e_dists.items():
                p, a, b = 1 / (1 + exp(-zp)), exp(za), exp(zb)
                dzp, dza, dzb = 0, 0, 0
                for i, emit in enumerate(emit_seq):
                    # Equations 2.15 and 2.16 (emission parameter phi only)
                    mn = mis[s][i] - nis[s][i]
                    dzp -= mn / bernoulli_pmf(emit[0], p) * bernoulli_pmf_prime(emit[0], p) * p / (1 + exp(zp))
                    dza -= mn / beta_binom_pmf(emit[1], n-1, a, b) * beta_binom_pmf_prime1(emit[1], n-1, a, b) * a
                    dzb -= mn / beta_binom_pmf(emit[1], n-1, a, b) * beta_binom_pmf_prime2(emit[1], n-1, a, b) * b
                e_dists[s] = (zp - eta * dzp, za - eta * dza, zb - eta * dzb)

        j += 1

    # Save parameters
    if not os.path.exists('out/'):
        os.mkdir('out/')
    with open('out/model.json', 'w') as file:
        _, t_dists, e_dists = max(models)
        json.dump({'t_dists': t_dists, 'e_dists': e_dists, 'start_dist': start_dist}, file)

"""
NOTES
This HMM uses a beta-binomial emission distribution on the gap counts. It also uses a Bernoulli distribution on if the
pattern of gaps is the same as in the previous column. The parameters are trained discriminatively using gradient
descent.

The gradients are calculated using the formulas in:
Krogh A, Riis SK. Hidden Neural Networks. Neural Computation. 11, 541-563. 1999.

DEPENDENCIES
../../ortho_MSA/realign_hmmer1/realign_hmmer1.py
    ../../ortho_MSA/realign_hmmer1/out/*.mfa
../config/segments.tsv
../insertion_hmm2/fit.py
    ../insertion_hmm2/out/model.json
"""