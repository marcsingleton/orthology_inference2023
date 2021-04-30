"""Fit HMM to segmented alignments."""

import json
import os

import numpy as np
import scipy.optimize as opt
import scipy.stats as stats
from numpy import log


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


def create_ar1_betabinom_likelihood(data, n):

    def betabinom_likelihood(x):
        a0, b0, a1, b1 = x
        cache = {}
        lls = []
        for i, j in data:
            try:
                lls.append(cache[(i, j)])
            except KeyError:
                ll = log(ar1_betabinom_pmf(i, j, n, a0, b0, a1, b1))
                cache[(i, j)] = ll
                lls.append(ll)
        return sum(lls)

    return betabinom_likelihood


def ar1_betabinom_pmf(i, j, n, a0, b0, a1, b1):
    pmfs = [stats.betabinom.pmf(j-k, n-i, a0, b0)*stats.betabinom.pmf(k, i, a1, b1)
            for k in range(max(0, i+j-n), min(i, j)+1)]
    return sum(pmfs)


def create_betabinom_likelihood(data, n):

    def betabinom_likelihood(x):
        u, v = x
        s = 0
        for x in data:
            s += sum([log(u+i*v) for i in range(x)])
            s += sum([log(1-u+i*v) for i in range(n-x)])
        s -= len(data)*sum([log(1+i*v) for i in range(n)])

        return s

    return betabinom_likelihood


def mm_betabinom(data, n):
    ests = {}

    # Moments
    m1 = data.sum() / len(data)
    m2 = (data ** 2).sum() / len(data)

    # Estimators
    ests['a'] = (n*m1-m2) / (n*(m2/m1-m1-1)+m1)
    ests['b'] = (n-m1)*(n-m2/m1) / (n*(m2/m1-m1-1)+m1)

    return ests


# Load regions
OGid2regions = {}
states = set()
with open('../config/segments.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, state = line.split()
        if state != '0':  # Skip terminal insertions as actual state
            states.add(state)
        try:
            OGid2regions[OGid].append((int(start), int(stop), state))
        except KeyError:
            OGid2regions[OGid] = [(int(start), int(stop), state)]

# Initialize counts with pseudocounts
t_counts = {state: {s: 1 for s in states} for state in states}
e_counts = {state: {} for state in states}
start_t_count = {state: 1 for state in states}

# Get observed counts
for OGid, regions in OGid2regions.items():
    # Load MSA and trim terminal insertions
    msa = load_msa(f'../realign_hmmer/out/{OGid}.mfa')
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
    emits = []
    for j in range(len(msa[0][1])):
        col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
        emits.append(sum(col))

    # Create state sequence
    states = []
    for (start, stop, state) in regions:
        states.extend((stop-start)*[state])

    # Count emissions and transitions
    n, state0, emit0 = len(msa), None, None
    for i, (emit, state) in enumerate(zip(emits, states)):
        if emit0 is not None:
            try:
                e_counts[state][n].append((emit0, emit))
            except KeyError:
                e_counts[state][n] = [(emit0, emit)]

        if state0 is not None:
            t_counts[state0][state] = t_counts[state0].get(state, 0) + 1
        else:
            start_t_count[state] += 1
        state0 = state
        emit0 = emit

# Convert counts to distributions
t_dists = {}
for state, t_count in t_counts.items():
    total = sum(t_count.values())
    t_dists[state] = {s: count/total for s, count in t_count.items()}

e_dists = {}
for state in ['1A', '3']:
    # Compile counts from different models
    lls = []
    a, b, w = [], [], []
    for n, count in e_counts[state].items():
        lls.append(create_betabinom_likelihood([c[1] for c in count], n-1))
        ests = mm_betabinom(np.asarray([c[1] for c in count]), n-1)
        if ests['a'] > 0 and ests['b'] > 0:  # Under-dispersed data can yield negative mm estimates
            a.append(len(count)*ests['a'])
            b.append(len(count)*ests['b'])
            w.append(len(count))

    # Merge likelihoods and initial parameter estimates
    ll = lambda x: -sum([ll(x) for ll in lls])
    a = sum(a) / sum(w)
    b = sum(b) / sum(w)
    v = 1 / (a + b)
    u = a * v

    # Compute ML estimates
    result = opt.minimize(ll, [u, v], bounds=[(0, 1), (0, None)])
    u, v = result.x
    a = u / v
    b = (1 - u) / v
    e_dists[state] = (a, b)

    # Print results
    print(f'STATE {state}')
    print(result)
for state in ['1B', '2']:
    # Compile counts from different ns and merge likelihoods
    lls = []
    for n, count in e_counts[state].items():
        lls.append(create_ar1_betabinom_likelihood(count, n-1))
    ll = lambda x: -sum([ll(x) for ll in lls])

    # Compute ML estimates
    result = opt.minimize(ll, 4*[1], bounds=4*[(0, None)], method='Powell')
    a0, b0, a1, b1 = result.x
    e_dists[state] = (a0, b0, a1, b1)

    # Print results
    print(f'STATE {state}')
    print(result)

total = sum(start_t_count.values())
start_t_dist = {s: count/total for s, count in start_t_count.items()}

# Save parameters
if not os.path.exists('out/'):
    os.mkdir('out/')
with open('out/model.json', 'w') as file:
    json.dump({'t_dists': t_dists, 'e_dists': e_dists, 'start_t_dist': start_t_dist}, file)

"""
DEPENDENCIES
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/*.mfa
./segments.tsv
"""