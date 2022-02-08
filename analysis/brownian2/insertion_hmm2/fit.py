"""Fit HMM to segmented alignments."""

import json
import os

import numpy as np
import scipy.optimize as opt
from numpy import log
from src.utils import read_fasta


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
e_counts = {state: ([], {}) for state in states}
start_count = {state: 1 for state in states}

# Get observed counts
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
    emits = []
    for j in range(len(msa[0][1])):
        col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
        emit0 = all([c0 == c for c0, c in zip(col0, col)])
        emit1 = sum(col)
        emits.append((emit0, emit1))
        col0 = col

    # Create state sequence
    states = []
    for (start, stop, state) in regions:
        states.extend((stop-start)*[state])

    # Count emissions and transitions
    n, state0 = len(msa), None
    for i, ((emit0, emit1), state) in enumerate(zip(emits, states)):
        e_counts[state][0].append(emit0)
        try:
            e_counts[state][1][n].append(emit1)
        except KeyError:
            e_counts[state][1][n] = [emit1]

        if state0 is not None:
            t_counts[state0][state] = t_counts[state0].get(state, 0) + 1
        else:
            start_count[state] += 1
        state0 = state

# Convert counts to distributions
t_dists = {}
for state, t_count in t_counts.items():
    total = sum(t_count.values())
    t_dists[state] = {s: count/total for s, count in t_count.items()}

e_dists = {}
for state, (e_count0, e_count1) in e_counts.items():
    # Compile counts from different ns
    lls = []
    a, b, w = [], [], []
    for n, count1 in e_count1.items():
        lls.append(create_betabinom_likelihood(count1, n-1))
        ests = mm_betabinom(np.asarray(count1), n-1)
        if ests['a'] > 0 and ests['b'] > 0:  # Under-dispersed data can yield negative mm estimates
            a.append(len(count1)*ests['a'])
            b.append(len(count1)*ests['b'])
            w.append(len(count1))

    # Merge likelihoods and initial parameter estimates
    ll = lambda x: -sum([ll(x) for ll in lls])
    a = sum(a) / sum(w) if a else 0.5
    b = sum(b) / sum(w) if b else 0.5
    v = 1 / (a + b)
    u = a * v

    # Compute ML estimates
    result = opt.minimize(ll, [u, v], bounds=[(0, 1), (0, None)])
    u, v = result.x
    a = u / v
    b = (1 - u) / v
    p = sum(e_count0) / len(e_count0)
    e_dists[state] = (p, a, b)

    # Print results
    print(f'STATE {state}')
    print(result)

total = sum(start_count.values())
start_dist = {state: count/total for state, count in start_count.items()}

# Save parameters
if not os.path.exists('out/'):
    os.mkdir('out/')
with open('out/model.json', 'w') as file:
    json.dump({'t_dists': t_dists, 'e_dists': e_dists, 'start_dist': start_dist}, file)

"""
NOTES
This HMM uses a beta-binomial emission distribution on the gap counts. It also uses a Bernoulli distribution on if the
pattern of gaps is the same as in the previous column.

DEPENDENCIES
../../ortho_MSA/realign_hmmer1/realign_hmmer1.py
    ../../ortho_MSA/realign_hmmer1/out/*.mfa
../config/segments.tsv
"""