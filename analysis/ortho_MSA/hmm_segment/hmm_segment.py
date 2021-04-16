"""Fit HMM to segmented alignments."""

import json
import os

import numpy as np
import scipy.optimize as opt
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
with open('segments.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, state = line.split()
        if state != '0':  # Skip terminal insertions as actual state
            states.add(state)
        try:
            OGid2regions[OGid].append((int(start), int(stop), state))
        except KeyError:
            OGid2regions[OGid] = [(int(start), int(stop), state)]

# Get counts
t_counts = {state: {s: 1 for s in states} for state in states}  # Add pseudocounts to transitions
e_counts = {state: {} for state in states}
start_count = {state: 1 for state in states}  # Add pseudocounts to starts
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

    # Create state and emission sequences
    n = len(msa)
    gaps = len(msa[0][1]) * [0]
    for _, seq in msa:
        for i, sym in enumerate(seq):
            if sym in ['-', '.']:
                gaps[i] += 1
    states = []
    for (start, stop, state) in regions:
        states.extend((stop-start)*[state])

    # Count transitions and emissions
    state0 = None
    for i, (gap, state) in enumerate(zip(gaps, states)):
        try:
            e_counts[state][n].append(gap)
        except KeyError:
            e_counts[state][n] = [gap]

        if state0 is None:
            start_count[state] += 1
        else:
            t_counts[state0][state] = t_counts[state0].get(state, 0) + 1
        state0 = state

# Convert counts to distributions
t_dists = {}
for state, t_count in t_counts.items():
    total = sum(t_count.values())
    t_dists[state] = {s: count/total for s, count in t_count.items()}

e_dists = {}
for state, e_count in e_counts.items():
    # Compile counts from different models
    lls = []
    a, b, w = [], [], []
    for n, count in e_count.items():
        lls.append(create_betabinom_likelihood(count, n-1))
        ests = mm_betabinom(np.asarray(count), n-1)
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

total = sum(start_count.values())
start_dist = {s: count/total for s, count in start_count.items()}

# Save parameters
if not os.path.exists('out/'):
    os.mkdir('out/')
with open('out/model.json', 'w') as file:
    json.dump({'t_dists': t_dists, 'e_dists': e_dists, 'start_dist': start_dist}, file)

"""
DEPENDENCIES
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/*.mfa
./segments.tsv
"""