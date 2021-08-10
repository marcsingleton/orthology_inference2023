"""Plot example alignments segmented via posterior decoding."""

import json

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import src.hmm as hmm
import src.draw as draw
from scipy.linalg import solve
from src.brownian2.trim import trim_terminals


class ar1_betabinom_gen:
    def pmf(self, i, j, n, a0, b0, a1, b1):
        pmfs = [stats.betabinom.pmf(j - k, n - i, a0, b0) * stats.betabinom.pmf(k, i, a1, b1)
                for k in range(max(0, i + j - n), min(i, j) + 1)]
        return sum(pmfs)

    def rvs(self, i, n, a0, b0, a1, b1, size=None, random_state=None):
        js = list(range(n+1))
        ps = [self.pmf(i, j, n, a0, b0, a1, b1) for j in js]
        rvs = stats.rv_discrete(values=(js, ps)).rvs(size=size, random_state=random_state)
        return rvs


ar1_betabinom = ar1_betabinom_gen()


class ar1_betabinom_frozen:
    def __init__(self, n, a0, b0, a1, b1):
        self.n = n
        self.a0 = a0
        self.b0 = b0
        self.a1 = a1
        self.b1 = b1

        dist = {}
        for i in range(n+1):
            for j in range(n+1):
                dist[(i, j)] = ar1_betabinom.pmf(i, j, n, a0, b0, a1, b1)
        self.dist = dist

    def pmf(self, i, j):
        return self.dist[(i, j)]

    def rvs(self, i, size=None, random_state=None):
        js = list(range(self.n+1))
        ps = [self.pmf(i, j) for j in js]
        rvs = stats.rv_discrete(values=(js, ps)).rvs(size=size, random_state=random_state)
        return rvs


def get_stationary_dist(n, a0, b0, a1, b1):
    a = np.zeros((n+1, n+1))
    for i in range(n+1):
        for j in range(n+1):
            a[i, j] = ar1_betabinom.pmf(i, j, n, a0, b0, a1, b1)

    a = a.transpose() - np.identity(n+1)
    a[n] = 1
    b = np.zeros(n+1)
    b[n] = 1
    pi = solve(a, b)
    pi = np.where(pi < 0, 0, pi)  # Remove negative values
    pi = pi / pi.sum()  # Rescale to sum to 1

    return stats.rv_discrete(values=(np.arange(n+1), pi))


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line.rstrip()
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


# Load model parameters
with open('out/model.json') as file:
    params = json.load(file)

# Load OGids
OGids = set()
with open('../config/segments.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid = line.split()[0]
        OGids.add(OGid)

# Load msa and trim terminal insertions
for OGid in OGids:
    msa = trim_terminals(load_msa(f'../../ortho_MSA/realign_hmmer/out/{OGid}.mfa'))

    # Create emission sequence
    emits = []
    for j in range(len(msa[0][1])):
        col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
        emits.append(sum(col))

    # Instantiate model
    e_dists_rv = {state: ar1_betabinom_frozen(len(msa)-1, a0, b0, a1, b1) for state, (a0, b0, a1, b1) in params['e_dists'].items()}
    start_e_dists_rv = {state: get_stationary_dist(len(msa)-1, a0, b0, a1, b1) for state, (a0, b0, a1, b1) in params['e_dists'].items()}
    model = hmm.ARHMM(params['t_dists'], e_dists_rv, params['start_t_dist'], start_e_dists_rv)

    # Decode states and plot
    fbs = model.forward_backward(emits)
    draw.plot_msa_lines([record[1].upper() for record in msa], [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']], figsize=(15, 6))
    plt.savefig(f'out/{OGid}_wide.png', bbox_inches='tight')
    plt.close()

    draw.plot_msa_lines([record[1].upper() for record in msa], [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']], figsize=(8, 8))
    plt.savefig(f'out/{OGid}_tall.png', bbox_inches='tight')
    plt.close()

"""
DEPENDENCIES
../../ortho_MSA/realign_hmmer/realign_hmmer.py
    ../../ortho_MSA/realign_hmmer/out/*.mfa
./fit.py
    ./out/model.json
"""