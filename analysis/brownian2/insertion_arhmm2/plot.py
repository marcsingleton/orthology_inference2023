"""Plot example alignments segmented via posterior decoding."""

import json
import re

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import skbio
import src.hmm as hmm
import src.draw as draw
from scipy.linalg import solve
from src.brownian2.trim import trim_terminals
from src.utils import read_fasta


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


class ar1_binom_gen:
    def pmf(self, i, j, n, p0, p1):
        pmfs = [stats.binom.pmf(j - k, n - i, p0) * stats.betabinom.pmf(k, i, p1)
                for k in range(max(0, i + j - n), min(i, j) + 1)]
        return sum(pmfs)

    def rvs(self, i, n, p0, p1, size=None, random_state=None):
        js = list(range(n+1))
        ps = [self.pmf(i, j, n, p0, p1) for j in js]
        rvs = stats.rv_discrete(values=(js, ps)).rvs(size=size, random_state=random_state)
        return rvs


ar1_binom = ar1_binom_gen()


class ar1_binom_frozen:
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


class betabinom_frozen:
    def __init__(self, n, a, b):
        self.n = n
        self.a = a
        self.b = b
        self.dist = stats.betabinom(n, a, b)

    def pmf(self, i, j):
        return self.dist.pmf(j)

    def rvs(self, i, size=None, random_state=None):
        return self.dist.rvs(size=size, random_state=random_state)


class rv_discrete_frozen:
    def __init__(self, values):
        self.values = values
        self.dist = stats.rv_discrete(values=values)

    def pmf(self, i, j):
        return self.dist.pmf(j)

    def rvs(self, i, size=None, random_state=None):
        return self.dist.rvs(size=size, random_state=random_state)


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

# Load tree
tree = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree.tips())}

# Load msa and trim terminal insertions
for OGid in OGids:
    msa = trim_terminals(read_fasta(f'../../ortho_MSA/realign_hmmer1/out/{OGid}.mfa'))
    msa = [(re.search(r'spid=([a-z]+)', header).group(1), seq) for header, seq in msa]

    # Create emission sequence
    emits = []
    for j in range(len(msa[0][1])):
        col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
        emits.append(sum(col))

    # Instantiate model
    e_dists_rv = {}
    start_e_dists_rv = {}
    for state in ['1A', '3']:
        a, b = params['e_dists'][state]
        e_dists_rv[state] = betabinom_frozen(len(msa)-1, a, b)
        start_e_dists_rv[state] = stats.betabinom(len(msa)-1, a, b)
    for state in ['1B', '2']:
        a0, b0, a1, b1 = params['e_dists'][state]
        e_dists_rv[state] = ar1_betabinom_frozen(len(msa)-1, a0, b0, a1, b1)
        start_e_dists_rv[state] = get_stationary_dist(len(msa)-1, a0, b0, a1, b1)
    model = hmm.ARHMM(params['t_dists'], e_dists_rv, params['start_t_dist'], start_e_dists_rv)

    # Decode states and plot
    fbs = model.forward_backward(emits)
    msa = [seq.upper() for _, seq in sorted(msa, key=lambda x: tip_order[x[0]])]  # Re-order sequences and extract seq only

    draw.plot_msa_lines(msa, [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']], figsize=(15, 6))
    plt.savefig(f'out/{OGid}_wide.png', bbox_inches='tight')
    plt.close()

    draw.plot_msa_lines(msa, [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']], figsize=(8, 8))
    plt.savefig(f'out/{OGid}_tall.png', bbox_inches='tight')
    plt.close()

"""
DEPENDENCIES
../../ortho_MSA/realign_hmmer1/realign_hmmer1.py
    ../../ortho_MSA/realign_hmmer/out/*.mfa
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
../config/segments.tsv
./fit.py
    ./out/model.json
"""