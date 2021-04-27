"""Plot example alignments segmented via posterior decoding."""

import json

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import src.ortho_MSA.hmm as hmm
import src.draw as draw
from scipy.linalg import solve


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
                header = line
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

# Load msa and trim terminal insertions
OGid = '0e80'
msa = load_msa(f'../realign_hmmer/out/{OGid}.mfa')

idx = 0
for j in range(len(msa[0][1])):
    for i in range(len(msa)):
        sym = msa[i][1][j]
        if sym == '.' or sym.islower():
            break
    else:
        idx = j
        break  # if no break exit
msa = [(header, seq[idx:]) for header, seq in msa]

idx = len(msa[0][1])
for j in range(len(msa[0][1]), 0, -1):
    for i in range(len(msa)):
        sym = msa[i][1][j-1]
        if sym == '.' or sym.islower():
            break
    else:
        idx = j
        break  # if no break exit
msa = [(header, seq[:idx]) for header, seq in msa]

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
fb = model.forward_backward(emits)
draw.plot_msa_lines([record[1].upper() for record in msa], [[d['1A'] for d in fb], [d['2'] for d in fb], [d['3'] for d in fb], [d['1B'] for d in fb]])
plt.savefig(f'out/{OGid}.png', bbox_inches='tight')

"""
DEPENDENCIES
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/*.mfa
./hmm_segment_calc.py
    ./out/model.json
"""