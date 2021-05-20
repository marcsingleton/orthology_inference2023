"""Plot example alignments segmented via posterior decoding."""

import json

import matplotlib.pyplot as plt
import scipy.stats as stats
import src.hmm as hmm
import src.draw as draw


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
OGid = '20d6'
msa = load_msa(f'../../ortho_MSA/realign_hmmer/out/{OGid}.mfa')

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
e_dists_rv = {state: stats.betabinom(len(msa)-1, a, b) for state, (a, b) in params['e_dists'].items()}
model = hmm.HMM(params['t_dists'], e_dists_rv, params['start_dist'])

# Decode states and plot
fbs = model.forward_backward(emits)
draw.plot_msa_lines([record[1].upper() for record in msa], [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']])
plt.savefig(f'out/{OGid}.png', bbox_inches='tight')

"""
DEPENDENCIES
../../ortho_MSA/realign_hmmer/realign_hmmer.py
    ../../ortho_MSA/realign_hmmer/out/*.mfa
./segment_hmm1_calc.py
    ./out/model.json
"""