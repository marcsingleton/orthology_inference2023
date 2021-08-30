"""Plot example alignments segmented via posterior decoding."""

import json

import matplotlib.pyplot as plt
import scipy.stats as stats
import src.hmm as hmm
import src.draw as draw
from src.brownian2.trim import trim_terminals


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

# Plot alignments
for OGid in OGids:
    # Load msa and trim terminal insertions
    msa = trim_terminals(load_msa(f'../../ortho_MSA/realign_hmmer/out/{OGid}.mfa'))

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
    draw.plot_msa_lines([seq.upper() for _, seq in msa], [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']], figsize=(15, 6))
    plt.savefig(f'out/{OGid}_wide.png', bbox_inches='tight')
    plt.close()

    draw.plot_msa_lines([seq.upper() for _, seq in msa], [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']], figsize=(8, 8))
    plt.savefig(f'out/{OGid}_tall.png', bbox_inches='tight')
    plt.close()

"""
DEPENDENCIES
../../ortho_MSA/realign_hmmer/realign_hmmer.py
    ../../ortho_MSA/realign_hmmer/out/*.mfa
./fit.py
    ./out/model.json
"""