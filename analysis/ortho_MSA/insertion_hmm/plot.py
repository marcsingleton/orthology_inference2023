"""Plot example alignments segmented via posterior decoding."""

import json
import re

import homomorph
import matplotlib.pyplot as plt
import numpy as np
import skbio
import utils
from src.draw import plot_msa_data
from src.utils import read_fasta

# Load model parameters
with open('out/model.json') as file:
    params = json.load(file)


OGid2labels = {}
with open('labels.tsv') as file:
    file.readline()
    for line in file:
        OGid, start, stop, state = line.rstrip('\n').split('\t')
        if OGid in OGid2labels:
            OGid2labels[OGid][state].append((int(start), int(stop)))
        else:
            labels = {'0': [], '1A': [], '1B': [], '2': [], '3': []}
            labels[state].append((int(start), int(stop)))
            OGid2labels[OGid] = labels

# Load tree
tree_template = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree_template.tips())}

# Plot alignments
for OGid, labels in OGid2labels.items():
    # Load MSA
    msa = read_fasta(f'../../ortho_MSA/realign_hmmer/out/mafft/{OGid}.afa')
    msa = [(re.search(r'spid=([a-z]+)', header).group(1), seq) for header, seq in msa]

    # Create emission sequence
    col0 = []
    emit_seq = []
    for j in range(len(msa[0][1])):
        col = [1 if msa[i][1][j] in ['-', '.'] else 0 for i in range(len(msa))]
        emit0 = all([c0 == c for c0, c in zip(col0, col)])
        emit_seq.append((emit0, j))  # The tree probabilities are pre-calculated, so emission value is its index
        col0 = col

    # Load tree and convert to vectors at tips
    tree = tree_template.deepcopy().shear([spid for spid, _ in msa])
    tips = {tip.name: tip for tip in tree.tips()}
    for spid, seq in msa:
        tip = tips[spid]
        conditional = np.zeros((2, len(seq)))
        for j, sym in enumerate(seq):
            if sym in ['-', '.']:
                conditional[0, j] = 1
            else:
                conditional[1, j] = 1
        tip.conditional = conditional

    # Instantiate model
    e_dists_rv = {}
    for state, (p, pi, q0, q1) in params['e_dists'].items():
        array = utils.get_tree_probability(tree, pi, q0, q1)
        e_dists_rv[state] = utils.BinomialArrayRV(p, array)
    model = homomorph.HMM(params['t_dists'], e_dists_rv, params['start_dist'])

    # Plot labels
    msa = [seq for _, seq in sorted(msa, key=lambda x: tip_order[x[0]])]  # Re-order sequences and extract seq only

    lines = {}
    for state in ['1A', '1B', '2', '3']:
        line = np.zeros(len(msa[0][1]))
        for start, stop in labels[state]:
            line[start:stop] = 1
        lines[state] = line

    plot_msa_data(msa, [lines['1A'], lines['2'], lines['3'], lines['1B']], figsize=(15, 6))
    plt.savefig(f'out/{OGid}_wide_labels.png', bbox_inches='tight')
    plt.close()

    plot_msa_data(msa, [lines['1A'], lines['2'], lines['3'], lines['1B']], figsize=(15, 6))
    plt.savefig(f'out/{OGid}_tall_labels.png', bbox_inches='tight')
    plt.close()

    # Decode states and plot
    fbs = model.forward_backward(emit_seq)

    plot_msa_data(msa, [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']], figsize=(15, 6))
    plt.savefig(f'out/{OGid}_wide.png', bbox_inches='tight')
    plt.close()

    plot_msa_data(msa, [fbs['1A'], fbs['2'], fbs['3'], fbs['1B']], figsize=(8, 8))
    plt.savefig(f'out/{OGid}_tall.png', bbox_inches='tight')
    plt.close()

"""
DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/mafft/*.afa
./labels.tsv
./fit.py
    ./out/model.json
"""