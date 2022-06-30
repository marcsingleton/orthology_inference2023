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

ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'

# Load model parameters
with open('out/model.json') as file:
    model_json = json.load(file)

# Load labels
OGid2labels = {}
state_set = set()
with open('labels.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        OGid, start, stop, label = fields['OGid'], int(fields['start']), int(fields['stop']), fields['label']
        state_set.add(label)
        try:
            OGid2labels[OGid].append((start, stop, label))
        except KeyError:
            OGid2labels[OGid] = [(start, stop, label)]

# Load tree
tree_template = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree_template.tips())}

# Plot alignments
for OGid, labels in OGid2labels.items():
    # Load MSA
    msa = []
    for header, seq in read_fasta(f'../realign_hmmer/out/mafft/{OGid}.afa'):
        spid = re.search(spid_regex, header).group(1)
        msa.append({'spid': spid, 'seq': seq})

    # Create emission sequence
    col0 = []
    emit_seq = []
    for j in range(len(msa[0]['seq'])):
        col = [1 if msa[i]['seq'][j] in ['-', '.'] else 0 for i in range(len(msa))]
        emit0 = sum([c0 == c for c0, c in zip(col0, col)])
        emit_seq.append(emit0)  # The tree probabilities are pre-calculated, so emission value is its index
        col0 = col
    emit_seq = np.array(emit_seq)

    # Load tree and convert to vectors at tips
    tree = tree_template.deepcopy().shear([record['spid'] for record in msa])
    tips = {tip.name: tip for tip in tree.tips()}
    for record in msa:
        spid, seq = record['spid'], record['seq']
        conditional = np.zeros((2, len(seq)))
        for j, sym in enumerate(seq):
            if sym in ['-', '.']:
                conditional[0, j] = 1
            else:
                conditional[1, j] = 1
        tip = tips[spid]
        tip.conditional = conditional

    # Instantiate model
    e_dists_rv = {}
    for s, params in model_json['e_dists'].items():
        a, b, pi, q0, q1, r = params
        array1 = utils.get_betabinom_pmf(emit_seq, len(msa), a, b)
        array2 = utils.get_tree_pmf(tree, pi, q0, q1, r)
        e_dists_rv[s] = utils.ArrayRV(array1 * array2)
    model = homomorph.HMM(model_json['t_dists'], e_dists_rv, model_json['start_dist'])

    # Make plotting parameters
    msa = sorted(msa, key=lambda x: tip_order[x['spid']])  # Re-order sequences
    data_labels = ['1A', '1B', '2', '3']
    data_colors = ['C0', 'C3', 'C1', 'C2']

    kwargs_wide = {'height_ratio': 0.5, 'hspace': 0.2, 'data_max': 1.1, 'data_min': -0.1, 'data_labels': data_labels, 'data_colors': data_colors,
                   'msa_legend': True, 'legend_kwargs': {'bbox_to_anchor': (0.945, 0.5), 'loc': 'center left', 'fontsize': 8,
                                                         'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1}}
    adjust_wide = {'left': 0.015, 'bottom': 0.01, 'right': 0.94, 'top': 0.99}
    kwargs_tall = {'height_ratio': 0.5, 'hspace': 0.2, 'data_max': 1.1, 'data_min': -0.1, 'data_labels': data_labels, 'data_colors': data_colors,
                   'msa_legend': True, 'legend_kwargs': {'bbox_to_anchor': (0.90, 0.5), 'loc': 'center left', 'fontsize': 8,
                                                         'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1}}
    adjust_tall = {'left': 0.025, 'bottom': 0.01, 'right': 0.89, 'top': 0.99}

    # Plot labels
    lines = {s: np.zeros(len(msa[0]['seq'])) for s in state_set}
    for start, stop, label in labels:
        lines[label][start:stop] = 1
    data = [lines[label] for label in data_labels]

    plot_msa_data([record['seq'] for record in msa], data, figsize=(15, 6), **kwargs_wide)
    plt.subplots_adjust(**adjust_wide)
    plt.savefig(f'out/{OGid}_wide_labels.png')
    plt.close()

    plot_msa_data([record['seq'] for record in msa], data, figsize=(8, 8), **kwargs_tall)
    plt.subplots_adjust(**adjust_tall)
    plt.savefig(f'out/{OGid}_tall_labels.png')
    plt.close()

    # Decode states and plot
    emit_seq = list(range(len(emit_seq)))  # Everything is pre-calculated, so emit_seq is the emit index
    fbs = model.forward_backward(emit_seq)
    data = [fbs[label] for label in data_labels]

    plot_msa_data([record['seq'] for record in msa], data, figsize=(15, 6), **kwargs_wide)
    plt.subplots_adjust(**adjust_wide)
    plt.savefig(f'out/{OGid}_wide.png')
    plt.close()

    plot_msa_data([record['seq'] for record in msa], data, figsize=(8, 8), **kwargs_tall)
    plt.subplots_adjust(**adjust_tall)
    plt.savefig(f'out/{OGid}_tall.png')
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