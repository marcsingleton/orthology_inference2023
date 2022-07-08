"""Plot example alignments segmented via posterior decoding."""

import os
import json
import re

import homomorph
import matplotlib.pyplot as plt
import numpy as np
import skbio
from src.ortho_MSA import utils
from matplotlib.lines import Line2D
from src.draw import plot_msa_data
from src.utils import read_fasta

spid_regex = r'spid=([a-z]+)'
state_labels = ['1A', '1B', '2', '3']
state_colors = ['C0', 'C3', 'C1', 'C2']

tree_template = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree_template.tips())}

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

if set(state_labels) != state_set:
    raise RuntimeError('state_labels is not equal to state_set')

# Load history and model parameters
with open('out/history.json') as file:
    history = json.load(file)
with open('out/model.json') as file:
    model_json = json.load(file)

# Plot loss curve
xs = [record['iter_num'] for record in history]
ys = [record['ll'] for record in history]
plt.plot(xs, ys)
plt.xlabel('Iteration')
plt.ylabel('Conditional log-likelihood')
plt.savefig('out/line_ll-iter.png')
plt.close()

# Plot model parameters
fig, axs = plt.subplots(3, 1)
for label, color in zip(state_labels, state_colors):
    for ax, param in zip(axs, ['pi', 'q0', 'q1']):
        xs = [record['iter_num'] for record in history]
        ys = [record['e_dists_norm'][label][param] for record in history]
        ax.plot(xs, ys, label=label, color=color)
        ax.set_ylabel(param)
axs[2].set_xlabel('Iteration')
handles = [Line2D([], [], label=label, color=color) for label, color in zip(state_labels, state_colors)]
fig.legend(handles=handles, bbox_to_anchor=(0.875, 0.5), loc='center left')
plt.subplots_adjust(right=0.875)
plt.savefig('out/line_rate-iter.png')
plt.close()

fig, axs = plt.subplots(2, 1)
for label, color in zip(state_labels, state_colors):
    for ax, param in zip(axs, ['a', 'b']):
        xs = [record['iter_num'] for record in history]
        ys = [record['e_dists_norm'][label][param] for record in history]
        ax.plot(xs, ys, label=label, color=color)
        ax.set_ylabel(param)
axs[1].set_xlabel('Iteration')
handles = [Line2D([], [], label=label, color=color) for label, color in zip(state_labels, state_colors)]
fig.legend(handles=handles, bbox_to_anchor=(0.875, 0.5), loc='center left')
plt.subplots_adjust(right=0.875)
plt.savefig('out/line_betabinom-iter.png')
plt.close()

fig, axs = plt.subplots(len(state_labels), 1, figsize=(6.4, 1.6 * len(state_labels)))
for label, color in zip(state_labels, state_colors):
    for ax, param in zip(axs, state_labels):
        if param == label:
            continue
        xs = [record['iter_num'] for record in history]
        ys = [record['t_dists_norm'][label][param] for record in history]
        ax.plot(xs, ys, label=label, color=color)
        ax.set_ylabel(param)
axs[3].set_xlabel('Iteration')
handles = [Line2D([], [], label=label, color=color) for label, color in zip(state_labels, state_colors)]
fig.legend(handles=handles, bbox_to_anchor=(0.875, 0.5), loc='center left')
plt.subplots_adjust(left=0.15, right=0.875)
plt.savefig('out/line_jump-iter.png')
plt.close()

# Plot alignments
if not os.path.exists('out/traces/'):
    os.mkdir('out/traces/')

for OGid, labels in OGid2labels.items():
    # Load MSA
    msa = []
    for header, seq in read_fasta(f'../realign_hmmer/out/mafft/{OGid}.afa'):
        spid = re.search(spid_regex, header).group(1)
        msa.append({'spid': spid, 'seq': seq})
    msa = sorted(msa, key=lambda x: tip_order[x['spid']])  # Re-order sequences

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
    tree = tree_template.shear([record['spid'] for record in msa])
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
    for s, e_dist in model_json['e_dists'].items():
        a, b, pi, q0, q1, p0, p1 = [e_dist[param] for param in ['a', 'b', 'pi', 'q0', 'q1', 'p0', 'p1']]
        array1 = utils.get_betabinom_pmf(emit_seq, len(msa), a, b)
        array2 = utils.get_tree_pmf(tree, pi, q0, q1, p0, p1)
        e_dists_rv[s] = utils.ArrayRV(array1 * array2)
    model = homomorph.HMM(model_json['t_dists'], e_dists_rv, model_json['start_dist'])

    # Make plotting parameters
    kwargs_wide = {'figsize': (15, 6), 'height_ratio': 0.5, 'hspace': 0.2,
                   'data_max': 1.1, 'data_min': -0.1, 'data_labels': state_labels, 'data_colors': state_colors,
                   'msa_legend': True, 'legend_kwargs': {'bbox_to_anchor': (0.945, 0.5), 'loc': 'center left', 'fontsize': 8,
                                                         'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1}}
    adjust_wide = {'left': 0.015, 'bottom': 0.01, 'right': 0.94, 'top': 0.99}

    kwargs_tall = {'figsize': (8, 8), 'height_ratio': 0.5, 'hspace': 0.2,
                   'data_max': 1.1, 'data_min': -0.1, 'data_labels': state_labels, 'data_colors': state_colors,
                   'msa_legend': True, 'legend_kwargs': {'bbox_to_anchor': (0.90, 0.5), 'loc': 'center left', 'fontsize': 8,
                                                         'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1}}
    adjust_tall = {'left': 0.025, 'bottom': 0.01, 'right': 0.89, 'top': 0.99}

    # Plot labels
    lines = {label: np.zeros(len(msa[0]['seq'])) for label in state_labels}
    for start, stop, label in labels:
        lines[label][start:stop] = 1
    data = [lines[label] for label in state_labels]

    plot_msa_data([record['seq'] for record in msa], data, **kwargs_wide)
    plt.subplots_adjust(**adjust_wide)
    plt.savefig(f'out/traces/{OGid}_wide_labels.png')
    plt.close()

    plot_msa_data([record['seq'] for record in msa], data, **kwargs_tall)
    plt.subplots_adjust(**adjust_tall)
    plt.savefig(f'out/traces/{OGid}_tall_labels.png')
    plt.close()

    # Decode states and plot
    idx_seq = list(range(len(emit_seq)))  # Everything is pre-calculated, so emit_seq is the emit index
    fbs = model.forward_backward(idx_seq)
    data = [fbs[label] for label in state_labels]

    plot_msa_data([record['seq'] for record in msa], data, **kwargs_wide)
    plt.subplots_adjust(**adjust_wide)
    plt.savefig(f'out/traces/{OGid}_wide.png')
    plt.close()

    plot_msa_data([record['seq'] for record in msa], data, **kwargs_tall)
    plt.subplots_adjust(**adjust_tall)
    plt.savefig(f'out/traces/{OGid}_tall.png')
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