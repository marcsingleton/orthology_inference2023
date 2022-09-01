"""Plot example alignments segmented via posterior decoding."""

import os
import json
import re

import homomorph
import matplotlib.pyplot as plt
import numpy as np
import skbio
from matplotlib.lines import Line2D
from src.ortho_MSA import utils
from src.draw import plot_msa_data
from src.utils import read_fasta

ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'
state_labels = ['1', '2']
state_colors = ['C0', 'C1']

tree_template = skbio.read('../../ortho_tree/consensus_GTR2/out/NI.nwk', 'newick', skbio.TreeNode)
tree_order = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree_order.tips())}

# Load labels
ids2labels = {}
label_set = set()
with open('labels.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        OGid, ppid, start, stop, label = fields['OGid'], fields['ppid'], int(fields['start']), int(fields['stop']), fields['label']
        label_set.add(label)
        try:
            ids2labels[(OGid, ppid)].append((start, stop, label))
        except KeyError:
            ids2labels[(OGid, ppid)] = [(start, stop, label)]

if set(state_labels) != label_set:
    raise RuntimeError('state_labels is not equal to set of state_labels')

# Load history and model parameters
with open('out/history.json') as file:
    history = json.load(file)
with open('out/model.json') as file:
    model_json = json.load(file)

# Plot state distribution
counts = {label: 0 for label in state_labels}
for labels in ids2labels.values():
    for start, stop, label in labels:
        counts[label] += stop - start
values = [counts[label] for label in state_labels]
labels = [f'{label}\n{value:,}' for label, value in zip(state_labels, values)]
plt.pie(values, colors=state_colors, labels=labels, labeldistance=1.15, textprops={'ha': 'center'})
plt.title('Distribution of position labels across states')
plt.savefig('out/pie_labels.png')
plt.close()

# Write some metrics to file
output = f"""\
Number of sequences: {len(ids2labels)}
Number of labelled positions: {sum(counts.values()):,}
"""
with open('out/output.txt', 'w') as file:
    file.write(output)

# Plot loss curve
xs = [record['iter_num'] for record in history]
ys = [record['ll'] for record in history]
plt.plot(xs, ys)
plt.xlabel('Iteration')
plt.ylabel('Conditional log-likelihood')
plt.savefig('out/line_ll-iter.png')
plt.close()

# Plot model parameters
params = ['pi', 'q0', 'q1']
fig, axs = plt.subplots(len(params), 1)
label, color = state_labels[0], state_colors[0]
for ax, param in zip(axs, params):
    xs = [record['iter_num'] for record in history]
    ys = [record['e_dists_norm'][label][param] for record in history]
    ax.plot(xs, ys, label=label, color=color)
    ax.set_ylabel(param)
axs[-1].set_xlabel('Iteration')
handles = [Line2D([], [], label=label, color=color) for label, color in zip(state_labels[:1], state_colors[:1])]
fig.legend(handles=handles, bbox_to_anchor=(0.875, 0.5), loc='center left')
plt.subplots_adjust(right=0.875)
plt.savefig('out/line_rate-iter.png')
plt.close()

param = 'p'
fig, ax = plt.subplots()
label, color = state_labels[1], state_colors[1]
xs = [record['iter_num'] for record in history]
ys = [record['e_dists_norm'][label][param] for record in history]
ax.plot(xs, ys, label=label, color=color)
ax.set_ylabel(param)
ax.set_xlabel('Iteration')
handles = [Line2D([], [], label=label, color=color) for label, color in zip(state_labels[1:], state_colors[1:])]
fig.legend(handles=handles, bbox_to_anchor=(0.875, 0.5), loc='center left')
plt.subplots_adjust(right=0.875)
plt.savefig('out/line_p-iter.png')
plt.close()

fig, axs = plt.subplots(len(state_labels), 1)
for label, color in zip(state_labels, state_colors):
    for ax, param in zip(axs, state_labels):
        if param == label:
            continue
        xs = [record['iter_num'] for record in history]
        ys = [record['t_dists_norm'][label][param] for record in history]
        ax.plot(xs, ys, label=label, color=color)
        ax.set_ylabel(param)
for ax in axs:
    ax.set_yscale('log')
axs[-1].set_xlabel('Iteration')
handles = [Line2D([], [], label=label, color=color) for label, color in zip(state_labels, state_colors)]
fig.legend(handles=handles, bbox_to_anchor=(0.875, 0.5), loc='center left')
plt.subplots_adjust(left=0.15, right=0.875)
plt.savefig('out/line_transition-iter.png')
plt.close()

# Plot alignments
if not os.path.exists('out/traces/'):
    os.mkdir('out/traces/')

for (OGid, ppid), labels in ids2labels.items():
    # Load MSA
    msa, ppid2spid = [], {}
    for header, seq in read_fasta(f'../realign_hmmer/out/mafft/{OGid}.afa'):
        msa_ppid = re.search(ppid_regex, header).group(1)
        msa_spid = re.search(spid_regex, header).group(1)
        msa.append({'ppid': msa_ppid, 'spid': msa_spid, 'seq': seq})
        ppid2spid[msa_ppid] = msa_spid
    msa = sorted(msa, key=lambda x: tip_order[x['spid']])

    # Load tree and convert to vectors at tips
    tree = tree_template.shear([record['spid'] for record in msa])
    tips = {tip.name: tip for tip in tree.tips()}
    tree.tip_dict = tips
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
        if s == '1':
            pi, q0, q1 = [e_dist[param] for param in ['pi', 'q0', 'q1']]
            pmf = utils.get_tip_pmf(tree, ppid2spid[ppid], pi, q0, q1, 0, 0)
            e_dists_rv[s] = utils.ArrayRV(pmf)
        elif s == '2':
            p = e_dist['p']
            conditional = tree.tip_dict[ppid2spid[ppid]].conditional[1]  # Second row is gaps=0, non-gaps=1
            pmf = utils.get_bernoulli_pmf(conditional, p)
            e_dists_rv[s] = utils.ArrayRV(pmf)
    model = homomorph.HMM(model_json['t_dists'], e_dists_rv, model_json['start_dist'])

    # Make plotting parameters
    msa_labels = [record['ppid'] if record['ppid'] == ppid else '' for record in msa]

    kwargs_wide = {'figsize': (15, 6), 'height_ratio': 0.5, 'hspace': 0.2,
                   'msa_labels': msa_labels, 'msa_labelsize': 4,
                   'data_max': 1.1, 'data_min': -0.1, 'data_labels': state_labels, 'data_colors': state_colors,
                   'msa_legend': True, 'legend_kwargs': {'bbox_to_anchor': (0.945, 0.5), 'loc': 'center left', 'fontsize': 8,
                                                         'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1}}
    adjust_wide = {'left': 0.04, 'bottom': 0.01, 'right': 0.94, 'top': 0.99}

    kwargs_tall = {'figsize': (8, 8), 'height_ratio': 0.5, 'hspace': 0.2,
                   'msa_labels': msa_labels, 'msa_labelsize': 4,
                   'data_max': 1.1, 'data_min': -0.1, 'data_labels': state_labels, 'data_colors': state_colors,
                   'msa_legend': True, 'legend_kwargs': {'bbox_to_anchor': (0.90, 0.5), 'loc': 'center left', 'fontsize': 8,
                                                         'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1}}
    adjust_tall = {'left': 0.06, 'bottom': 0.01, 'right': 0.89, 'top': 0.99}

    # Plot labels
    lines = {label: np.zeros(len(msa[0]['seq'])) for label in state_labels}
    for start, stop, label in labels:
        lines[label][start:stop] = 1
    data = [lines[label] for label in state_labels]

    plot_msa_data([record['seq'] for record in msa], data, **kwargs_wide)
    plt.subplots_adjust(**adjust_wide)
    plt.savefig(f'out/traces/{OGid}-{ppid}_wide_labels.png')
    plt.close()

    plot_msa_data([record['seq'] for record in msa], data, **kwargs_tall)
    plt.subplots_adjust(**adjust_tall)
    plt.savefig(f'out/traces/{OGid}-{ppid}_tall_labels.png')
    plt.close()

    # Decode states and plot
    idx_seq = list(range(len(msa[0]['seq'])))  # Everything is pre-calculated, so emit_seq is the emit index
    fbs = model.forward_backward(idx_seq)
    data = [fbs[label] for label in state_labels]

    plot_msa_data([record['seq'] for record in msa], data, **kwargs_wide)
    plt.subplots_adjust(**adjust_wide)
    plt.savefig(f'out/traces/{OGid}-{ppid}_wide.png')
    plt.close()

    plot_msa_data([record['seq'] for record in msa], data, **kwargs_tall)
    plt.subplots_adjust(**adjust_tall)
    plt.savefig(f'out/traces/{OGid}-{ppid}_tall.png')
    plt.close()

"""
DEPENDENCIES
../../ortho_tree/consensus_GTR2/consensus_GTR2.py
    ../../ortho_tree/consensus_GTR2/out/NI.nwk
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../insertion_trim/trim.py
    ../trim_extract/out/*.afa
./fit.py
    ./out/model.json
./labels.tsv
"""