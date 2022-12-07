"""Plot decoded insertions in alignments with large gap contrasts."""

import json
import os
import re

import homomorph
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skbio
from src.draw import plot_msa_data
from src.ortho_MSA.trim import get_hull_slices, get_trim_slices
from src.ortho_MSA import utils
from src.utils import get_brownian_weights, read_fasta

ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'
state_labels = ['1A', '1B', '2', '3']
state_colors = ['C0', 'C3', 'C1', 'C2']

# Cutoffs for state 3 trims
posterior_high1 = 0.75
posterior_low1 = 0.01
profile_low1 = 0.1

# Cutoffs for state 2+3 trims
posterior_high2 = 0.9
posterior_low2 = 0.01

# Common cutoffs for gradients
gradient_high = 0.001
gradient_low = 0.001

alpha = 0.01
mean_min = 2
mean_trim = 5  # Number of "outliers" to remove before calculating mean

tree_template = skbio.read('../../ortho_tree/consensus_GTR2/out/NI.nwk', 'newick', skbio.TreeNode)
tree_order = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree_order.tips())}

df = pd.read_table('../gap_contrasts/out/total_sums.tsv')
df['norm1'] = df['total'] / df['gnidnum']
df['norm2'] = df['total'] / (df['gnidnum'] * df['len2'])

with open('../insertion_hmm/out/model.json') as file:
    model_json = json.load(file)

for label in ['norm1', 'norm2']:
    if not os.path.exists(f'out/{label}/'):
        os.makedirs(f'out/{label}/')

    head = df.sort_values(by=label, ascending=False).head(150)
    for rank, row in enumerate(head.itertuples()):
        # Load MSA
        input_msa = []
        for header, seq in read_fasta(f'../realign_fastas/out/{row.OGid}.afa'):
            ppid = re.search(ppid_regex, header).group(1)
            spid = re.search(spid_regex, header).group(1)
            input_msa.append({'header': header, 'ppid': ppid, 'spid': spid, 'seq': seq})
        input_msa = sorted(input_msa, key=lambda x: tip_order[x['spid']])

        # Create emission sequence
        col0 = []
        emit_seq = []
        for j in range(len(input_msa[0]['seq'])):
            col = [1 if input_msa[i]['seq'][j] in ['-', '.'] else 0 for i in range(len(input_msa))]
            emit0 = sum([c0 == c for c0, c in zip(col0, col)])
            emit_seq.append(emit0)  # The tree probabilities are pre-calculated, so emission value is its index
            col0 = col
        emit_seq = np.array(emit_seq)

        # Load tree and convert to vectors at tips
        tree = tree_template.shear([record['spid'] for record in input_msa])
        tips = {tip.name: tip for tip in tree.tips()}
        for record in input_msa:
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
            pmf1 = utils.get_betabinom_pmf(emit_seq, len(input_msa), a, b)
            pmf2 = utils.get_tree_pmf(tree, pi, q0, q1, p0, p1)
            e_dists_rv[s] = utils.ArrayRV(pmf1 * pmf2)
        model = homomorph.HMM(model_json['t_dists'], e_dists_rv, model_json['start_dist'])

        # Make plotting parameters
        kwargs = {'figsize': (15, 6),
                  'height_ratio': 0.5, 'hspace': 0.45, 'left': 0.015, 'right': 0.94, 'top': 0.99, 'bottom': 0.03,
                  'data_max': 1.05, 'data_min': -0.05,
                  'msa_legend': True,
                  'legend_kwargs': {'bbox_to_anchor': (0.945, 0.5), 'loc': 'center left', 'fontsize': 8,
                                    'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1}}

        # Decode states and plot
        idx_seq = list(range(len(input_msa[0]['seq'])))  # Everything is pre-calculated, so emit_seq is the emit index
        fbs = model.forward_backward(idx_seq)
        data = [fbs[label] for label in state_labels]

        plot_msa_data([record['seq'] for record in input_msa], data,
                      data_labels=state_labels, data_colors=state_colors, **kwargs)
        plt.savefig(f'out/{label}/{rank:03}_{row.OGid}.png')
        plt.close()

        # Calculate weights
        spid_set = set([record['spid'] for record in input_msa])
        spid2ppids = {spid: [] for spid in spid_set}
        for record in input_msa:
            ppid, spid = record['ppid'], record['spid']
            spid2ppids[spid].append(ppid)

        tree = tree_template.shear(spid_set)
        spid_weights = {tip.name: weight for tip, weight in get_brownian_weights(tree)}
        ppid_weights = {}
        for spid, ppids in spid2ppids.items():
            weight = spid_weights[spid] / len(ppids)  # Species weight is distributed evenly over all associated proteins
            for ppid in ppids:
                ppid_weights[ppid] = weight

        # Calculate profile
        weight_msa = np.zeros((len(input_msa), len(input_msa[0]['seq'])))
        for i, record in enumerate(input_msa):
            ppid, seq = record['ppid'], record['seq']
            weight = ppid_weights[ppid]
            for j, sym in enumerate(seq):
                if sym in ['-', '.']:
                    weight_msa[i, j] = weight
        profile = weight_msa.sum(axis=0)

        # Identify state 3 trims
        posterior = np.array(fbs['3'])
        posterior[profile <= profile_low1] = 0
        gradient = np.gradient(posterior)
        slices = get_trim_slices(profile, posterior, gradient, posterior_high1, posterior_low1, gradient_high, gradient_low)

        seq_slices = {record['ppid']: [] for record in input_msa}
        for s in slices:
            count_records = []
            for record in input_msa:
                ppid, seq = record['ppid'], record['seq']
                subseq = seq[s]
                count = len(subseq) - subseq.count('-') - subseq.count('.')
                count_records.append((ppid, count))
            count_records = sorted(count_records, key=lambda x: x[1])

            count_sum = sum([ppid_weights[ppid] * count for ppid, count in count_records[:-mean_trim]])
            weight_sum = sum([ppid_weights[ppid] for ppid, _ in count_records[:-mean_trim]])
            mean = max(count_sum / weight_sum, mean_min)
            p = 1 / (mean + 1)  # Geometric distribution on support 0, 1, ...
            k = np.log(alpha) / np.log(1 - p) - 1  # Expression for minimum k to achieve alpha significance

            for ppid, count in count_records:
                if count >= k:
                    seq_slices[ppid].append(s)

        # Identify state 2+3 trims
        posterior = np.array(fbs['3'])
        slices = get_hull_slices(posterior, gradient, posterior_high1, posterior_low1, gradient_low)
        for s in slices:
            posterior[s] = 0
        posterior += np.array(fbs['2'])
        gradient = np.gradient(posterior)
        slices = get_trim_slices(profile, posterior, gradient, posterior_high2, posterior_low2, gradient_high, gradient_low)

        # Create trimmed MSA with seq_slices
        trimmed_msa = []
        for record in input_msa:
            header, ppid = record['header'], record['ppid'],
            seq = list(record['seq'])
            for s in seq_slices[ppid]:
                seq[s] = (s.stop - s.start) * ['.']
            trimmed_msa.append({'header': header, 'seq': ''.join(seq)})

        trims = np.zeros(len(posterior))
        for s in slices:
            trims[s] = 1

        plot_msa_data([record['seq'] for record in trimmed_msa], [posterior, trims],
                      data_labels=['2+3', 'trim'], data_colors=['C1', 'C0'], **kwargs)
        plt.savefig(f'out/{label}/{rank:03}_{row.OGid}_trim.png')
        plt.close()

"""
DEPENDENCIES
../../ortho_tree/consensus_GTR2/consensus_GTR2.py
    ../../ortho_tree/consensus_GTR2/out/NI.nwk
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../gap_contrasts/gap_contrasts.py
    ../gap_contrasts/out/total_sums.tsv
../insertion_trim/decode.py
    ./insertion_trim/out/*.tsv
../realign_fastas/realign_fastas.py
    ../realign_fastas/out/*.afa
"""