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
from src.ortho_MSA.trim import get_slices
from src.ortho_MSA import utils
from src.utils import read_fasta

spid_regex = r'spid=([a-z]+)'
state_labels = ['1A', '1B', '2', '3']
state_colors = ['C0', 'C3', 'C1', 'C2']

posterior_high = 0.75
posterior_low = 0.5
gradient_high = 0.02
gradient_low = 0.001

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
    for i, row in enumerate(head.itertuples()):
        # Load MSA
        msa = []
        for header, seq in read_fasta(f'../realign_hmmer/out/mafft/{row.OGid}.afa'):
            spid = re.search(spid_regex, header).group(1)
            msa.append({'spid': spid, 'seq': seq})
        msa = sorted(msa, key=lambda x: tip_order[x['spid']])

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
            pmf1 = utils.get_betabinom_pmf(emit_seq, len(msa), a, b)
            pmf2 = utils.get_tree_pmf(tree, pi, q0, q1, p0, p1)
            e_dists_rv[s] = utils.ArrayRV(pmf1 * pmf2)
        model = homomorph.HMM(model_json['t_dists'], e_dists_rv, model_json['start_dist'])

        # Make plotting parameters
        kwargs = {'figsize': (15, 6), 'height_ratio': 0.5, 'hspace': 0.2,
                  'data_max': 1.1, 'data_min': -0.1, 'data_labels': state_labels, 'data_colors': state_colors,
                  'msa_legend': True,
                  'legend_kwargs': {'bbox_to_anchor': (0.945, 0.5), 'loc': 'center left', 'fontsize': 8,
                                    'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1}}
        adjust = {'left': 0.015, 'bottom': 0.01, 'right': 0.94, 'top': 0.99}

        # Decode states and plot
        idx_seq = list(range(len(msa[0]['seq'])))  # Everything is pre-calculated, so emit_seq is the emit index
        fbs = model.forward_backward(idx_seq)
        data = [fbs[label] for label in state_labels]

        plot_msa_data([record['seq'] for record in msa], data,
                      figsize=(15, 6), height_ratio=0.5, hspace=0.2,
                      data_max=1.1, data_min=-0.1, data_labels=state_labels, data_colors=state_colors,
                      msa_legend=True,
                      legend_kwargs={'bbox_to_anchor': (0.945, 0.5), 'loc': 'center left', 'fontsize': 8,
                                     'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1})
        plt.subplots_adjust(left=0.015, bottom=0.01, right=0.94, top=0.99)
        plt.savefig(f'out/{label}/{i:03}_{row.OGid}.png')
        plt.close()

        # Get trims and plot
        posterior = np.array(fbs['2']) + np.array(fbs['3'])
        gradient = np.gradient(posterior)
        slices = get_slices([(record['spid'], record['seq']) for record in msa],
                            posterior, gradient,
                            posterior_high, posterior_low,
                            gradient_high, gradient_low)
        trims = np.zeros(len(posterior))
        for s in slices:
            trims[s] = 1

        plot_msa_data([record['seq'] for record in msa], [posterior, trims],
                      figsize=(15, 6), height_ratio=0.5, hspace=0.2,
                      data_max=1.1, data_min=-0.1, data_labels=['2+3', 'trim'], data_colors=['C1', 'C0'],
                      msa_legend=True,
                      legend_kwargs={'bbox_to_anchor': (0.945, 0.5), 'loc': 'center left', 'fontsize': 8,
                                     'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1})
        plt.subplots_adjust(left=0.015, bottom=0.01, right=0.94, top=0.99)
        plt.savefig(f'out/{label}/{i:03}_{row.OGid}_trim.png')
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
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/mafft/*.afa
"""