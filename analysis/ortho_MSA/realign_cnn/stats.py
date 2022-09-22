"""Plot metrics associated with fitting CNN."""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skbio
import tensorflow as tf
from matplotlib.gridspec import GridSpec
from src.draw import plot_msa_data, default_sym2color
from src.utils import read_fasta

alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X', '-']
sym2idx = {sym: i for i, sym in enumerate(alphabet)}
ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'

tree = skbio.read('../../ortho_tree/consensus_LG/out/100R_NI.nwk', 'newick', skbio.TreeNode)
tip_order = {tip.name: i for i, tip in enumerate(tree.tips())}

# Load labels
OGid2ppids = {}
ppid2labels = {}
with open('labels.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        if line.startswith('#'):
            continue
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        if fields['active'] == 'False':
            continue
        OGid, ppid = fields['OGid'], fields['ppid']
        try:
            OGid2ppids[OGid].add(ppid)
        except KeyError:
            OGid2ppids[OGid] = {ppid}
        start, stop, label = int(fields['start']), int(fields['stop']), int(fields['label'])
        try:
            ppid2labels[ppid].append((start, stop, label))
        except KeyError:
            ppid2labels[ppid] = [(start, stop, label)]

# Create records
records = []
for OGid, ppids in OGid2ppids.items():
    msa1 = read_fasta(f'../get_repseqs/out/{OGid}.afa')

    # Convert alignment to indices and make profile
    ppid2idx, msa2 = {}, []
    for i, (header, seq1) in enumerate(msa1):
        ppid = re.search(ppid_regex, header).group(1)
        ppid2idx[ppid] = i
        seq2 = [sym2idx.get(sym, sym2idx['X']) for sym in seq1]  # All non-standard symbols mapped to X
        msa2.append(seq2)
    msa2 = tf.keras.utils.to_categorical(msa2, len(alphabet))
    profile = msa2.sum(axis=0) / len(msa2)

    # Create labels
    for ppid in ppids:
        labels = np.zeros(msa2.shape[1])
        weights = np.zeros(msa2.shape[1])
        for start, stop, label in ppid2labels[ppid]:
            labels[start:stop] = label
            weights[start:stop] = 1
        records.append((OGid, ppid, profile, msa2[ppid2idx[ppid]], labels, weights))

# Load model and history
df = pd.read_table('out/history.tsv')
model = tf.keras.models.load_model('out/model.h5')
layers = {layer.name: layer for layer in model.layers}

# Count residue labels
positive, negative = 0, 0
for record in records:
    labels, weights = record[4], record[5]
    positive += labels.sum()
    negative += weights.sum() - labels.sum()

# Write some metrics to file
output = f"""\
Number of train examples: {len(records)}
Number of positive train residues: {positive}
Number of negative train residues: {negative}
"""
with open('out/output.txt', 'w') as file:
    file.write(output)


values = [negative, positive]
labels = [f'{label}\n{value:,g}' for label, value in zip(['negative', 'positive'], values)]
plt.pie(values, labels=labels, labeldistance=1.25, textprops={'ha': 'center'})
plt.title('Distribution of position labels across states')
plt.savefig('out/bar_residues-data.png')
plt.close()

# Plot embedding
fig, ax = plt.subplots()
weights = layers['embedding1'].get_weights()[0]
ax.scatter(weights[:, 0], weights[:, 1], s=75,
           c=[f'#{default_sym2color[sym]}' for sym in alphabet],
           edgecolors=['black' if sym == '-' else 'none' for sym in alphabet])
for sym, weight in zip(alphabet, weights):
    ax.annotate(sym, xy=weight, va='center', ha='center', size=5, fontweight='bold', fontfamily='monospace')
ax.set_xlabel('Embedding axis 1')
ax.set_ylabel('Embedding axis 2')
plt.savefig('out/scatter_embedding.png')
plt.close()

# Plot weights
fig = plt.figure(layout='tight')
gs = GridSpec(4, 2, figure=fig)
positions = {'conv1_1_0': gs[0, :],
             'conv1_1_1': gs[1, :],
             'conv1_2_0': gs[2, 0], 'conv2_1_0': gs[2, 1],
             'conv1_2_1': gs[3, 0], 'conv2_1_1': gs[3, 1]}
for layer_name in ['conv1_1', 'conv1_2', 'conv2_1']:
    layer = layers[layer_name]
    weights = layer.get_weights()[0]
    for i in range(weights.shape[2]):
        ax = fig.add_subplot(positions[f'{layer.name}_{i}'])
        ax.imshow(weights[..., i].transpose())
        ax.set_title(f'{layer.name}_{i}')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
fig.savefig('out/weights.png')
plt.close()

# Plot curves
plt.plot(df['loss'], label='train')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.savefig('out/line_loss-epoch.png')
plt.close()

plt.plot(df['binary_accuracy'], label='accuracy')
plt.plot(df['recall'], label='recall')
plt.plot(df['precision'], label='precision')
plt.xlabel('Epoch')
plt.ylabel('Value')
plt.legend()
plt.savefig('out/line_metrics-epoch.png')
plt.close()

# Show decoding curves with MSA
if not os.path.exists(f'out/traces/'):
    os.mkdir(f'out/traces/')

for record in records:
    OGid, ppid, profile, seq, labels, weights = record
    output = tf.squeeze(model([np.expand_dims(profile, 0), np.expand_dims(seq, 0)]))  # Expand and contract dims

    msa = []
    for header, seq in read_fasta(f'../get_repseqs/out/{OGid}.afa'):
        msa_ppid = re.search(ppid_regex, header).group(1)
        msa_spid = re.search(spid_regex, header).group(1)
        msa.append({'ppid': msa_ppid, 'spid': msa_spid, 'seq': seq})
    msa = sorted(msa, key=lambda x: tip_order[x['spid']])

    data = [output, labels, weights]
    msa_labels = [msa_record['ppid'] if msa_record['ppid'] == ppid else '' for msa_record in msa]
    plot_msa_data([msa_record['seq'] for msa_record in msa], data, figsize=(15, 6),
                  msa_labels=msa_labels, msa_ticklength=1, msa_tickwidth=0.25, msa_tickpad=1.1, msa_labelsize=5,
                  height_ratio=0.5, hspace=0.2, data_max=1.1, data_min=-0.1, data_labels=['output', 'label', 'weight'],
                  msa_legend=True, legend_kwargs={'bbox_to_anchor': (0.945, 0.5), 'loc': 'center left', 'fontsize': 8, 'handletextpad': 0.5, 'markerscale': 1.25, 'handlelength': 1})
    plt.subplots_adjust(left=0.04, bottom=0.01, right=0.94, top=0.99)
    plt.savefig(f'out/traces/{OGid}_{ppid}.png', dpi=500)
    plt.close()

"""
DEPENDENCIES
../../ortho_tree/consensus_LG/consensus_LG.py
    ../../ortho_tree/consensus_LG/out/100R_NI.nwk
../get_repseqs/get_repseqs.py
    ../get_repseqs/out/*.afa
./labels.tsv
./realign_cnn.py
    ./out/history.tsv
    ./out/model.h5
"""