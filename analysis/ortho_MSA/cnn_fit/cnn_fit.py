"""Fit CNN to identify mis-aligned regions."""

import os
import random
import re
from math import ceil, floor

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tensorflow as tf
from matplotlib.gridspec import GridSpec
from src.draw import plot_msa_data, default_sym2color
from src.utils import read_fasta


class BatchGenerator(tf.keras.utils.Sequence):
    """Label, batch, and pad protein sequence data.

    Only complete batches are returned, so a single epoch may not train on every example.
    """
    def __init__(self, records, batch_size, shuffle=True):
        self.records = records
        self.indices = np.arange(len(self.records))
        self.batch_size = batch_size
        self.shuffle = shuffle
        self.on_epoch_end()

    def __len__(self):
        """Return number of batches."""
        return floor(len(self.records) / self.batch_size)

    def __getitem__(self, index):
        """Generate one batch of data."""
        indices = self.indices[index*self.batch_size:(index+1)*self.batch_size]
        records = [self.records[i] for i in indices]
        length = max([len(record[2]) for record in records])

        x1 = np.zeros((self.batch_size, length, len(alphabet)))  # Input1 array
        x2 = np.zeros((self.batch_size, length, len(alphabet)))  # Input2 array
        y = np.zeros((self.batch_size, length))  # Output array
        w = np.zeros((self.batch_size, length))  # Weight array
        for i, record in enumerate(records):
            _, _, profile, seq, labels, weights = record
            x1[i, :len(profile)] = profile
            x2[i, :len(seq)] = seq
            y[i, :len(labels)] = labels
            w[i, :len(weights)] = weights
        y = np.expand_dims(y, -1)  # Add additional axis so output shape matches

        return [x1, x2], y, w

    def on_epoch_end(self):
        """Shuffles data after each epoch."""
        if self.shuffle:
            np.random.shuffle(self.indices)


alphabet = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X', '-']
sym2idx = {sym: i for i, sym in enumerate(alphabet)}
ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
epochs = 300
batch_size = 10
validation_split = 0.2
embedding_dim = 2
regularizer = tf.keras.regularizers.L2(0.0025)
tf.keras.utils.set_random_seed(930715)  # Make validation split and all TensorFlow operations consistent

# Load labels
OGid2ppids = {}
ppid2labels = {}
with open('labels.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        if line.startswith('#'):
            continue
        OGid, ppid, start, stop, label, active = line.rstrip('\n').split('\t')
        if active == 'False':
            continue
        try:
            OGid2ppids[OGid].add(ppid)
        except KeyError:
            OGid2ppids[OGid] = {ppid}
        record = (int(start), int(stop), int(label))
        try:
            ppid2labels[ppid].append(record)
        except KeyError:
            ppid2labels[ppid] = [record]

# Create records
records = []
for OGid, ppids in OGid2ppids.items():
    try:
        msa1 = read_fasta(f'../align_fastas1/out/{OGid}.afa')
    except FileNotFoundError:
        msa1 = read_fasta(f'../align_fastas2/out/{OGid}.afa')

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

# Split data
random.shuffle(records)
split_idx = ceil(len(records) * (1 - validation_split))
train_records = records[:split_idx]
validation_records = records[split_idx:]

# Batch data
train_batches = BatchGenerator(train_records, batch_size)
validation_batches = BatchGenerator(validation_records, batch_size)

# Build model
embedding = tf.keras.layers.Dense(embedding_dim, activation='linear', use_bias=False, name='embedding1')
conv1_1 = tf.keras.layers.Conv1D(2, 60, dilation_rate=3, padding='same', kernel_regularizer=regularizer, activation='tanh', name='conv1-1')
conv1_2 = tf.keras.layers.Conv1D(2, 10, dilation_rate=1, padding='same', kernel_regularizer=regularizer, activation='tanh', name='conv1-2')
conv2_1 = tf.keras.layers.Conv1D(2, 10, dilation_rate=1, kernel_regularizer=regularizer, padding='same', activation='tanh', name='conv2-1')

input1 = tf.keras.layers.Input(shape=(None, len(alphabet)), name='input1')
input2 = tf.keras.layers.Input(shape=(None, len(alphabet)), name='input2')
x = tf.keras.layers.concatenate([embedding(input1), embedding(input2)], name='concatenate1')
x = tf.keras.layers.concatenate([conv1_1(x), conv1_2(x)], name='concatenate2')
x = conv2_1(x)
outputs = tf.keras.layers.Dense(1, activation='sigmoid', name='dense1')(x)

model = tf.keras.Model(inputs=[input1, input2], outputs=outputs)
model.compile(loss='binary_crossentropy', optimizer='adam',
              metrics=[tf.keras.metrics.BinaryAccuracy(), tf.keras.metrics.Recall(), tf.keras.metrics.Precision()])
model.summary()

# Train model
history = model.fit(train_batches, epochs=epochs, validation_data=validation_batches)

# Save model and history
if not os.path.exists('out/'):
    os.mkdir('out/')

df = pd.DataFrame(history.history)
df.to_csv('out/history.tsv', sep='\t', index=False)
model.save('out/model.h5')

# Plot embedding
fig, ax = plt.subplots()
weights = embedding.get_weights()[0]
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
positions = {'conv1-1-0': gs[0, :],
             'conv1-1-1': gs[1, :],
             'conv1-2-0': gs[2, 0], 'conv2-1-0': gs[2, 1],
             'conv1-2-1': gs[3, 0], 'conv2-1-1': gs[3, 1]}
for layer in [conv1_1, conv1_2, conv2_1]:
    weights = layer.get_weights()[0]
    for i in range(weights.shape[2]):
        ax = fig.add_subplot(positions[f'{layer.name}-{i}'])
        ax.imshow(weights[..., i].transpose())
        ax.set_title(f'{layer.name}-{i}')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
fig.savefig('out/weights.png')
plt.close()

# Plot curves
plt.plot(df['loss'], label='train')
plt.plot(df['val_loss'], label='validation')
plt.xlabel('Epoch')
plt.ylabel('Loss')
plt.legend()
plt.savefig('out/line_loss-epoch.png')
plt.close()

plt.plot(df['binary_accuracy'], label='accuracy')
plt.plot(df['recall'], label='recall')
plt.plot(df['precision'], label='precision')
plt.plot(df['val_binary_accuracy'], label='validation accuracy')
plt.plot(df['val_recall'], label='validation recall')
plt.plot(df['val_precision'], label='validation precision')
plt.xlabel('Epoch')
plt.ylabel('Value')
plt.legend()
plt.savefig('out/line_metrics-epoch.png')
plt.close()

# Show decoding curves with MSA
for label, records in [('train', train_records), ('validation', validation_records)]:
    if not os.path.exists(f'out/{label}/'):
        os.mkdir(f'out/{label}/')

    for record in records:
        OGid, ppid, profile, seq, labels, weights = record
        output = tf.squeeze(model([np.expand_dims(profile, 0), np.expand_dims(seq, 0)]))  # Expand and contract dims

        try:
            msa = read_fasta(f'../align_fastas1/out/{OGid}.afa')
        except FileNotFoundError:
            msa = read_fasta(f'../align_fastas2/out/{OGid}.afa')
        msa = [(re.search(ppid_regex, header).group(1), seq) for header, seq in msa]

        data = [output, labels, weights]
        msa_labels = [header if header == ppid else '' for header, _ in msa]
        plot_msa_data([seq for _, seq in msa], data,
                      msa_labels=msa_labels, msa_length=1, msa_width=0.25, msa_pad=1.1, msa_labelsize=2.5,
                      height_ratio=0.5, hspace=0.2, data_max=1.1, data_min=-0.1,
                      legend=True, legend_kwargs={'bbox_to_anchor': (0.925, 0.5), 'loc': 'center left', 'fontsize': 8, 'handletextpad': 0, 'markerscale': 1.25})
        plt.subplots_adjust(left=0.05, bottom=0.01, right=0.925, top=0.99)
        plt.savefig(f'out/{label}/{OGid}_{ppid}.png', dpi=100)
        plt.close()

"""
DEPENDENCIES
../align_fastas1/align_fastas1.py
    ../align_fastas1/out/*.afa
../align_fastas2/align_fastas2.py
    ../align_fastas2/out/*.afa
./labels.tsv
"""