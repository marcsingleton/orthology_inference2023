"""Fit CNN to identify mis-aligned regions."""

import os
import random
import re
from math import floor

import numpy as np
import pandas as pd
import tensorflow as tf
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
spid_regex = r'spid=([a-z]+)'

epochs = 300
batch_size = 30
validation_split = 0.1
embedding_dim = 2
regularizer = tf.keras.regularizers.L2(0.0075)
tf.keras.utils.set_random_seed(930715)  # Make validation split and all TensorFlow operations consistent

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

# Split data
random.shuffle(records)

# Batch data
batches = BatchGenerator(records, batch_size)

# Build model
embedding = tf.keras.layers.Dense(embedding_dim, activation='linear', use_bias=False, name='embedding1')
conv1_1 = tf.keras.layers.Conv1D(2, 60, dilation_rate=3, padding='same', kernel_regularizer=regularizer, activation='tanh', name='conv1_1')
conv1_2 = tf.keras.layers.Conv1D(2, 10, dilation_rate=1, padding='same', kernel_regularizer=regularizer, activation='tanh', name='conv1_2')
conv2_1 = tf.keras.layers.Conv1D(2, 10, dilation_rate=1, kernel_regularizer=regularizer, padding='same', activation='tanh', name='conv2_1')

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
history = model.fit(batches, epochs=epochs)

# Save model and history
if not os.path.exists('out/'):
    os.mkdir('out/')

df = pd.DataFrame(history.history)
df.to_csv('out/history.tsv', sep='\t', index=False)
model.save('out/model.h5')

"""
DEPENDENCIES
../get_repseqs/get_repseqs.py
    ../get_repseqs/out/*.afa
./labels.tsv
"""