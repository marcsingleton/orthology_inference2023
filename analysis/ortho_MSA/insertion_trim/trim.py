"""Trim state 1 segments to yield trimmed alignments."""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from src.ortho_MSA.trim import get_slices
from src.utils import read_fasta

posterior_high = 0.9
posterior_low = 0.05
gradient_high = 0.02
gradient_low = 0.001

# Load OGids
OGids = []
with open('../realign_hmmer/out/errors.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        OGid, error_flag = fields['OGid'], fields['error_flag']
        if error_flag == 'False':
            OGids.append(OGid)

if not os.path.exists('out/trims/'):
    os.mkdir('out/trims/')

rows = []
for OGid in OGids:
    msa = read_fasta(f'../realign_hmmer/out/mafft/{OGid}.afa')

    # Load decoded states and calculate derivative
    df = pd.read_table(f'out/posteriors/{OGid}.tsv')
    posterior = (df['2'] + df['3']).to_numpy()
    gradient = np.gradient(posterior)

    # Find trimmed regions
    slices = get_slices(msa, posterior, gradient, posterior_high, posterior_low, gradient_high, gradient_low)

    # Invert slices
    inverts, stop = [], 0
    for s in slices:
        inverts.append(slice(stop, s.start))
        stop = s.stop
    inverts.append(slice(stop, len(msa[0][1])))

    # Write trimmed MSA
    with open(f'out/trims/{OGid}.afa', 'w') as file:
        for header, seq1 in msa:
            seq2 = ''.join([seq1[s] for s in inverts])
            seqstring = '\n'.join([seq2[i:i+80] for i in range(0, len(seq2), 80)])
            file.write(f'{header}\n{seqstring}\n')

    # Store data about trims and OGs
    for s in slices:
        rows.append({'OGid': OGid, 'colnum': len(msa[0][1]), 'start': s.start, 'stop': s.stop,
                     'posterior2': df.loc[s, '2'].sum(), 'posterior3': df.loc[s, '3'].sum()})

# Plot stats
df = pd.DataFrame(rows)
df.to_csv('out/trim_stats.tsv', sep='\t', index=False)

df['length'] = df['stop'] - df['start']
df['length_ratio'] = df['length'] / df['colnum']
df['norm2'] = df['posterior2'] / df['length']
df['norm3'] = df['posterior3'] / df['length']
groups = df.groupby('OGid')

# Pie chart by presence of trims
values = [len(set(OGids) - set(df['OGid'])), len(set(df['OGid']))]
labels = [f'{label}\n{value:,}' for label, value in zip(['w/o trims', 'w/ trims'], values)]
fig, ax = plt.subplots()
ax.pie(values, labels=labels, labeldistance=1.15)
ax.set_title('OGs w/ and w/o trims')
fig.savefig('out/pie_trims.png')
plt.close()

# Distribution of number of trims
counts = groups.size().value_counts()
fig, ax = plt.subplots()
ax.bar(counts.index, counts.values, width=1)
ax.set_xlabel('Number of trims in OG')
ax.set_ylabel('Number of OGs')
fig.savefig('out/hist_OGnum-trimnum.png')
plt.close()

# Distribution of length of trims
fig, ax = plt.subplots()
ax.hist(df['length'], bins=100)
ax.set_xlabel('Length of trim')
ax.set_ylabel('Number of trims')
fig.savefig('out/hist_trimnum-length.png')
ax.set_yscale('log')
fig.savefig('out/hist_trimnum-length_log.png')
plt.close()

# Distribution of length ratio of trims
fig, ax = plt.subplots()
ax.hist(df['length_ratio'], bins=50)
ax.set_xlabel('Length ratio of trim')
ax.set_ylabel('Number of trims')
fig.savefig('out/hist_trimnum-ratio.png')
plt.close()

# Distribution of length of total trims
fig, ax = plt.subplots()
ax.hist(groups['length'].sum(), bins=100)
ax.set_xlabel('Length of trims in OG')
ax.set_ylabel('Number of OGs')
fig.savefig('out/hist_OGnum-length.png')
ax.set_yscale('log')
fig.savefig('out/hist_OGnum-length_log.png')
plt.close()

# Distribution of length ratio of total trims
fig, ax = plt.subplots()
ax.hist(groups['length_ratio'].sum(), bins=100)
ax.set_xlabel('Total length ratio of trims in OG')
ax.set_ylabel('Number of OGs')
fig.savefig('out/hist_OGnum-ratio.png')
plt.close()

# Hexbin of length ratios vs number of trims
fig = plt.figure(figsize=(6, 6), layout='constrained')
gs = fig.add_gridspec(4, 2, height_ratios=(1, 2, 0.15, 0.1), width_ratios=(4, 1))
ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

hb = ax.hexbin(groups.size(), groups['length_ratio'].sum(), bins='log', gridsize=50, mincnt=1, linewidth=0)
cax = fig.add_subplot(gs[3, 0])
fig.colorbar(hb, cax=cax, orientation='horizontal')

counts = groups.size().value_counts()
ax_histx.bar(counts.index, counts.values, width=1)
ax_histy.hist(groups['length_ratio'].sum(), bins=100, orientation='horizontal')

ax.set_xlabel('Number of trims in OG')
ax.set_ylabel('Total length ratio of trims in OG')
fig.savefig('out/hexbin_ratio-trimnum.png')
plt.close()

# Hexbin of posterior2 vs posterior3
fig, ax = plt.subplots()
hb = ax.hexbin(df['norm2'], df['norm3'], bins='log', gridsize=25, mincnt=1, linewidth=0)
ax.set_xlabel('Average state 2 posterior in trim')
ax.set_ylabel('Average state 3 posterior in trim')
fig.colorbar(hb)
fig.savefig('out/hexbin_norm3-norm2.png')
plt.close()

"""
DEPENDENCIES
../realign_hmmer/realign_hmmer.py
    ../realign_hmmer/out/errors.tsv
    ../realign_hmmer/out/mafft/*.afa
./decode.py
    ./out/posteriors/*.tsv
"""