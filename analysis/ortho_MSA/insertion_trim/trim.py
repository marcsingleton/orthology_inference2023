"""Trim state 2 and 3 regions to yield trimmed alignments."""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skbio
from src.ortho_MSA.trim import get_hull_slices, get_trim_slices
from src.utils import get_brownian_weights, read_fasta

ppid_regex = r'ppid=([A-Za-z0-9_.]+)'
spid_regex = r'spid=([a-z]+)'

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
tip_order = {tip.name: i for i, tip in enumerate(tree_template.tips())}

# Load OGids
OGids = []
with open('../realign_fastas/out/errors.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        OGid, error_flag1, error_flag2 = fields['OGid'], fields['error_flag1'], fields['error_flag2']
        if error_flag1 == 'False' and error_flag2 == 'False':
            OGids.append(OGid)

if not os.path.exists('out/trims/'):
    os.mkdir('out/trims/')

rows1, rows2 = [], []
for OGid in OGids:
    # Load MSA
    input_msa = []
    for header, seq in read_fasta(f'../realign_fastas/out/{OGid}.afa'):
        ppid = re.search(ppid_regex, header).group(1)
        spid = re.search(spid_regex, header).group(1)
        input_msa.append({'header': header, 'ppid': ppid, 'spid': spid, 'seq': seq})

    df = pd.read_table(f'out/posteriors/{OGid}.tsv')

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
    posterior = df['3'].to_numpy()
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
                rows1.append({'OGid': OGid, 'ppid': ppid, 'start': s.start, 'stop': s.stop, 'count': count})

    # Identify state 2+3 trims
    posterior = df['3'].to_numpy(copy=True)
    slices = get_hull_slices(posterior, gradient, posterior_high1, posterior_low1, gradient_low)
    for s in slices:
        posterior[s] = 0
    posterior += df['2'].to_numpy()
    gradient = np.gradient(posterior)
    slices = get_trim_slices(profile, posterior, gradient, posterior_high2, posterior_low2, gradient_high, gradient_low)
    for s in slices:
        rows2.append({'OGid': OGid, 'colnum': len(input_msa[0][1]), 'start': s.start, 'stop': s.stop,
                      'posterior2': df.loc[s, '2'].sum(), 'posterior3': df.loc[s, '3'].sum()})

    # Invert slices to get untrimmed regions
    region_slices, stop = [], 0
    for s in slices:
        region_slices.append(slice(stop, s.start))
        stop = s.stop
    region_slices.append(slice(stop, len(input_msa[0]['seq'])))

    # Create trimmed MSA from two sets of slices
    trimmed_msa = []
    for record in input_msa:
        header, ppid = record['header'], record['ppid'],
        seq1 = list(record['seq'])
        for s in seq_slices[ppid]:
            seq1[s] = (s.stop - s.start) * ['.']
        seq2 = []
        for s in region_slices:
            seq2.extend(seq1[s])
        trimmed_msa.append({'header': header, 'seq': ''.join(seq2)})

    # Remove excess gaps
    slices, idx = [], None
    for j in range(len(trimmed_msa[0]['seq'])):
        for i in range(len(trimmed_msa)):
            sym = trimmed_msa[i]['seq'][j]
            if sym not in ['-', '.']:
                if idx is None:  # Store position only if new slice is not started
                    idx = j
                break
        else:
            if idx is not None:
                slices.append(slice(idx, j))
                idx = None
    if idx is not None:  # Add final slice to end
        slices.append(slice(idx, len(trimmed_msa[0]['seq'])))

    with open(f'out/trims/{OGid}.afa', 'w') as file:
        for record in trimmed_msa:
            header, seq1 = record['header'], record['seq']
            seq2 = ''.join([seq1[s] for s in slices])  # Remove gap only columns
            seqstring = '\n'.join([seq2[i:i+80] for i in range(0, len(seq2), 80)])
            file.write(f'{header}\n{seqstring}\n')

# Plot stats
df = pd.DataFrame(rows1)
df.to_csv('out/trim_stats.tsv', sep='\t', index=False)

df['length'] = df['stop'] - df['start']
df['length_ratio'] = df['length'] / df['colnum']
df['norm2'] = df['posterior2'] / df['length']
df['norm3'] = df['posterior3'] / df['length']
groups = df.groupby('OGid')

# Pie chart by presence of trims
values = [len(set(OGids)) - df['OGid'].nunique(), df['OGid'].nunique()]
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
../realign_fastas/realign_fastas.py
    ../realign_fastas/out/errors.tsv
    ../realign_fastas/out/*.afa
./decode.py
    ./out/posteriors/*.tsv
"""