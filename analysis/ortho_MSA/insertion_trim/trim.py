"""Trim state 2 and 3 regions to yield trimmed alignments."""

import os
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import skbio
from src.ortho_MSA.trim import get_complement_slices, get_hull_slices, get_trim_slices
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
        subseq_records = []
        for i, record in enumerate(input_msa):
            ppid, seq = record['ppid'], record['seq']
            subseq = weight_msa[i, s] == 0
            sym_count = subseq.sum()
            align_count = ((weight_msa[:, s] == 0) * subseq).sum() - sym_count
            subseq_records.append({'ppid': ppid, 'sym_count': sym_count, 'align_count': align_count})
        subseq_records = sorted(subseq_records, key=lambda x: x['sym_count'])

        count_sum = sum([ppid_weights[record['ppid']] * record['sym_count'] for record in subseq_records[:-mean_trim]])
        weight_sum = sum([ppid_weights[record['ppid']] for record in subseq_records[:-mean_trim]])
        mean = max(count_sum / weight_sum, mean_min)
        p = 1 / (mean + 1)  # Geometric distribution on support 0, 1, ...
        k = np.log(alpha) / np.log(1 - p) - 1  # Expression for minimum k to achieve alpha significance

        for record in subseq_records:
            ppid, sym_count, align_count = record['ppid'], record['sym_count'], record['align_count']
            if sym_count >= k:
                seq_slices[ppid].append(s)
                rows1.append({'OGid': OGid, 'ppid': ppid, 'start': s.start, 'stop': s.stop,
                              'sym_count': sym_count, 'align_count': align_count})

    # Identify state 2+3 trims
    posterior = df['3'].to_numpy(copy=True)
    slices = get_hull_slices(posterior, gradient, posterior_high1, posterior_low1, gradient_low)
    for s in slices:
        posterior[s] = 0
    posterior += df['2'].to_numpy()
    gradient = np.gradient(posterior)
    slices = get_trim_slices(profile, posterior, gradient, posterior_high2, posterior_low2, gradient_high, gradient_low)
    for s in slices:
        rows2.append({'OGid': OGid, 'colnum': len(input_msa[0]['seq']), 'start': s.start, 'stop': s.stop,
                      'posterior2': df.loc[s, '2'].sum(), 'posterior3': df.loc[s, '3'].sum()})
    region_slices = get_complement_slices(slices)  # Complement slices to get untrimmed regions

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
df1 = pd.DataFrame(rows1)
df1.to_csv('out/seq_stats.tsv', sep='\t', index=False)

df2 = pd.DataFrame(rows2)
df2.to_csv('out/region_stats.tsv', sep='\t', index=False)

df2['length'] = df2['stop'] - df2['start']
df2['length_ratio'] = df2['length'] / df2['colnum']
df2['norm2'] = df2['posterior2'] / df2['length']
df2['norm3'] = df2['posterior3'] / df2['length']

groups1 = df1.groupby('OGid')
groups2 = df2.groupby('OGid')

# Pie chart by presence of trims
union = set(df1['OGid']) | set(df2['OGid'])
intersection = set(df1['OGid']) & set(df2['OGid'])
OGids1 = set(df1['OGid']) - intersection
OGids2 = set(df2['OGid']) - intersection

values = [len(set(OGids)) - len(union), len(OGids1), len(OGids2), len(intersection)]
labels = [f'{label}\n({value:,})' for label, value in zip(['no trims', 'sequence trims only', 'region trims only', 'sequence and region trims'], values)]
colors = ['C0', 'C2', 'C1', 'C3']
fig, ax = plt.subplots(gridspec_kw={'bottom': 0.2})
ax.pie(values, labels=labels, labeldistance=None, colors=colors)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0), ncol=2)
ax.set_title('OGs by type of trims')
fig.savefig('out/pie_trims.png')
plt.close()

# Distribution of number of sequence trims in OGs
counts = groups1.size().value_counts()
fig, ax = plt.subplots()
ax.bar(counts.index, counts.values, width=1)
ax.set_xlabel('Number of sequence trims in OG')
ax.set_ylabel('Number of OGs')
fig.savefig('out/hist_OGnum-seqtrimnum.png')
plt.close()

# Distribution of number of sequences with trims in OGs
counts = groups1['ppid'].nunique().value_counts()
fig, ax = plt.subplots()
ax.bar(counts.index, counts.values, width=1)
ax.set_xlabel('Number of sequences with trims in OG')
ax.set_ylabel('Number of OGs')
fig.savefig('out/hist_OGnum-seqnum.png')
plt.close()

# Distribution of number of sequence trims by sequence
counts = df1.groupby(['OGid', 'ppid']).size().value_counts()
fig, ax = plt.subplots()
ax.bar(counts.index, counts.values, width=1)
ax.set_xlabel('Number of sequence trims within a sequence')
ax.set_ylabel('Number of sequences')
fig.savefig('out/hist_seqnum-seqtrimnum.png')
plt.close()

# Distribution of number of removed symbols in sequence trims
counts = df1['sym_count'].value_counts()
fig, ax = plt.subplots()
ax.bar(counts.index, counts.values, width=1)
ax.set_xlabel('Number of non-gap symbols in sequence trim')
ax.set_ylabel('Number of sequence trims')
fig.savefig('out/hist_seqtrimnum-symnum.png')
plt.close()

# Distribution of number of removed symbols in sequence trims (truncated)
idx = int(np.ceil(len(df1) * 0.95))
counts = df1['sym_count'].sort_values(ignore_index=True)[:idx].value_counts()
fig, ax = plt.subplots()
ax.bar(counts.index, counts.values, width=1)
ax.set_xlabel('Number of non-gap symbols in sequence trim')
ax.set_ylabel('Number of sequence trims')
fig.savefig('out/hist_seqtrimnum-symnum_95.png')
plt.close()

# Distribution of number of aligned symbols per non-gap symbols in sequence trims
fig, ax = plt.subplots(gridspec_kw={'bottom': 0.15})
ax.hist(df1['align_count'] / df1['sym_count'], bins=100)
ax.set_xlabel('Number of aligned non-gap symbols\nper non-gap symbol in trim')
ax.set_ylabel('Number of sequence trims')
fig.savefig('out/hist_seqtrimnum-symratio.png')
plt.close()

# Correlation of non-gap symbols with aligned symbols
idx = int(np.ceil(len(df1) * 0.95))
x = df1.sort_values(by='sym_count', ignore_index=True)[:idx]
fig, ax = plt.subplots()
hb = ax.hexbin(x['sym_count'], x['align_count'], gridsize=50, mincnt=1, linewidth=0)
ax.set_xlabel('Number of non-gap symbols in sequence trim')
ax.set_ylabel('Number of aligned non-gap symbols')
fig.colorbar(hb)
fig.savefig('out/hist_alignnum-symnum.png')
plt.close()

# Distribution of number of region trims in OGs
counts = groups2.size().value_counts()
fig, ax = plt.subplots()
ax.bar(counts.index, counts.values, width=1)
ax.set_xlabel('Number of region trims in OG')
ax.set_ylabel('Number of OGs')
fig.savefig('out/hist_OGnum-regionnum.png')
plt.close()

# Distribution of lengths of region trims
fig, ax = plt.subplots()
ax.hist(df2['length'], bins=100)
ax.set_xlabel('Length of region trim')
ax.set_ylabel('Number of region trims')
fig.savefig('out/hist_regionnum-regionlength.png')
plt.close()

# Distribution of lengths ratio of region trims
fig, ax = plt.subplots()
ax.hist(df2['length_ratio'], bins=50)
ax.set_xlabel('Length ratio of region trim')
ax.set_ylabel('Number of region trims')
fig.savefig('out/hist_regionnum-regionratio.png')
plt.close()

# Distribution of lengths of total region trims
fig, ax = plt.subplots()
ax.hist(groups2['length'].sum(), bins=100)
ax.set_xlabel('Total length of region trims in OG')
ax.set_ylabel('Number of OGs')
fig.savefig('out/hist_OGnum-regionlength.png')
plt.close()

# Distribution of length ratios of total region trims
fig, ax = plt.subplots()
ax.hist(groups2['length_ratio'].sum(), bins=100)
ax.set_xlabel('Total length ratio of region trims in OG')
ax.set_ylabel('Number of OGs')
fig.savefig('out/hist_OGnum-regionratio.png')
plt.close()

# Hexbin of length ratios vs number of region trims
fig = plt.figure(figsize=(6, 6), layout='constrained')
gs = fig.add_gridspec(4, 2, height_ratios=(1, 2, 0.15, 0.1), width_ratios=(4, 1))
ax = fig.add_subplot(gs[1, 0])
ax_histx = fig.add_subplot(gs[0, 0], sharex=ax)
ax_histy = fig.add_subplot(gs[1, 1], sharey=ax)

hb = ax.hexbin(groups2.size(), groups2['length_ratio'].sum(), bins='log', gridsize=50, mincnt=1, linewidth=0)
cax = fig.add_subplot(gs[3, 0])
fig.colorbar(hb, cax=cax, orientation='horizontal')

counts = groups2.size().value_counts()
ax_histx.bar(counts.index, counts.values, width=1)
ax_histy.hist(groups2['length_ratio'].sum(), bins=100, orientation='horizontal')

ax.set_xlabel('Number of region trims in OG')
ax.set_ylabel('Total length ratio of region trims in OG')
fig.savefig('out/hexbin_regionratio-regionnum.png')
plt.close()

# Hexbin of posterior2 vs posterior3 in region trims
fig, ax = plt.subplots()
hb = ax.hexbin(df2['norm2'], df2['norm3'], bins='log', gridsize=25, mincnt=1, linewidth=0)
ax.set_xlabel('Average state 2 posterior in region trim')
ax.set_ylabel('Average state 3 posterior in region trim')
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