"""Plot statistics of filtered OGs."""

import os
import re

import matplotlib.pyplot as plt
import pandas as pd
from numpy import linspace


def load_msa(path):
    msa = []
    with open(path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                header = line.rstrip()
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq = ''.join(seqlines)
            msa.append((header, seq))
    return msa


ppid_regex = r'ppid=([A-Za-z0-9_]+)'

# Parse segments
rows = []
with open('../aucpred_filter/out/segments.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, start, stop, disorder, ppids = line.split()

        msa = load_msa(f'../insertion_trim/out/{OGid}.mfa')
        msa = {re.search(ppid_regex, header).group(1): seq for header, seq in msa}

        for ppid in ppids.split(','):
            length = len([sym for sym in msa[ppid] if sym not in ['-', '.']])
            rows.append({'OGid': OGid, 'start': int(start), 'stop': int(stop), 'disorder': disorder == 'True', 'ppid': ppid, 'length': length})

# Convert to pandas objects
df = pd.DataFrame(rows)
OGs = df.groupby(['OGid', 'start', 'stop', 'disorder'])
idx = pd.IndexSlice

if not os.path.exists('out/'):
    os.mkdir('out/')


# Mean region length histogram
lengths = OGs['length'].aggregate(['mean', 'std'])

fig, axs = plt.subplots(2, 1, sharex=True)
xmin, xmax = lengths['mean'].min(), lengths['mean'].max()
axs[0].hist(lengths.loc[idx[:, :, :, True], 'mean'], bins=linspace(xmin, xmax, 100), color='C0', label='disorder')
axs[1].hist(lengths.loc[idx[:, :, :, False], 'mean'], bins=linspace(xmin, xmax, 100), color='C1', label='order')
axs[1].set_xlabel('Average length of region')
for i in range(2):
    axs[i].set_ylabel('Number of regions')
    axs[i].legend()
plt.savefig('out/hist_numregions-length.png')
plt.close()

# Number of sequences in region bar plot
fig, ax = plt.subplots()
counts1 = OGs.size()[idx[:, :, :, True]].value_counts()
counts2 = OGs.size()[idx[:, :, :, False]].value_counts()
ax.bar(counts1.index - 0.35/2, counts1.values, label='disorder', width=0.35)
ax.bar(counts2.index + 0.35/2, counts2.values, label='order', width=0.35)
ax.set_xlabel('Number of sequences in region')
ax.set_ylabel('Number of regions')
ax.legend()
plt.savefig('out/bar_numregions-numseqs.png')
plt.close()

# Counts of regions and unique OGs in each class
disorder = df[df['disorder']]
order = df[~df['disorder']]

plt.bar([0, 1], [len(disorder[['OGid', 'start', 'stop']].drop_duplicates()), len(order[['OGid', 'start', 'stop']].drop_duplicates())],
        tick_label=['disorder', 'order'], color=['C0', 'C1'], width=0.35)
plt.xlim((-0.5, 1.5))
plt.ylabel('Number of regions')
plt.savefig('out/bar_numregions-DO.png')
plt.close()

plt.bar([0, 1], [len(disorder['OGid'].drop_duplicates()), len(order['OGid'].drop_duplicates())],
        tick_label=['disorder', 'order'], color=['C0', 'C1'], width=0.35)
plt.xlim((-0.5, 1.5))
plt.ylabel('Number of unique OGs')
plt.savefig('out/bar_numOGs-DO.png')
plt.close()

"""
DEPENDENCIES
../aucpred_filter/aucpred_filter.py
    ../aucpred_filter/out/segments.tsv
../insertion_trim/extract.py
    ../insertion_trim/out/*.mfa
"""