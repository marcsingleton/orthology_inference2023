"""Plot various statistics of the input genomes."""

import os
import re
from itertools import groupby

import matplotlib.pyplot as plt
import pandas as pd
from numpy import linspace
from src.utils import read_fasta


def get_Xmax(seq):
    if 'X' in seq:
        Xlens = [len(list(group)) for key, group in groupby(seq) if key == 'X']
        return max(Xlens)
    else:
        return 0


ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+)'}
alphabet = {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X'}

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, _, sqid = line.split()
        ppid2meta[ppid] = (gnid, sqid)

# Parse genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, prot_path, _ = line.split()
        genomes.append((spid, source, prot_path))

# Parse proteins
rows, syms = [], set()
for spid, source, prot_path in genomes:
    fasta = read_fasta(prot_path)
    for header, seq in fasta:
        # Extract data
        seq = seq.upper()
        ppid = re.search(ppid_regex[source], header).group(1)
        gnid, sqid = ppid2meta[ppid]
        rows.append({'ppid': ppid, 'gnid': gnid, 'spid': spid, 'sqid': sqid, 'source': source,
                     'seqlen': len(seq), 'Xnum': seq.count('X'), 'Xmax': get_Xmax(seq)})

        # Check for ambiguous symbols
        syms |= set(seq) - alphabet
print('Out-of-alphabet symbols detected:', syms)
print()

# Make plots output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Define dataframe and groups
df = pd.DataFrame(rows)

spids = [spid for spid, _, _ in sorted(genomes)]
ncbi_spids = [spid for spid, source, _ in genomes if source == 'NCBI']
ncbi_positions = [i for i, (_, source, _) in enumerate(sorted(genomes)) if source == 'NCBI']
flybase_spids = [spid for spid, source, _ in genomes if source == 'FlyBase']
flybase_positions = [i for i, (_, source, _) in enumerate(sorted(genomes)) if source == 'FlyBase']

# 1.1 Distribution of proteins across species
spid2ppidnum = df['spid'].value_counts()
ncbi_values = [spid2ppidnum[spid] for spid in ncbi_spids]
flybase_values = [spid2ppidnum[spid] for spid in flybase_spids]

plt.bar(ncbi_positions, ncbi_values, label='NCBI', color='C0', width=0.75)
plt.bar(flybase_positions, flybase_values, label='FlyBase', color='C1', width=0.75)
plt.xticks(range(len(spids)), spids, rotation=60, fontsize=8)
plt.xlabel('Species')
plt.ylabel('Number of proteins')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_ppidnum-species.png')
plt.close()

# 1.2 Distribution of unique sequences across species
spid2sqidnum = df.groupby('spid')['sqid'].nunique()
ncbi_values = [spid2sqidnum[spid] for spid in ncbi_spids]
flybase_values = [spid2sqidnum[spid] for spid in flybase_spids]

plt.bar(ncbi_positions, ncbi_values, label='NCBI', color='C0', width=0.75)
plt.bar(flybase_positions, flybase_values, label='FlyBase', color='C1', width=0.75)
plt.xticks(range(len(spids)), spids, rotation=60, fontsize=8)
plt.xlabel('Species')
plt.ylabel('Number of unique sequences')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_sqidnum-species.png')
plt.close()

# 1.3 Distribution of redundant sequences across species
spid2redundant = spid2ppidnum - spid2sqidnum
ncbi_values = [spid2redundant.get(spid, 0) for spid in ncbi_spids]
flybase_values = [spid2redundant.get(spid, 0) for spid in flybase_spids]

plt.bar(ncbi_positions, ncbi_values, label='NCBI', color='C0', width=0.75)
plt.bar(flybase_positions, flybase_values, label='FlyBase', color='C1', width=0.75)
plt.xticks(range(len(spids)), spids, rotation=60, fontsize=8)
plt.xlabel('Species')
plt.ylabel('Number of redundant sequences')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_redundant-species.png')
plt.close()

# 1.4 Distribution of number of genes across species
spid2gnidnum = df.groupby('spid')['gnid'].nunique()
ncbi_values = [spid2gnidnum[spid] for spid in ncbi_spids]
flybase_values = [spid2gnidnum[spid] for spid in flybase_spids]

plt.bar(ncbi_positions, ncbi_values, label='NCBI', color='C0', width=0.75)
plt.bar(flybase_positions, flybase_values, label='FlyBase', color='C1', width=0.75)
plt.xticks(range(len(spids)), spids, rotation=60, fontsize=8)
plt.xlabel('Species')
plt.ylabel('Number of genes')
plt.ylim((0, 1.1 * plt.ylim()[1]))  # Increase y-axis span to prevent overlap with legend
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_gnidnum-species.png')
plt.close()

# 2 Distribution of proteins across length
ncbi = df.loc[df['source'] == 'NCBI', 'seqlen']
flybase = df.loc[df['source'] == 'FlyBase', 'seqlen']

# NCBI and FlyBase
bins = linspace(min(ncbi.min(), flybase.min()), max(ncbi.max(), flybase.max()), 100, endpoint=True)
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].hist(ncbi, bins=bins, label='NCBI', color='C0')
axs[1].hist(flybase, bins=bins, label='FlyBase', color='C1')
axs[1].set_xlabel('Length of protein')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.set_ylabel('Number of proteins')
    ax.legend()
fig.savefig('out/hist_ppidnum-seqlen.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (flybase, 'FlyBase', 'C1', 'FB')]:
    plt.hist(data, bins=100, label=data_label, color=color)
    plt.xlabel('Length of protein')
    plt.ylabel('Number of proteins')
    plt.legend()
    plt.savefig(f'out/hist_ppidnum-seqlen_{file_label}.png')
    plt.close()

# Calculate counts grouped by gene
gnidnums = df.groupby(['spid', 'gnid'])[['ppid', 'sqid']].nunique().rename(columns={'ppid': 'ppidnum', 'sqid': 'sqidnum'}).reset_index()
gnidnums.to_csv('out/gnidnums.tsv', sep='\t', index=False)

# 3.1 Distribution of genes across number of proteins
ncbi = gnidnums.loc[gnidnums['spid'].isin(set(ncbi_spids)), 'ppidnum'].value_counts()
flybase = gnidnums.loc[gnidnums['spid'].isin(set(flybase_spids)), 'ppidnum'].value_counts()

# NCBI and FlyBase
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].bar(ncbi.index, ncbi.values, label='NCBI', color='C0', width=1)
axs[1].bar(flybase.index, flybase.values, label='FlyBase', color='C1', width=1)
axs[1].set_xlabel('Number of proteins')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.set_ylabel('Number of genes')
    ax.legend()
fig.savefig('out/hist_gnidnum-ppidnum.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (flybase, 'FlyBase', 'C1', 'FB')]:
    plt.bar(data.index, data.values, label=data_label, color=color, width=1)
    plt.xlabel('Number of proteins')
    plt.ylabel('Number of genes')
    plt.legend()
    plt.savefig(f'out/hist_gnidnum-ppidnum_{file_label}.png')
    plt.close()

# 3.2 Distribution of genes across number of unique sequences
ncbi = gnidnums.loc[gnidnums['spid'].isin(set(ncbi_spids)), 'sqidnum'].value_counts()
flybase = gnidnums.loc[gnidnums['spid'].isin(set(flybase_spids)), 'sqidnum'].value_counts()

# NCBI and FlyBase
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].bar(ncbi.index, ncbi.values, label='NCBI', color='C0', width=1)
axs[1].bar(flybase.index, flybase.values, label='FlyBase', color='C1', width=1)
axs[1].set_xlabel('Number of unique sequences')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.set_ylabel('Number of genes')
    ax.legend()
fig.savefig('out/hist_gnidnum-sqidnum.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (flybase, 'FlyBase', 'C1', 'FB')]:
    plt.bar(data.index, data.values, label=data_label, color=color, width=1)
    plt.xlabel('Number of unique sequences')
    plt.ylabel('Number of genes')
    plt.legend()
    plt.savefig(f'out/hist_gnidnum-sqidnum_{file_label}.png')
    plt.close()

# 3.3 Distribution of genes across number of redundant sequences
gnidnums['redundant'] = gnidnums['ppidnum'] - gnidnums['sqidnum']
ncbi = gnidnums.loc[gnidnums['spid'].isin(set(ncbi_spids)), 'redundant'].value_counts()
flybase = gnidnums.loc[gnidnums['spid'].isin(set(flybase_spids)), 'redundant'].value_counts()

# NCBI and FlyBase
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].bar(ncbi.index, ncbi.values, label='NCBI', color='C0', width=1)
axs[1].bar(flybase.index, flybase.values, label='FlyBase', color='C1', width=1)
axs[1].set_xlabel('Number of redundant sequences')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.set_ylabel('Number of genes')
    ax.legend()
fig.savefig('out/hist_gnidnum-redundant.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (flybase, 'FlyBase', 'C1', 'FB')]:
    plt.bar(data.index, data.values, label=data_label, color=color, width=1)
    plt.xlabel('Number of redundant sequences')
    plt.ylabel('Number of genes')
    plt.legend()
    plt.savefig(f'out/hist_gnidnum-redundant_{file_label}.png')
    plt.close()

# 4 Unknown amino acid plots
dfX = df[df['Xnum'] > 0].drop('ppid', axis=1).drop_duplicates()  # Use only unique sequences

# 4.1 Distribution of number of unknown amino acids
counts = dfX['Xnum'].value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of unknown amino acids in sequence')
plt.ylabel('Number of unique sequences')
plt.savefig('out/hist_sqidnum-Xnum.png')
plt.close()

# 4.2 Distribution of length of largest segment of unknown amino acids
counts = dfX['Xmax'].value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Length of largest segment of unknown amino acids in sequence')
plt.ylabel('Number of unique sequences')
plt.savefig('out/hist_sqidnum-Xmax.png')
plt.close()

# 4.3 Distribution of fraction of unknown amino acids
plt.hist(dfX['Xmax'] / dfX['seqlen'], bins=50)
plt.xlabel('Fraction of unknown amino acids in sequence')
plt.ylabel('Number of unique sequences')
plt.savefig('out/hist_sqidnum-Xfrac.png')
plt.close()

# 4.4 Distribution of unique sequences with unknown amino acids across species
spid2Xsqidnum = dfX['spid'].value_counts()
ncbi_values = [spid2Xsqidnum.get(spid, 0) for spid in ncbi_spids]
flybase_values = [spid2Xsqidnum.get(spid, 0) for spid in flybase_spids]

plt.bar(ncbi_positions, ncbi_values, label='NCBI', color='C0', width=0.75)
plt.bar(flybase_positions, flybase_values, label='FlyBase', color='C1', width=0.75)
plt.xticks(range(len(spids)), spids, rotation=60, fontsize=8)
plt.xlabel('Species')
plt.ylabel('Number of unique sequences with unknown amino acids')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_Xsqid-species.png')
plt.close()

# 4.5 Distribution of genes with unknown amino acids across species
spid2Xgnidnum = dfX.groupby('spid')['gnid'].nunique()
ncbi_values = [spid2Xgnidnum.get(spid, 0) for spid in ncbi_spids]
flybase_values = [spid2Xgnidnum.get(spid, 0) for spid in flybase_spids]

plt.bar(ncbi_positions, ncbi_values, label='NCBI', color='C0', width=0.75)
plt.bar(flybase_positions, flybase_values, label='FlyBase', color='C1', width=0.75)
plt.xticks(range(len(spids)), spids, rotation=60, fontsize=8)
plt.xlabel('Species')
plt.ylabel('Number of genes with unknown amino acids')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_Xgnid-species.png')
plt.close()

# 4.5 Correlation of number with length of largest segment
plt.scatter(dfX['Xnum'], dfX['Xmax'], alpha=0.5, s=10, edgecolors='none')
plt.xlabel('Number of unknown amino acids')
plt.ylabel('Length of largest segment of unknown amino acids')
plt.savefig('out/scatter_Xmax-Xnum.png')
plt.close()

# 4.6 Counts
gnid2Xval = df.groupby('gnid')['Xnum'].agg(['min', 'max'])
s1 = sum(gnid2Xval['min'] == 0)
s2 = sum(gnid2Xval['max'] > 0)
s3 = sum((gnid2Xval['min'] == 0) & (gnid2Xval['max'] > 0))

print('Fraction of sequences with unknown amino acids:', round(len(dfX) / len(df), 3))
print()
print('Genes with at least one sequence without unknown amino acids')
print('Number:', s1)
print('Fraction:', round(s1 / len(gnid2Xval), 3))
print()
print('Genes with at least one sequence with unknown amino acids')
print('Number:', s2)
print('Fraction:', round(s2 / len(gnid2Xval), 3))
print()
print('Genes with at least one sequence without unknown amino acids and at least one sequence with unknown amino acids')
print('Number:', s3)
print('Fraction (all genes):', round(s3 / len(gnid2Xval), 3))
print('Fraction (genes with unknown amino acids):', round(s3 / s2, 3))

"""
Out-of-alphabet symbols detected: {'U'}

Fraction of sequences with unknown amino acids: 0.014

Genes with at least one sequence without unknown amino acids
Number: 447551
Fraction: 0.979

Genes with at least one sequence with unknown amino acids
Number: 9948
Fraction: 0.022

Genes with at least one sequence without unknown amino acids and at least one sequence with unknown amino acids
Number: 429
Fraction (all genes): 0.001
Fraction (genes with unknown amino acids): 0.043

DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../config/genomes.tsv
"""