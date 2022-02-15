"""Plot various statistics of the input genomes."""

import matplotlib.pyplot as plt
import os
import pandas as pd
import re
from itertools import groupby
from numpy import linspace


def get_Xmax(seq):
    if 'X' in seq.upper():
        Xlens = [len(list(group)) for key, group in groupby(seq.upper()) if key == 'X']
        return max(Xlens)
    else:
        return 0


ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+)'}

# Load seq metadata
ppid2gnid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, _, _ = line.split()
        ppid2gnid[ppid] = gnid

# Parse genomes
genomes = {}
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, prot_path = line.split()
        genomes[spid] = (source, prot_path)

# Parse polypeptides
sqid0 = 0
sqids = {}
rows = []
for spid, (source, prot_path) in genomes.items():
    with open(prot_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid0 = re.search(ppid_regex[source], line).group(1)
                gnid = ppid2gnid[ppid0]
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq0 = ''.join(seqlines)

            try:
                for sqid1, seq1 in sqids[gnid]:
                    if seq0 == seq1:
                        rows.append({'ppid': ppid0, 'gnid': gnid, 'spid': spid, 'sqid': sqid1,
                                     'sqlen': len(seq1), 'Xnum': seq1.upper().count('X'), 'Xmax': get_Xmax(seq1)})
                        break
                else:
                    sqids[gnid].append((str(sqid0).zfill(8), seq0))
                    rows.append({'ppid': ppid0, 'gnid': gnid, 'spid': spid, 'sqid': str(sqid0).zfill(8),
                                 'sqlen': len(seq0), 'Xnum': seq0.upper().count('X'), 'Xmax': get_Xmax(seq0)})
                    sqid0 += 1
            except KeyError:
                sqids[gnid] = [(str(sqid0).zfill(8), seq0)]
                rows.append({'ppid': ppid0, 'gnid': gnid, 'spid': spid, 'sqid': str(sqid0).zfill(8),
                             'sqlen': len(seq0), 'Xnum': seq0.upper().count('X'), 'Xmax': get_Xmax(seq0)})
                sqid0 += 1

# Make plots output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Define dataframe and groups
df = pd.DataFrame(rows)
ncbi_spids = [spid for spid in genomes if genomes[spid][0] == 'NCBI']
flybase_spids = [spid for spid in genomes if genomes[spid][0] == 'FlyBase']

# 1.1 Distribution of polypeptides across associated species
spid_ppidnum = df['spid'].value_counts()
ncbi = {spid: spid_ppidnum[spid] for spid in ncbi_spids}
flybase = {spid: spid_ppidnum[spid] for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI', color='C0', width=0.75)
plt.bar(range(len(ncbi), len(ncbi) + len(flybase)),
        flybase.values(), label=r'FlyBase', color='C1', width=0.75)
plt.xticks(list(range(len(ncbi) + len(flybase))), list(ncbi.keys()) + list(flybase.keys()), rotation=60, fontsize=8)
plt.xlabel('Associated species')
plt.ylabel('Number of polypeptides')
plt.title('Distribution of polypeptides across associated species')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_ppidnum-species.png')
plt.close()

# 1.2 Distribution of unique sequences across associated species
spid_sqidnum = df.groupby('spid')['sqid'].nunique()
ncbi = {spid: spid_sqidnum[spid] for spid in ncbi_spids}
flybase = {spid: spid_sqidnum[spid] for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI', color='C0', width=0.75)
plt.bar(range(len(ncbi), len(ncbi) + len(flybase)),
        flybase.values(), label=r'FlyBase', color='C1', width=0.75)
plt.xticks(list(range(len(ncbi) + len(flybase))), list(ncbi.keys()) + list(flybase.keys()), rotation=60, fontsize=8)
plt.xlabel('Associated species')
plt.ylabel('Number of unique sequences')
plt.title('Distribution of unique sequences across associated species')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_sqidnum-species.png')
plt.close()

# 1.3 Distribution of redundant sequences across associated species
spid_rsqidnum = spid_ppidnum - spid_sqidnum
ncbi = {spid: spid_rsqidnum.get(spid, 0) for spid in ncbi_spids}
flybase = {spid: spid_rsqidnum.get(spid, 0) for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI', color='C0', width=0.75)
plt.bar(range(len(ncbi), len(ncbi) + len(flybase)),
        flybase.values(), label=r'FlyBase', color='C1', width=0.75)
plt.xticks(list(range(len(ncbi) + len(flybase))), list(ncbi.keys()) + list(flybase.keys()), rotation=60, fontsize=8)
plt.xlabel('Associated species')
plt.ylabel('Number of redundant sequences')
plt.title('Distribution of redundant sequences across associated species')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_rsqidnum-species.png')
plt.close()

# 1.4 Distribution of number of genes across associated species
spid_gnidnum = df.groupby('spid')['gnid'].nunique()
ncbi = {spid: spid_gnidnum[spid] for spid in ncbi_spids}
flybase = {spid: spid_gnidnum[spid] for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI', color='C0', width=0.75)
plt.bar(range(len(ncbi), len(ncbi) + len(flybase)),
        flybase.values(), label='FlyBase', color='C1', width=0.75)
plt.xticks(list(range(len(ncbi) + len(flybase))), list(ncbi.keys()) + list(flybase.keys()), rotation=60, fontsize=8)
plt.xlabel('Associated species')
plt.ylabel('Number of genes')
plt.ylim((0, 1.1 * plt.ylim()[1]))  # Increase y-axis span to prevent overlap with legend
plt.title('Distribution of genes across associated species')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_gnidnum-species.png')
plt.close()

# 2 Distribution of polypeptides across length
ncbi = df.loc[df['spid'].isin(ncbi_spids), 'sqlen']
flybase = df.loc[df['spid'].isin(flybase_spids), 'sqlen']

# NCBI and FlyBase
bins = linspace(min(ncbi.min(), flybase.min()), max(ncbi.max(), flybase.max()), 100, endpoint=True)
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].hist(ncbi, bins=bins, label='NCBI', color='C0')
axs[1].hist(flybase, bins=bins, label='FlyBase', color='C1')
axs[1].set_xlabel('Length')
fig.suptitle('Distribution of polypeptides across length')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.set_ylabel('Number of polypeptides')
    ax.legend()
fig.savefig('out/hist_ppidnum-sqlen.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (flybase, 'FlyBase', 'C1', 'FB')]:
    plt.hist(data, bins=100, label=data_label, color=color)
    plt.xlabel('Length')
    plt.ylabel('Number of polypeptides')
    plt.title('Distribution of polypeptides across length')
    plt.legend()
    plt.savefig(f'out/hist_ppidnum-sqlen_{file_label}.png')
    plt.close()

# Calculate counts grouped by gene
spidgnid_pairs = df[['spid', 'gnid']].drop_duplicates()
ppidnum = df.groupby('gnid').size().rename('ppidnum')
sqidnum = df.groupby('gnid')['sqid'].nunique().rename('sqidnum')
gnid_nums = spidgnid_pairs.join(ppidnum, on='gnid').join(sqidnum, on='gnid')
gnid_nums.to_csv('out/gnid_nums.tsv', sep='\t', index=False)

# 3.1 Distribution of genes across number of associated polypeptides
ncbi = gnid_nums.loc[gnid_nums['spid'].isin(ncbi_spids), 'ppidnum'].value_counts()
flybase = gnid_nums.loc[gnid_nums['spid'].isin(flybase_spids), 'ppidnum'].value_counts()

# NCBI and FlyBase
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].bar(ncbi.index, ncbi.values, label='NCBI', color='C0', width=1)
axs[1].bar(flybase.index, flybase.values, label='FlyBase', color='C1', width=1)
axs[1].set_xlabel('Number of associated polypeptides')
fig.suptitle('Distribution of genes across\nnumber of associated polypeptides')
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
    plt.xlabel('Number of associated polypeptides')
    plt.ylabel('Number of genes')
    plt.title('Distribution of genes across\nnumber of associated polypeptides')
    plt.legend()
    plt.savefig(f'out/hist_gnidnum-ppidnum_{file_label}.png')
    plt.close()

# 3.2 Distribution of genes across number of associated unique sequences
ncbi = gnid_nums.loc[gnid_nums['spid'].isin(ncbi_spids), 'sqidnum'].value_counts()
flybase = gnid_nums.loc[gnid_nums['spid'].isin(flybase_spids), 'sqidnum'].value_counts()

# NCBI and FlyBase
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].bar(ncbi.index, ncbi.values, label='NCBI', color='C0', width=1)
axs[1].bar(flybase.index, flybase.values, label='FlyBase', color='C1', width=1)
axs[1].set_xlabel('Number of unique sequences')
fig.suptitle('Distribution of genes across\nnumber of associated unique sequences')
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
    plt.xlabel('Number of associated unique sequences')
    plt.ylabel('Number of genes')
    plt.title('Distribution of genes across\nnumber of associated unique sequences')
    plt.legend()
    plt.savefig(f'out/hist_gnidnum-sqidnum_{file_label}.png')
    plt.close()

# 3.3 Distribution of genes across number of associated redundant sequences
gnid_nums['rsqidnum'] = gnid_nums['ppidnum'] - gnid_nums['sqidnum']
ncbi = gnid_nums.loc[gnid_nums['spid'].isin(ncbi_spids), 'rsqidnum'].value_counts()
flybase = gnid_nums.loc[gnid_nums['spid'].isin(flybase_spids), 'rsqidnum'].value_counts()

# NCBI and FlyBase
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].bar(ncbi.index, ncbi.values, label='NCBI', color='C0', width=1)
axs[1].bar(flybase.index, flybase.values, label='FlyBase', color='C1', width=1)
axs[1].set_xlabel('Number of redundant sequences')
fig.suptitle('Distribution of genes across\nnumber of associated redundant sequences')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.set_ylabel('Number of genes')
    ax.legend()
fig.savefig('out/hist_gnidnum-rsqidnum.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (flybase, 'FlyBase', 'C1', 'FB')]:
    plt.bar(data.index, data.values, label=data_label, color=color, width=1)
    plt.xlabel('Number of associated redundant sequences')
    plt.ylabel('Number of genes')
    plt.title('Distribution of genes across\nnumber of associated redundant sequences')
    plt.legend()
    plt.savefig(f'out/hist_gnidnum-rsqidnum_{file_label}.png')
    plt.close()

# 4 Unknown amino acid plots
dfX = df[df['Xnum'] > 0].drop('ppid', axis=1).drop_duplicates()

# 4.1 Distribution of number of unknown amino acids
counts = dfX['Xnum'].value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of unknown amino acids')
plt.ylabel('Number of unique sequences')
plt.title('Distribution of unique sequences across\nnumber of unknown amino acids')
plt.savefig('out/hist_sqidnum-Xnum.png')
plt.close()

# 4.2 Distribution of length of largest segment of unknown amino acids
counts = dfX['Xmax'].value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Length of largest segment of unknown amino acids')
plt.ylabel('Number of unique sequences')
plt.title('Distribution of unique sequences across\nlength of largest segment of unknown amino acids')
plt.savefig('out/hist_sqidnum-Xmax.png')
plt.close()

# 4.3 Distribution of fraction of unknown amino acids
plt.hist(dfX['Xmax'] / dfX['sqlen'], bins=50)
plt.xlabel('Fraction of unknown amino acids')
plt.ylabel('Number of unique sequences')
plt.title('Distribution of unique sequences across\nfraction of unknown amino acids')
plt.savefig('out/hist_sqidnum-Xfrac.png')
plt.close()

# 4.4 Distribution of unique sequences with unknown amino acids across associated species
spid_Xsqidnum = dfX['spid'].value_counts()
ncbi = {spid: spid_Xsqidnum.get(spid, 0) for spid in ncbi_spids}
flybase = {spid: spid_Xsqidnum.get(spid, 0) for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI', color='C0', width=0.75)
plt.bar(range(len(ncbi), len(ncbi) + len(flybase)),
        flybase.values(), label=r'FlyBase', color='C1', width=0.75)
plt.xticks(list(range(len(ncbi) + len(flybase))), list(ncbi.keys()) + list(flybase.keys()), rotation=60, fontsize=8)
plt.xlabel('Associated species')
plt.ylabel('Number of unique sequences with unknown amino acids')
plt.title('Distribution of unique sequences with\nunknown amino acids across associated species')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_Xsqid-species.png')
plt.close()

# 4.5 Distribution of genes with unknown amino acids across associated species
spid_Xgnidnum = dfX.groupby('spid')['gnid'].nunique()
ncbi = {spid: spid_Xgnidnum.get(spid, 0) for spid in ncbi_spids}
flybase = {spid: spid_Xgnidnum.get(spid, 0) for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI', color='C0', width=0.75)
plt.bar(range(len(ncbi), len(ncbi) + len(flybase)),
        flybase.values(), label=r'FlyBase', color='C1', width=0.75)
plt.xticks(list(range(len(ncbi) + len(flybase))), list(ncbi.keys()) + list(flybase.keys()), rotation=60, fontsize=8)
plt.xlabel('Associated species')
plt.ylabel('Number of genes with unknown amino acids')
plt.title('Distribution of genes with \nunknown amino acids across associated species')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_Xgnid-species1.png')
plt.close()

# 4.5 Distribution of genes with unknown amino acids across associated species
gnid_Xval = df.groupby('gnid')['Xnum'].agg(['min', 'max'])
spidgnid_pairs = df[['spid', 'gnid']].drop_duplicates()
spid_Xgnidnum = gnid_Xval[(gnid_Xval['min'] > 0) & (gnid_Xval['max'] > 0)].join(spidgnid_pairs.set_index('gnid'))['spid'].value_counts()
ncbi = {spid: spid_Xgnidnum.get(spid, 0) for spid in ncbi_spids}
flybase = {spid: spid_Xgnidnum.get(spid, 0) for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI', color='C0', width=0.75)
plt.bar(range(len(ncbi), len(ncbi) + len(flybase)),
        flybase.values(), label=r'FlyBase', color='C1', width=0.75)
plt.xticks(list(range(len(ncbi) + len(flybase))), list(ncbi.keys()) + list(flybase.keys()), rotation=60, fontsize=8)
plt.xlabel('Associated species')
plt.ylabel('Number of genes with unknown amino acids')
plt.title('Distribution of genes with \nunknown amino acids across associated species')
plt.tight_layout()
plt.legend()
plt.savefig('out/bar_Xgnid-species2.png')
plt.close()

# 4.6 Correlation of number with length of largest segment
plt.scatter(dfX['Xnum'], dfX['Xmax'], alpha=0.5, s=10, edgecolors='none')
plt.xlabel('Number of unknown amino acids')
plt.ylabel('Length of largest segment of unknown amino acids')
plt.savefig('out/scatter_Xmax-Xnum.png')
plt.close()

# 4.7 Counts
gnid_Xval = df.groupby('gnid')['Xnum'].agg(['min', 'max'])
s1 = sum(gnid_Xval['min'] == 0)
s2 = sum(gnid_Xval['max'] > 0)
s3 = sum((gnid_Xval['min'] == 0) & (gnid_Xval['max'] > 0))

print('Fraction of sequences with unknown amino acids:', len(dfX) / len(df))
print()
print('Genes with at least one sequence without unknown amino acids')
print('Number:', s1)
print('Fraction:', s1 / len(gnid_Xval))
print()
print('Genes with at least one sequence with unknown amino acids')
print('Number:', s2)
print('Fraction:', s2 / len(gnid_Xval))
print()
print('Genes with at least one sequence without unknown amino acids and at least one sequence with unknown amino acids')
print('Number:', s3)
print('Fraction (all genes):', s3 / len(gnid_Xval))
print('Fraction (genes with unknown amino acids):', s3 / s2)

"""
OUTPUT
Fraction of sequences with unknown amino acids: 0.013849222288397346

Genes with at least one sequence without unknown amino acids
Number: 434011
Fraction: 0.9787918766843702

Genes with at least one sequence with unknown amino acids
Number: 9833
Fraction: 0.022175614266544883

Genes with at least one sequence without unknown amino acids and at least one sequence with unknown amino acids
Number: 429
Fraction (all genes): 0.000967490950915057
Fraction (genes with unknown amino acids): 0.043628597579578966

DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../config/genomes.tsv
"""