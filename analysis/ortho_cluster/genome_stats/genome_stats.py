"""Plot various statistics of the input genomes."""

import matplotlib.pyplot as plt
import os
import pandas as pd
import re
from numpy import linspace

pp_regex = {'FlyBase': r'>(FBpp[0-9]+)',
            'NCBI': r'(XP_[0-9]+(\.[0-9]+)?)',
            'YO': r'(YOtr[A-Z]{2}[0-9]+\|orf[0-9]+)'}
gn_regex = {'FlyBase': r'parent=(FBgn[0-9]+)',
            'NCBI': r'gene=LOC([0-9]+)',
            'YO': r'>(YOgn[A-Z]{2}[0-9]+)'}
colors = {'FlyBase': 'C1',
          'NCBI': 'C0',
          'YO': 'C2'}

# Parse parameters
params = {}
with open('params.tsv') as infile:
    fields = infile.readline().split()  # Skip header
    for line in infile:
        spid, _, source, tcds_path = line.split()
        params[spid] = (source, tcds_path)

# Parse polypeptides
rows = []
for spid, (source, tcds_path) in params.items():
    with open(tcds_path) as infile:
        line = infile.readline()
        while line:
            if line.startswith('>'):
                ppid = re.search(pp_regex[source], line).group(1)
                gnid = re.search(gn_regex[source], line).group(1)
                line = infile.readline()

            pplen = 0
            while line and not line.startswith('>'):
                pplen += len(line.rstrip())
                line = infile.readline()

            rows.append({'ppid': ppid, 'gnid': gnid, 'spid': spid, 'pplen': pplen})

# Make plots output directory
if not os.path.exists(f'out/'):
    os.makedirs(f'out/')  # Recursive folder creation

# Define dataframe and groups
df = pd.DataFrame(rows)
ncbi_spids = [spid for spid in params if params[spid][0] == 'NCBI']
yo_spids = [spid for spid in params if params[spid][0] == 'YO']
flybase_spids = [spid for spid in params if params[spid][0] == 'FlyBase']

# Distribution of polypeptides across associated species
spid_ppidnum = df['spid'].value_counts()
ncbi = {spid: spid_ppidnum[spid] for spid in ncbi_spids}
yo = {spid: spid_ppidnum[spid] for spid in yo_spids}
flybase = {spid: spid_ppidnum[spid] for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI')
plt.bar(range(len(ncbi), len(ncbi) + len(yo)),
        yo.values(), label=r'Yang $et\ al.$')
plt.bar(range(len(ncbi) + len(yo), len(ncbi) + len(yo) + len(flybase)),
        flybase.values(), label=r'FlyBase')
plt.xticks(list(range(len(ncbi) + len(yo) + len(flybase))), list(ncbi.keys()) + list(yo.keys()) + list(flybase.keys()))
plt.xlabel('Associated species')
plt.ylabel('Number of polypeptides')
plt.title('Distribution of polypeptides across associated species')
plt.legend()
plt.savefig('out/dist_ppidnum-species_all.png')
plt.close()

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI')
plt.bar(range(len(ncbi), len(ncbi) + len(flybase)),
        flybase.values(), label=r'FlyBase', color='C2')
plt.xticks(list(range(len(ncbi) + len(flybase))), list(ncbi.keys()) + list(flybase.keys()))
plt.xlabel('Associated species')
plt.ylabel('Number of polypeptides')
plt.title('Distribution of polypeptides across associated species')
plt.legend()
plt.savefig('out/dist_ppidnum-species_NCBI-FB.png')
plt.close()

# Distribution of number of genes across associated species
spid_gnidnum = df[['spid', 'gnid']].drop_duplicates().groupby('spid').size()
ncbi = {spid: spid_gnidnum[spid] for spid in ncbi_spids}
yo = {spid: spid_gnidnum[spid] for spid in yo_spids}
flybase = {spid: spid_gnidnum[spid] for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI')
plt.bar(range(len(ncbi), len(ncbi) + len(yo)),
        yo.values(), label=r'Yang $et\ al.$')
plt.bar(range(len(ncbi) + len(yo), len(ncbi) + len(yo) + len(flybase)),
        flybase.values(), label=r'FlyBase')
plt.xticks(list(range(len(ncbi) + len(yo) + len(flybase))), list(ncbi.keys()) + list(yo.keys()) + list(flybase.keys()))
plt.xlabel('Associated species')
plt.ylabel('Number of genes')
plt.ylim((0, 1.1 * plt.ylim()[1]))  # Increase y-axis span to prevent overlap with legend
plt.title('Distribution of genes across associated species')
plt.legend()
plt.savefig('out/dist_gnidnum-species.png')
plt.close()

# Distribution of polypeptides across length
ncbi = df.loc[df['spid'].isin(ncbi_spids), 'pplen']
yo = df.loc[df['spid'].isin(yo_spids), 'pplen']
flybase = df.loc[df['spid'].isin(flybase_spids), 'pplen']

# All
bins = linspace(min(df['pplen']), max(df['pplen']), 150, endpoint=True)
fig, axs = plt.subplots(3, 1, sharex=True, figsize=(4.8, 6))
axs[0].hist(ncbi, bins=bins, label='NCBI', color='C0')
axs[1].hist(yo, bins=bins, label=r'Yang $et\ al.$', color='C1')
axs[2].hist(flybase, bins=bins, label='FlyBase', color='C2')
axs[2].set_xlabel('Length')
axs[1].set_ylabel('Number of polypeptides', labelpad=15)  # Adjust padding since center ylabels are narrow
fig.suptitle('Distribution of polypeptides across length')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.legend()
fig.savefig('out/dist_ppidnum-pplen_all.png')
plt.close()

# NCBI and FlyBase
bins = linspace(min(ncbi.min(), flybase.min()), max(ncbi.max(), flybase.max()), 100, endpoint=True)
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].hist(ncbi, bins=bins, label='NCBI', color='C0')
axs[1].hist(flybase, bins=bins, label='FlyBase', color='C2')
axs[1].set_xlabel('Length')
fig.suptitle('Distribution of polypeptides across length')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.set_ylabel('Number of polypeptides')
    ax.legend()
fig.savefig('out/dist_ppidnum-pplen_NCBI-FB.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (yo, r'Yang $et\ al.$', 'C1', 'YO'),
                                            (flybase, 'FlyBase', 'C2', 'FB')]:
    plt.hist(ncbi, bins=100, label=data_label, color=color)
    plt.xlabel('Length')
    plt.ylabel('Number of polypeptides')
    plt.title('Distribution of polypeptides across length')
    plt.legend()
    plt.savefig(f'out/dist_ppidnum-pplen_{file_label}.png')
    plt.close()

# Distribution of genes across number of associated polypeptides
gnid_ppidnum = df[['spid', 'gnid']].drop_duplicates().join(df.groupby('gnid').size().rename('ppidnum'), on='gnid')
gnid_ppidnum.to_csv('out/gnid_ppidnum.tsv', sep='\t', index=False)
ncbi = gnid_ppidnum.loc[gnid_ppidnum['spid'].isin(ncbi_spids), 'ppidnum']
yo = gnid_ppidnum.loc[gnid_ppidnum['spid'].isin(yo_spids), 'ppidnum']
flybase = gnid_ppidnum.loc[gnid_ppidnum['spid'].isin(flybase_spids), 'ppidnum']

# All
bins = linspace(min(gnid_ppidnum['ppidnum']), max(gnid_ppidnum['ppidnum']), 150, endpoint=True)
fig, axs = plt.subplots(3, 1, sharex=True, figsize=(4.8, 6))
axs[0].hist(ncbi, bins=bins, label='NCBI', color='C0')
axs[1].hist(yo, bins=bins, label=r'Yang $et\ al.$', color='C1')
axs[2].hist(flybase, bins=bins, label='FlyBase', color='C2')
axs[2].set_xlabel('Number of associated polypeptides')
axs[1].set_ylabel('Number of genes')
fig.suptitle('Distribution of genes across\nnumber of associated polypeptides')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.legend()
fig.savefig('out/dist_gnidnum-ppidnum_all.png')
plt.close()

# NCBI and FlyBase
bins = linspace(min(ncbi.min(), flybase.min()), max(ncbi.max(), flybase.max()), 35, endpoint=True)
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].hist(ncbi, bins=bins, label='NCBI', color='C0')
axs[1].hist(flybase, bins=bins, label='FlyBase', color='C2')
axs[1].set_xlabel('Number of associated polypeptides')
fig.suptitle('Distribution of genes across\nnumber of associated polypeptides')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.set_ylabel('Number of genes')
    ax.legend()
fig.savefig('out/dist_gnidnum-ppidnum_NCBI-FB.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (yo, r'Yang $et\ al.$', 'C1', 'YO'),
                                            (flybase, 'FlyBase', 'C2', 'FB')]:
    plt.hist(data, bins=35, label=data_label, color=color)
    plt.xlabel('Number of associated polypeptides')
    plt.ylabel('Number of genes')
    plt.title('Distribution of genes across\nnumber of associated polypeptides')
    plt.legend()
    plt.savefig(f'out/dist_ppidnum-pplen_{file_label}.png')
    plt.close()

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.32_FB2020_01/fasta/dmel-all-translation-r6.32.fasta
../extract_orfs/extract_orfs.py
    ../extract_orfs/out/dpse_translated_orfs.faa
    ../extract_orfs/out/dyak_translated_orfs.faa
./params.tsv
"""