"""Plot various statistics of the input genomes."""

import matplotlib.pyplot as plt
import os
import pandas as pd
import re
from numpy import linspace

pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+(\.[0-9]+)?)'}
gn_regex = {'FlyBase': r'parent=(FBgn[0-9]+)',
            'NCBI': r'db_xref=GeneID:([0-9]+)'}

# Parse parameters
params = {}
with open('../blast_dbs/params.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, tcds_path = line.split()
        params[spid] = (source, tcds_path)

# Parse polypeptides
rows = []
for spid, (source, tcds_path) in params.items():
    with open(tcds_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid = re.search(pp_regex[source], line).group(1)
                gnid = re.search(gn_regex[source], line).group(1)
                line = file.readline()

            pplen = 0
            while line and not line.startswith('>'):
                pplen += len(line.rstrip())
                line = file.readline()

            rows.append({'ppid': ppid, 'gnid': gnid, 'spid': spid, 'pplen': pplen})

# Make plots output directory
if not os.path.exists('out/'):
    os.makedirs('out/')  # Recursive folder creation

# Define dataframe and groups
df = pd.DataFrame(rows)
ncbi_spids = [spid for spid in params if params[spid][0] == 'NCBI']
flybase_spids = [spid for spid in params if params[spid][0] == 'FlyBase']

# Distribution of polypeptides across associated species
spid_ppidnum = df['spid'].value_counts()
ncbi = {spid: spid_ppidnum[spid] for spid in ncbi_spids}
flybase = {spid: spid_ppidnum[spid] for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI', color='C0')
plt.bar(range(len(ncbi), len(ncbi) + len(flybase)),
        flybase.values(), label=r'FlyBase', color='C2')
plt.xticks(list(range(len(ncbi) + len(flybase))), list(ncbi.keys()) + list(flybase.keys()), rotation=60)
plt.xlabel('Associated species')
plt.ylabel('Number of polypeptides')
plt.title('Distribution of polypeptides across associated species')
plt.subplots_adjust(bottom=0.15)
plt.legend()
plt.savefig('out/dist_ppidnum-species.png')
plt.close()

# Distribution of number of genes across associated species
spid_gnidnum = df.groupby('spid')['gnid'].nunique()
ncbi = {spid: spid_gnidnum[spid] for spid in ncbi_spids}
flybase = {spid: spid_gnidnum[spid] for spid in flybase_spids}

plt.bar(range(0, len(ncbi)),
        ncbi.values(), label='NCBI', color='C0')
plt.bar(range(len(ncbi), len(ncbi) + len(flybase)),
        flybase.values(), label=r'FlyBase', color='C2')
plt.xticks(list(range(len(ncbi) + len(flybase))), list(ncbi.keys()) + list(flybase.keys()), rotation=60)
plt.xlabel('Associated species')
plt.ylabel('Number of genes')
plt.title('Distribution of genes across associated species')
plt.subplots_adjust(bottom=0.15)
plt.legend()
plt.savefig('out/dist_gnidnum-species.png')
plt.close()

# Distribution of polypeptides across length
ncbi = df.loc[df['spid'].isin(ncbi_spids), 'pplen']
flybase = df.loc[df['spid'].isin(flybase_spids), 'pplen']

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
fig.savefig('out/dist_ppidnum-pplen.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (flybase, 'FlyBase', 'C2', 'FB')]:
    plt.hist(data, bins=100, label=data_label, color=color)
    plt.xlabel('Length')
    plt.ylabel('Number of polypeptides')
    plt.title('Distribution of polypeptides across length')
    plt.legend()
    plt.savefig(f'out/dist_ppidnum-pplen_{file_label}.png')
    plt.close()

# Distribution of genes across number of associated polypeptides
spidgnid_pairs = df[['spid', 'gnid']].drop_duplicates()
gnid_ppidnum = spidgnid_pairs.join(df.groupby('gnid').size().rename('ppidnum'), on='gnid')
gnid_ppidnum.to_csv('out/gnid_ppidnum.tsv', sep='\t', index=False)
ncbi = gnid_ppidnum.loc[gnid_ppidnum['spid'].isin(ncbi_spids), 'ppidnum']
flybase = gnid_ppidnum.loc[gnid_ppidnum['spid'].isin(flybase_spids), 'ppidnum']

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
fig.savefig('out/dist_gnidnum-ppidnum.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (flybase, 'FlyBase', 'C2', 'FB')]:
    plt.hist(data, bins=35, label=data_label, color=color)
    plt.xlabel('Number of associated polypeptides')
    plt.ylabel('Number of genes')
    plt.title('Distribution of genes across\nnumber of associated polypeptides')
    plt.legend()
    plt.savefig(f'out/dist_gnidnum-ppidnum_{file_label}.png')
    plt.close()

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
./params.tsv
"""