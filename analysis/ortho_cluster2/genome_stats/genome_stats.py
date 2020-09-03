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
with open('params.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, tcds_path = line.split()
        params[spid] = (source, tcds_path)

# Parse polypeptides
sqid0 = 0
sqids = {}
rows = []
for spid, (source, tcds_path) in params.items():
    with open(tcds_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid = re.search(pp_regex[source], line).group(1)
                gnid = re.search(gn_regex[source], line).group(1)
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq0 = ''.join(seqlines)

            try:
                for sqid, seq in sqids[gnid]:
                    if seq0 == seq:
                        rows.append({'ppid': ppid, 'gnid': gnid, 'spid': spid, 'sqid': sqid, 'sqlen': len(seq)})
                        break
                else:
                    sqids[gnid].append((str(sqid0).zfill(8), seq0))
                    rows.append({'ppid': ppid, 'gnid': gnid, 'spid': spid, 'sqid': str(sqid0).zfill(8), 'sqlen': len(seq0)})
                    sqid0 += 1
            except KeyError:
                sqids[gnid] = [(str(sqid0).zfill(8), seq0)]
                rows.append({'ppid': ppid, 'gnid': gnid, 'spid': spid, 'sqid': str(sqid0).zfill(8), 'sqlen': len(seq0)})
                sqid0 += 1

# Make plots output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Define dataframe and groups
df = pd.DataFrame(rows)
ncbi_spids = [spid for spid in params if params[spid][0] == 'NCBI']
flybase_spids = [spid for spid in params if params[spid][0] == 'FlyBase']

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
spid_xsqidnum = spid_ppidnum - spid_sqidnum
ncbi = {spid: spid_xsqidnum[spid] for spid in ncbi_spids}
flybase = {spid: spid_xsqidnum[spid] for spid in flybase_spids}

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
plt.savefig('out/bar_xsqidnum-species.png')
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
ncbi = gnid_nums.loc[gnid_nums['spid'].isin(ncbi_spids), 'ppidnum']
flybase = gnid_nums.loc[gnid_nums['spid'].isin(flybase_spids), 'ppidnum']

# NCBI and FlyBase
bins = linspace(min(ncbi.min(), flybase.min()), max(ncbi.max(), flybase.max()), 35, endpoint=True)
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].hist(ncbi, bins=bins, label='NCBI', color='C0')
axs[1].hist(flybase, bins=bins, label='FlyBase', color='C1')
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
    plt.hist(data, bins=35, label=data_label, color=color)
    plt.xlabel('Number of associated polypeptides')
    plt.ylabel('Number of genes')
    plt.title('Distribution of genes across\nnumber of associated polypeptides')
    plt.legend()
    plt.savefig(f'out/hist_gnidnum-ppidnum_{file_label}.png')
    plt.close()

# 3.2 Distribution of genes across number of associated unique sequences
ncbi = gnid_nums.loc[gnid_nums['spid'].isin(ncbi_spids), 'sqidnum']
flybase = gnid_nums.loc[gnid_nums['spid'].isin(flybase_spids), 'sqidnum']

# NCBI and FlyBase
bins = linspace(min(ncbi.min(), flybase.min()), max(ncbi.max(), flybase.max()), 35, endpoint=True)
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].hist(ncbi, bins=bins, label='NCBI', color='C0')
axs[1].hist(flybase, bins=bins, label='FlyBase', color='C1')
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
    plt.hist(data, bins=35, label=data_label, color=color)
    plt.xlabel('Number of associated unique sequences')
    plt.ylabel('Number of genes')
    plt.title('Distribution of genes across\nnumber of associated unique sequences')
    plt.legend()
    plt.savefig(f'out/hist_gnidnum-sqidnum_{file_label}.png')
    plt.close()

# 3.3 Distribution of genes across number of associated redundant sequences
gnid_nums['xsqidnum'] = gnid_nums['ppidnum'] - gnid_nums['sqidnum']
ncbi = gnid_nums.loc[gnid_nums['spid'].isin(ncbi_spids), 'xsqidnum']
flybase = gnid_nums.loc[gnid_nums['spid'].isin(flybase_spids), 'xsqidnum']

# NCBI and FlyBase
bins = linspace(min(ncbi.min(), flybase.min()), max(ncbi.max(), flybase.max()), 35, endpoint=True)
fig, axs = plt.subplots(2, 1, sharex=True, figsize=(4.8, 6))
axs[0].hist(ncbi, bins=bins, label='NCBI', color='C0')
axs[1].hist(flybase, bins=bins, label='FlyBase', color='C1')
axs[1].set_xlabel('Number of redundant sequences')
fig.suptitle('Distribution of genes across\nnumber of associated redundant sequences')
fig.subplots_adjust(left=0.175)
for ax in axs:
    ax.set_ylabel('Number of genes')
    ax.legend()
fig.savefig('out/hist_gnidnum-xsqidnum.png')
plt.close()

# Individual
for data, data_label, color, file_label in [(ncbi, 'NCBI', 'C0', 'NCBI'),
                                            (flybase, 'FlyBase', 'C1', 'FB')]:
    plt.hist(data, bins=35, label=data_label, color=color)
    plt.xlabel('Number of associated redundant sequences')
    plt.ylabel('Number of genes')
    plt.title('Distribution of genes across\nnumber of associated redundant sequences')
    plt.legend()
    plt.savefig(f'out/hist_gnidnum-xsqidnum_{file_label}.png')
    plt.close()

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
./params.tsv
"""