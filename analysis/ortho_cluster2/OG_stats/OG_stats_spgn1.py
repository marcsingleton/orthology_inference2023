"""Plot various statistics of OGs relating to their counts of polypeptides and species."""

import matplotlib.pyplot as plt
import os
import pandas as pd

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, spid, _ = line.split()
        ppid2meta[ppid] = gnid, spid

# Load OGs
rows = []
with open('../subcluster_pgraph/out/pgraph1/pclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        for ppid in ppids:
            gnid, spid = ppid2meta[ppid]
            rows.append({'CCid': CCid, 'OGid': OGid, 'ppid': ppid, 'gnid': gnid, 'spid': spid})
OGs = pd.DataFrame(rows)

spid_num = 35
groups = OGs.groupby('OGid')
OGnum = OGs['OGid'].nunique()
uppnum = OGs['ppid'].nunique()
ugnnum = OGs['gnid'].nunique()
OG_unums = groups.nunique()
unum1 = len(OG_unums[(OG_unums['gnid'] == spid_num) & (OG_unums['spid'] == spid_num)])
unum2 = len(OG_unums[(OG_unums['ppid'] == spid_num) & (OG_unums['spid'] == spid_num)])

# Make output directory
if not os.path.exists('out/pgraph1/spgn/'):
    os.makedirs('out/pgraph1/spgn/')  # Recursive folder creation

# Plots
# Number of associated OGs for species
spid_OGs = {}
for _, group in groups:
    for spid in group['spid'].unique():
        spid_OGs[spid] = spid_OGs.get(spid, 0) + 1

labels, h_OG = zip(*sorted(spid_OGs.items(), key=lambda i: i[0]))
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_OG)
ax1.set_xticks(x)
ax1.set_xticklabels(labels, rotation=60, fontsize=8)
ax1.set_xlabel('Species')
ax1.set_ylabel('Number of associated OGs')
ax1.set_title('Number of associated OGs for each species')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / OGnum, mx / OGnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/pgraph1/spgn/bar_OGnum-species.png')
plt.close()

# Distribution of polypeptides across species
spid_upps = OGs.groupby('spid')['ppid'].nunique().sort_index()
labels, h_upp = zip(*spid_upps.items())
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_upp)
ax1.set_xticks(x)
ax1.set_xticklabels(labels, rotation=60, fontsize=8)
ax1.set_xlabel('Associated species')
ax1.set_ylabel('Number of polypeptides')
ax1.set_title('Distribution of polypeptides across associated species')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / uppnum, mx / uppnum)
ax2.set_ylabel('Fraction of total polypeptides')

fig.tight_layout()
fig.savefig('out/pgraph1/spgn/bar_uppnum-species.png')
plt.close()

# Distribution of genes across species
spid_ugns = OGs.groupby('spid')['gnid'].nunique().sort_index()
labels, h_ugn = zip(*spid_ugns.items())
x = list(range(1, len(labels) + 1))
fig, ax1 = plt.subplots()
ax1.bar(x, h_ugn)
ax1.set_xticks(x)
ax1.set_xticklabels(labels, rotation=60, fontsize=8)
ax1.set_xlabel('Associated species')
ax1.set_ylabel('Number of genes')
ax1.set_title('Distribution of genes across associated species')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / ugnnum, mx / ugnnum)
ax2.set_ylabel('Fraction of total genes')

fig.tight_layout()
fig.savefig('out/pgraph1/spgn/bar_ugnnum-species.png')
plt.close()

# Number of exclusions for each species
for i in range(spid_num-10, spid_num):
    spids = set(OGs['spid'].drop_duplicates())
    spid_counts = {spid: 0 for spid in sorted(spids)}
    for spids in [spids - set(group['spid'].drop_duplicates()) for _, group in groups if group['spid'].nunique() == i]:
        for spid in spids:
            spid_counts[spid] += 1
    labels, h = zip(*spid_counts.items())
    x = list(range(1, len(labels) + 1))
    fig, ax1 = plt.subplots()
    ax1.bar(x, h)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels, rotation=60, fontsize=8)
    ax1.set_xlabel('Species')
    ax1.set_ylabel('Number of OGs')
    ax1.set_title(f'Number of exclusions for each species in OGs with {i} species')

    ax2 = ax1.twinx()
    mn, mx = ax1.get_ylim()
    OG_num = sum(h) / (spid_num - i)
    ax2.set_ylim(mn / OG_num, mx / OG_num)
    ax2.set_ylabel(f'Fraction of total OGs with {i} species')

    fig.tight_layout()
    fig.savefig(f'out/pgraph1/spgn/bar_OGexclusion{i}-species.png')
    plt.close()

# Correlation of number of polypeptides and associated OGs
fig, ax = plt.subplots()
ax.scatter(h_OG, h_upp, s=6)
ax.set_xlabel('Number of associated OGs')
ax.set_ylabel('Number of associated polypeptides')
ax.set_title('Correlation of numbers of associated\n polypeptides and OGs for each species')

fig.savefig('out/pgraph1/spgn/scatter_uppnum-OGnum.png')
plt.close()

# Correlation of number of genes and associated OGs
fig, ax = plt.subplots()
ax.scatter(h_OG, h_ugn, s=6)
ax.set_xlabel('Number of associated OGs')
ax.set_ylabel('Number of associated genes')
ax.set_title('Correlation of numbers of associated\n genes and OGs for each species')

fig.savefig('out/pgraph1/spgn/scatter_ugnnum-OGnum.png')
plt.close()

# Distribution of OGs across number of species
dist_sp = groups['spid'].nunique().value_counts()
sp, sp_count = zip(*dist_sp.items())
fig, ax1 = plt.subplots()
ax1.bar(sp, sp_count)
ax1.set_title('Distribution of OGs across number of species in OG')
ax1.set_xlabel('Number of species in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / OGnum, mx / OGnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/pgraph1/spgn/hist_OGnum-spnum.png')
plt.close()

# Distribution of OGs across number of genes
dist_gn = groups['gnid'].nunique().value_counts()
gn, gn_count = zip(*dist_gn.items())
fig, ax1 = plt.subplots()
ax1.bar(gn, gn_count, width=1)
ax1.set_title('Distribution of OGs across number of genes in OG')
ax1.set_xlabel('Number of genes in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / OGnum, mx / OGnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/pgraph1/spgn/hist_OGnum-gnnum.png')
plt.close()

# Distribution of OGs across number of polypeptides
dist_pp = groups['ppid'].nunique().value_counts()
pp, pp_count = zip(*dist_pp.items())
fig, ax1 = plt.subplots()
ax1.bar(pp, pp_count, width=1)
ax1.set_title('Distribution of OGs across number of polypeptides in OG')
ax1.set_xlabel('Number of polypeptides in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / OGnum, mx / OGnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/pgraph1/spgn/hist_OGnum-ppnum.png')
plt.close()

# Distribution of OGs across number of species duplicates
dist_dup = (groups.size() - groups['spid'].nunique()).value_counts()
seq, seq_count = zip(*dist_dup.drop(0, errors='ignore').items())  # Drop 0; ignore error if 0 does not exist
fig, ax1 = plt.subplots()
ax1.bar(seq, seq_count, width=1)
ax1.set_title('Distribution of OGs across number of species duplicates in OG')
ax1.set_xlabel('Number of species duplicates in OG')
ax1.set_ylabel('Number of OGs')

ax2 = ax1.twinx()
mn, mx = ax1.get_ylim()
ax2.set_ylim(mn / OGnum, mx / OGnum)
ax2.set_ylabel('Fraction of total OGs')

fig.tight_layout()
fig.savefig('out/pgraph1/spgn/hist_OGnum-spdup.png')
plt.close()

# Print counts
print('number of OGs:', OGnum)
print()
print(f'number of OGs with {spid_num} species:', dist_sp[spid_num])
print(f'fraction of OGs with {spid_num} species:', dist_sp[spid_num] / OGnum)
print()
print(f'number of OGs with {spid_num} genes:', dist_gn[spid_num])
print(f'fraction of OGs with {spid_num} genes:', dist_gn[spid_num] / OGnum)
print()
print(f'number of OGs with {spid_num} polypeptides:', dist_pp[spid_num])
print(f'fraction of OGs with {spid_num} polypeptides:', dist_pp[spid_num] / OGnum)
print()
print(f'number of OGs with {spid_num} species and {spid_num} genes:', unum1)
print(f'fraction of OGs with {spid_num} species and {spid_num} genes:', unum1 / OGnum)
print()
print(f'number of OGs with {spid_num} species and {spid_num} polypeptides:', unum2)
print(f'fraction of OGs with {spid_num} species and {spid_num} polypeptides:', unum2 / OGnum)
print()
print('number of OGs with duplicates:', OGnum - dist_dup[0])
print('fraction of OGs with duplicates', (OGnum - dist_dup[0]) / OGnum)

"""
OUTPUT
number of OGs: 23053

number of OGs with 35 species: 9455
fraction of OGs with 35 species: 0.4101418470481065

number of OGs with 35 genes: 7161
fraction of OGs with 35 genes: 0.3106320218626643

number of OGs with 35 polypeptides: 1777
fraction of OGs with 35 polypeptides: 0.07708324296187047

number of OGs with 35 species and 35 genes: 6815
fraction of OGs with 35 species and 35 genes: 0.29562312931071877

number of OGs with 35 species and 35 polypeptides: 1400
fraction of OGs with 35 species and 35 polypeptides: 0.060729623042554116

number of OGs with duplicates: 13920
fraction of OGs with duplicates 0.603825966251681

NOTES
These plots are largely based off those in analysis/EggNOGv5_validation/ali_stats/ali_stats.py

DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../subcluster_pgraph/subcluster_pgraph1.py
    ../subcluster_pgraph/pgraph1/out/pclusters.tsv
"""