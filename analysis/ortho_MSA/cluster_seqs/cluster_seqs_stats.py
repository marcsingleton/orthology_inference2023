"""Plot various statistics related to clustering of unique sequences in genes and OGs."""

import os

import matplotlib.pyplot as plt
import pandas as pd

# Load seq metadata
gnid2spid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        _, gnid, spid = line.split()
        gnid2spid[gnid] = spid

# Load OGs
rows = []
with open('../../ortho_cluster3/clique4+_community/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        for gnid in gnids:
            rows.append({'CCid': CCid, 'OGid': OGid, 'gnid': gnid, 'spid': gnid2spid[gnid]})
OGs = pd.DataFrame(rows)


df = pd.read_table('out/clusters.tsv', dtype={'gnid': str})

if not os.path.exists('out/'):
    os.mkdir('out/')

# 1.1 Distribution of genes over number of unique sequences
x = df[['gnid', 'repr_id']].drop_duplicates().groupby('gnid').size().value_counts()
plt.bar(x.index, x.values)
plt.xlabel('Number of unique sequences')
plt.ylabel('Number of genes')
plt.savefig('out/bar_gnnum-reprnum.png')
plt.yscale('log')
plt.savefig('out/bar_gnnum-reprnum_log.png')
plt.close()

# 1.2 Distribution of genes over number of clusters
x = df[['gnid', 'cluster_id']].drop_duplicates().groupby('gnid').size().value_counts()
plt.bar(x.index, x.values)
plt.xlabel('Number of sequence clusters')
plt.ylabel('Number of genes')
plt.savefig('out/bar_gnnum-clusternum.png')
plt.yscale('log')
plt.savefig('out/bar_gnnum-clusternum_log.png')
plt.close()

# 2 Correlation of number of unique sequences with number of clusters
x = df[['gnid', 'cluster_id']].drop_duplicates().groupby('gnid').size()
y = df[['gnid', 'repr_id']].drop_duplicates().groupby('gnid').size()
plt.scatter(x.values, y.values, s=10, alpha=0.5, edgecolors='none')
plt.xlabel('Number of sequence clusters')
plt.ylabel('Number of unique sequences')
plt.savefig('out/scatter_reprnum-clusternum.png')
plt.close()

# Aggregate counts at level of OGs
gnid_idnum = df.drop(columns=['ppid', 'pident']).groupby('gnid').nunique()
OGid_gnnum = OGs.groupby('OGid').size().rename('gnnum')

OGid_idsum = OGs.merge(gnid_idnum, on='gnid', how='left').groupby('OGid')[['repr_id', 'cluster_id']].sum()
sum_stats = OGid_idsum.merge(OGid_gnnum, on='OGid')

OGid_idmax = OGs.merge(gnid_idnum, on='gnid', how='left').groupby('OGid')[['repr_id', 'cluster_id']].max()
max_stats = OGid_idmax.merge(OGid_gnnum, on='OGid')

# 3.1.1 Distribution of OGs over number of unique sequences
x = sum_stats['repr_id']
plt.hist(x, bins=75)
plt.xlabel('Number of unique sequences in OG')
plt.ylabel('Number of OGs')
plt.savefig('out/hist_OGnum-reprnum.png')
plt.yscale('log')
plt.savefig('out/hist_OGnum-reprnum_log.png')
plt.close()

# 3.1.2 Distribution of OGs over number of unique sequences per gene
x = sum_stats['repr_id'] / sum_stats['gnnum']
plt.hist(x, bins=75)
plt.xlabel('Number of unique sequences per gene in OG')
plt.ylabel('Number of OGs')
plt.savefig('out/hist_OGnum-reprnum_gene.png')
plt.close()

# 3.1.3 Distribution of OGs over number of unique sequences (26 genes)
x = sum_stats.loc[sum_stats['gnnum'] == 26, 'repr_id']
plt.hist(x, bins=75)
plt.xlabel('Number of unique sequences in OG')
plt.ylabel('Number of OGs')
plt.savefig('out/hist_OGnum-reprnum_26.png')
plt.yscale('log')
plt.savefig('out/hist_OGnum-reprnum_26_log.png')
plt.close()

# 3.2.1 Distribution of OGs over number of clusters
x = sum_stats['cluster_id']
plt.hist(x, bins=75)
plt.xlabel('Number of sequence clusters in OG')
plt.ylabel('Number of OGs')
plt.savefig('out/hist_OGnum-clusternum.png')
plt.yscale('log')
plt.savefig('out/hist_OGnum-clusternum_log.png')
plt.close()

# 3.2.2 Distribution of OGs over number of clusters per gene
x = sum_stats['cluster_id'] / sum_stats['gnnum']
plt.hist(x, bins=75)
plt.xlabel('Number of sequence clusters per gene in OG')
plt.ylabel('Number of OGs')
plt.savefig('out/hist_OGnum-clusternum_gene.png')
plt.close()

# 3.2.3 Distribution of OGs over number of clusters (26 genes)
x = sum_stats.loc[sum_stats['gnnum'] == 26, 'cluster_id']
plt.hist(x, bins=75)
plt.xlabel('Number of sequence clusters in OG')
plt.ylabel('Number of OGs')
plt.savefig('out/hist_OGnum-clusternum_26.png')
plt.yscale('log')
plt.savefig('out/hist_OGnum-clusternum_26_log.png')
plt.close()

# 4.1 Correlation of number of unique sequences for gene with number of unique sequences in OG
x = max_stats['repr_id'] / sum_stats['repr_id']
plt.hist(x, bins=50)
plt.xlabel('Fraction of total unique sequences in OG associated with single gene')
plt.ylabel('Number of genes')
plt.savefig('out/hist_gnnum-repr_max_fraction.png')
plt.close()

# 4.2 Correlation of number of unique sequences for gene with number of unique sequences in OG (26 genes)
x = max_stats.loc[max_stats['gnnum'] == 26, 'repr_id'] / sum_stats.loc[sum_stats['gnnum'] == 26, 'repr_id']
plt.hist(x, bins=50)
plt.xlabel('Fraction of total unique sequences in OG associated with single gene')
plt.ylabel('Number of genes')
plt.savefig('out/hist_gnnum-repr_max_fraction_26.png')
plt.close()

# 4.3 Correlation of number of clusters for gene with number of clusters in OG
x = max_stats['cluster_id'] / sum_stats['cluster_id']
plt.hist(x, bins=50)
plt.xlabel('Fraction of total sequence clusters in OG associated with single gene')
plt.ylabel('Number of genes')
plt.savefig('out/hist_gnnum-cluster_max_fraction.png')
plt.close()

# 4.4 Correlation of number of clusters for gene with number of clusters in OG (26 genes)
x = max_stats.loc[max_stats['gnnum'] == 26, 'cluster_id'] / sum_stats.loc[sum_stats['gnnum'] == 26, 'cluster_id']
plt.hist(x, bins=50)
plt.xlabel('Fraction of total sequence clusters in OG associated with single gene')
plt.ylabel('Number of genes')
plt.savefig('out/hist_gnnum-cluster_max_fraction_26.png')
plt.close()

"""
DEPENDENCIES
../../ortho_cluster3/clique4+_community/clique4+_community2.py
    ../../ortho_cluster3/clique4+_community/out/ggraph2/5clique/gclusters.txt
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
./cluster_seqs_calc.py
    ./out/clusters.tsv
"""