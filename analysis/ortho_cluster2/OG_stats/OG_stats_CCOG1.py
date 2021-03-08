"""Plot various statistics of OGs relating to counts of CCs and OGs."""

import matplotlib.pyplot as plt
import os
import pandas as pd

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, spid, _ = line.split()
        ppid2meta[ppid] = gnid, spid

# Load CCs
rows = []
with open('../connect_pgraph/out/pconnect1.txt') as file:
    for line in file:
        CCid, nodes = line.rstrip().split(':')
        for ppid in nodes.split(','):
            rows.append({'CCid': CCid, 'ppid': ppid})
CCs = pd.DataFrame(rows)

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

# Make output directory
if not os.path.exists('out/pgraph1/CCOG/'):
    os.makedirs('out/pgraph1/CCOG/')  # Recursive folder creation

# Plots
# Distribution of polypeptides across number of associated OGs
groups = OGs.groupby('ppid')
dist1 = groups['OGid'].nunique().value_counts()
plt.bar(dist1.index, dist1.values, width=1)
plt.xlim((dist1.index.min()-0.75, dist1.index.max()+0.75))
plt.xlabel('Number of associated OGs')
plt.ylabel('Number of polypeptides')
plt.title('Distribution of polypeptides across\nnumber of associated OGs')
plt.savefig('out/pgraph1/CCOG/hist_ppnum-OGnum_1.png')
plt.close()

dist2 = dist1.drop(1)
plt.bar(dist2.index, dist2.values, width=1)
plt.xlim((dist1.index.min()-0.75, dist1.index.max()+0.75))
plt.xlabel('Number of associated OGs')
plt.ylabel('Number of polypeptides')
plt.title('Distribution of polypeptides across\nnumber of associated OGs')
plt.savefig('out/pgraph1/CCOG/hist_ppnum-OGnum_2.png')
plt.close()

# Distribution of genes across number of associated OGs
groups = OGs.groupby('gnid')
dist1 = groups['OGid'].nunique().value_counts()
plt.bar(dist1.index, dist1.values, width=1)
plt.xlim((dist1.index.min()-0.75, dist1.index.max()+0.75))
plt.xlabel('Number of associated OGs')
plt.ylabel('Number of genes')
plt.title('Distribution of genes across\nnumber of associated OGs')
plt.savefig('out/pgraph1/CCOG/hist_gnnum-OGnum_1.png')
plt.close()

dist2 = dist1.drop(1)
plt.bar(dist2.index, dist2.values, width=1)
plt.xlim((dist1.index.min()-0.75, dist1.index.max()+0.75))
plt.xlabel('Number of associated OGs')
plt.ylabel('Number of genes')
plt.title('Distribution of genes across\nnumber of associated OGs')
plt.savefig('out/pgraph1/CCOG/hist_gnnum-OGnum_2.png')
plt.close()

# Distribution of connected components across number of associated OGs
sizes = OGs.groupby('CCid')['OGid'].nunique()
dist = sizes.reindex(CCs['CCid'].unique(), fill_value=0).value_counts()
plt.bar(dist.index, dist.values, width=1)
plt.xlabel('Number of OGs in connected component')
plt.ylabel('Number of connected components')
plt.title('Distribution of connected components\nacross number of associated OGs')
plt.savefig('out/pgraph1/CCOG/hist_connectnum-OGnum.png')
plt.close()

CCOG_pairs = OGs[['CCid', 'OGid']].drop_duplicates()
OG_OGppnum = OGs.groupby('OGid').size().rename('OG_OGppnum')
CC_OGppnum = CCOG_pairs.join(OG_OGppnum, on='OGid').groupby('CCid')['OG_OGppnum']
CC_OGnum = OGs.groupby('CCid')['OGid'].nunique().rename('CC_OGnum')
CC_CCppnum = CCs.groupby('CCid').size()[OGs['CCid'].unique()].rename('CC_CCppnum')  # Filter out CCs w/o OGs

# Correlation of number of OGs associated with CC and number of polypeptides in CC
plt.scatter(CC_OGnum, CC_CCppnum, alpha=0.5, s=10, edgecolors='none')
plt.xlabel('Number of OGs associated with CC')
plt.ylabel('Number of polypeptides in CC')
plt.savefig('out/pgraph1/CCOG/scatter_CCppnum-CCOGnum.png')
plt.close()

# Correlation of aggregate number of polypeptides in OGs associated with CC with number of OGs associated with CC
plt.scatter(CC_OGppnum.max(), CC_OGnum, alpha=0.5, s=12, edgecolors='none')
plt.xlabel('Max number of polypeptides in OGs associated with CC')
plt.ylabel('Number of OGs associated with CC')
plt.savefig('out/pgraph1/CCOG/scatter_CCOGnum-CCOGppmax.png')
plt.close()

plt.scatter(CC_OGppnum.mean(), CC_OGnum, alpha=0.5, s=10, edgecolors='none')
plt.xlabel('Mean number of polypeptides in OGs associated with CC')
plt.ylabel('Number of OGs associated with CC')
plt.savefig('out/pgraph1/CCOG/scatter_CCOGnum-CCOGppmean.png')
plt.close()

# Correlation of aggregate number of polypeptides in OGs associated with CC with number of polypeptides in CC
plt.scatter(CC_OGppnum.max(), CC_CCppnum, alpha=0.5, s=10, edgecolors='none')
plt.xlabel('Max number of polypeptides in OGs associated with CC')
plt.ylabel('Number of polypeptides in CC')
plt.savefig('out/pgraph1/CCOG/scatter_CCppnum-CCOGppmax.png')
plt.close()

plt.scatter(CC_OGppnum.mean(), CC_CCppnum, alpha=0.5, s=10, edgecolors='none')
plt.xlabel('Mean number of polypeptides in OGs associated with CC')
plt.ylabel('Number of polypeptides in CC')
plt.savefig('out/pgraph1/CCOG/scatter_CCppnum-CCOGppmean.png')
plt.close()

CCpp_pairs = OGs[['CCid', 'ppid']].drop_duplicates()
pp_OGnum = OGs.groupby('ppid')['OGid'].nunique().rename('pp_OGnum')
hspnum = pd.read_csv('../hsp_stats/out/hsps_reciprocal/sppids.tsv', sep='\t', header=0, names=['ppid', 'hspnum', 'gnid', 'spid'], index_col=0)
df = CCpp_pairs.merge(pp_OGnum, on='ppid').merge(CC_CCppnum, on='CCid').merge(CC_OGnum, on='CCid').merge(hspnum, on='ppid')

# Correlation of number of associated OGs with number of polypeptides in CC for polypeptide
plt.scatter(df['pp_OGnum'], df['CC_CCppnum'], alpha=0.5, s=10, edgecolors='none')
plt.xlabel('Number of OGs associated with polypeptide')
plt.ylabel('Number of polypeptides in CC associated with polypeptide')
plt.savefig('out/pgraph1/CCOG/scatter_CCppnum-ppOGnum.png')
plt.close()

# Correlation of number of OGs associated with polypeptide with number of OGs in CC associated with polypeptide
plt.scatter(df['pp_OGnum'], df['CC_OGnum'], alpha=0.5, s=10, edgecolors='none')
plt.xlabel('Number of OGs associated with polypeptide')
plt.ylabel('Number of OGs in CC associated with polypeptide')
plt.savefig('out/pgraph1/CCOG/scatter_CCOGnum-ppOGnum.png')
plt.close()

# Correlation of number of unique reciprocal hits to polypeptide with number of OGs associated with polypeptide
plt.scatter(df['pp_OGnum'], df['hspnum'], alpha=0.5, s=10, edgecolors='none')
plt.xlabel('Number of OGs associated with polypeptide')
plt.ylabel('Number of reciprocal hits to polypeptide')
plt.savefig('out/pgraph1/CCOG/scatter_hitnum-ppOGnum.png')
plt.close()

"""
DEPENDENCIES
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../connect_pgraph/connect_pgraph1.py
    ../connect_pgraph/out/pconnect1.txt
../hsp_stats/hsp_stats.py
    ../hsp_stats/out/hits_reciprocal/sppids.tsv
../subcluster_pgraph/subcluster_pgraph1.py
    ../subcluster_pgraph/out/pgraph1/pclusters.txt
"""