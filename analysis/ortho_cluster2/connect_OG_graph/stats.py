"""Plot statistics associated with gene groups."""

import matplotlib.pyplot as plt
import pandas as pd

# Load sequence data
ppid2gnid = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppid2gnid[fields['ppid']] = fields['gnid']

# Load OGs
OGs = {}
with open('../cluster4+_graph/out/4clique/clusters.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppids = {node for edge in fields['edges'].split(',') for node in edge.split(':')}
        OGs[fields['OGid']] = ppids

# Load GGs
rows = []
ppid2GGids, gnid2GGids = {}, {}
with open('out/components.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        GGid = fields['GGid']
        OGids = fields['OGids'].split(',')
        ppids = {ppid for OGid in OGids for ppid in OGs[OGid]}
        gnids = {ppid2gnid[ppid] for ppid in ppids}

        rows.append({'GGid': GGid, 'OGidnum': len(OGids), 'ppidnum': len(ppids), 'gnidnum': len(gnids)})

        # Mappings of ppid and gnid to GGids
        for ppid in ppids:
            try:
                ppid2GGids[ppid].append(GGid)
            except KeyError:
                ppid2GGids[ppid] = [GGid]
        for gnid in gnids:
            try:
                gnid2GGids[gnid].append(GGid)
            except KeyError:
                gnid2GGids[gnid] = [GGid]
df = pd.DataFrame(rows)

counts = df['OGidnum'].value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of OGs in gene group')
plt.ylabel('Number of gene groups')
plt.savefig('out/bar_GGidnum-OGidnum.png')
plt.yscale('log')
plt.savefig('out/bar_GGidnum-OGidnum_log.png')
plt.close()

counts = df['ppidnum'].value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of proteins in gene group')
plt.ylabel('Number of gene groups')
plt.savefig('out/bar_GGidnum-ppidnum.png')
plt.yscale('log')
plt.savefig('out/bar_GGidnum-ppidnum_log.png')
plt.close()

counts = df['gnidnum'].value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of genes in gene group')
plt.ylabel('Number of gene groups')
plt.savefig('out/bar_GGidnum-gnidnum.png')
plt.yscale('log')
plt.savefig('out/bar_GGidnum-gnidnum_log.png')
plt.close()

counts = {}
for GGids in ppid2GGids.values():
    counts[len(GGids)] = counts.get(len(GGids), 0) + 1
index, values = list(zip(*counts.items()))
plt.bar(index, values, width=1)
plt.xlabel('Number of gene groups associated with protein')
plt.ylabel('Number of proteins')
plt.savefig('out/bar_ppidnum-GGidnum.png')
plt.yscale('log')
plt.savefig('out/bar_ppidnum-GGidnum_log.png')
plt.close()

counts = {}
for GGids in gnid2GGids.values():
    counts[len(GGids)] = counts.get(len(GGids), 0) + 1
index, values = list(zip(*counts.items()))
plt.bar(index, values, width=1)
plt.xlabel('Number of gene groups associated with gene')
plt.ylabel('Number of genes')
plt.savefig('out/bar_gnidnum-GGidnum.png')
plt.yscale('log')
plt.savefig('out/bar_gnidnum-GGidnum_log.png')
plt.close()

"""
DEPENDENCIES
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../cluster4+_graph/cluster.py
    ../cluster4+_graph/out/4clique/clusters.tsv
./connect_OG_graph.py
    ./out/components.tsv
"""