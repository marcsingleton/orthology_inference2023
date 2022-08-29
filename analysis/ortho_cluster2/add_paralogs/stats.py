"""Statistics related to the addition of in-paralogs to clusters."""

import matplotlib.pyplot as plt
import pandas as pd

# Load OGs
OGs1 = {}
with open('../cluster4+_graph/out/4clique/clusters.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppids = {node for edge in fields['edges'].split(',') for node in edge.split(':')}
        OGs1[fields['OGid']] = ppids

OGs2 = {}
with open('out/clusters.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppids = {node for edge in fields['edges'].split(',') for node in edge.split(':')}
        OGs2[fields['OGid']] = ppids

# Calculate stats
rows = []
for OGid, OG1 in OGs1.items():
    OG2 = OGs2[OGid]
    rows.append({'OGid': OGid, 'ppidnum1': len(OG1), 'ppidnum2': len(OG2), 'delta': len(OG2 - OG1)})
df = pd.DataFrame(rows)

# Make plots
plt.scatter(df['ppidnum1'], df['ppidnum2'], s=5, alpha=0.2)
plt.xlabel('Number of proteins in OG w/o in-paralogs')
plt.ylabel('Number of proteins in OG w/ in-paralogs')
plt.savefig('out/scatter_ppidnum2-ppidnum1.png')
plt.close()

counts = (df['delta'] == 0).value_counts()
labels = [('w/o in-paralogs' if idx else 'w/ in-paralogs') + f'\n({value:,})' for idx, value in zip(counts.index, counts.values)]
plt.pie(counts.values, labels=labels, labeldistance=1.4, textprops={'ha': 'center'})
plt.title('OGs w/ and w/o in-paralogs')
plt.savefig('out/pie_paralog.png')
plt.close()

counts = df.loc[df['delta'] != 0, 'delta'].value_counts()
plt.bar(counts.index, counts.values, width=1)
plt.xlabel('Number of in-paralogs')
plt.ylabel('Number of OGs')
plt.savefig('out/bar_OGnum-paralognum.png')
plt.yscale('log')
plt.savefig('out/bar_OGnum-paralognum_log.png')
plt.close()

"""
DEPENDENCIES
../cluster4+_graph/cluster4+_graph.py
    ../cluster4+_graph/out/4clique/clusters.tsv
./add_paralogs.py
    ../out/clusters.tsv
"""