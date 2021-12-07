"""Filter gene association file."""

import os
from functools import reduce
from operator import add

import matplotlib.pyplot as plt
import pandas as pd


def get_ancestors(GO, GOid):
    ancestors = set()
    for parent in GO[GOid]:
        ancestors.add(parent)
        ancestors.update(get_ancestors(GO, parent))
    return ancestors


# Load seq metadata
ppid2gnid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, _, sqid = line.split()
        ppid2gnid[ppid] = gnid

# Load OGs
rows = []
with open('../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt') as file:
    for line in file:
        CCid, OGid, edges = line.rstrip().split(':')
        ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        for ppid in ppids:
            gnid = ppid2gnid[ppid]
            rows.append({'CCid': CCid, 'OGid': OGid, 'gnid': gnid})
OG_gnids = pd.DataFrame(rows).drop_duplicates()  # All OGs with genes
OG_filter = pd.read_table('../../ortho_MSA/OG_filter/out/OG_filter.tsv', usecols=['CCid', 'OGid', 'gOGid'])  # OGs after filtering
OGs = OG_filter.merge(OG_gnids, how='left', on=['OGid', 'CCid'])  # Filtered OGs with genes

# Read raw table
df1 = pd.read_table('../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/precomputed_files/gene_association_v2.1.fb',
                    skiprows=5,
                    usecols=list(range(15)),  # File contains two spare tabs at end
                    names=['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO ID',  # Official column labels
                           'DB:Reference', 'Evidence', 'With (or) From', 'Aspect', 'DB_Object_Name',
                           'DB_Object_Synonym', 'DB_Object_Type', 'taxon', 'Date', 'Assigned_by'])

# Clean table
mapper = {'DB_Object_ID': 'gnid', 'DB_Object_Symbol': 'symbol', 'Qualifier': 'qualifier', 'GO ID': 'GOid',
          'Evidence': 'evidence', 'Aspect': 'aspect', 'Date': 'date', 'Assigned_by': 'origin'}
df2 = df1[['DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO ID', 'Evidence', 'Aspect', 'taxon', 'Date', 'Assigned_by']].rename(columns=mapper)

bool1 = df2['qualifier'].isna()  # Remove qualifiers from terms
bool2 = df2['evidence'].isin(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP',  # Remove lower quality annotations
                              'HTP', 'HDA', 'HMP', 'HGI', 'HEP',
                              'TAS', 'IC'])
bool3 = df2['taxon'] == 'taxon:7227'  # Keep only dmel annotations
df3 = df2[bool1 & bool2 & bool3].drop(['qualifier', 'taxon'], axis=1)
df4 = OGs.merge(df3, on='gnid')

# Read ontology
GO = {}
with open('../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/precomputed_files/go-basic.obo') as file:
    line = file.readline()
    while line:
        if line.startswith('[Term]'):
            parents = []
            line = file.readline()
            while line and not line.startswith('['):
                if line.startswith('id:'):
                    GOid = line[4:-1]
                if line.startswith('is_a:'):
                    parents.append(line.split('!')[0][6:-1])
                line = file.readline()
            GO[GOid] = parents
        else:  # Else suite is necessary so new [Term] line isn't skipped exiting if suite
            line = file.readline()


# Find term ancestors
rows = []
for GOid in GO:
    rows.append({'GOid': GOid, 'ancestor': GOid})  # Include self as ancestor so merge keeps original term
    for ancestor in get_ancestors(GO, GOid):
        rows.append({'GOid': GOid, 'ancestor': ancestor})
ancestors = pd.DataFrame(rows)

# Propagate ancestors to table and drop poorly represented annotations
df5 = df4.merge(ancestors, on='GOid').drop('GOid', axis=1).rename(columns={'ancestor': 'GOid'}).drop_duplicates()
df6 = df5.groupby('GOid').filter(lambda x: len(x) >= 30)

# Make plots
dfs = [df2, df3, df4, df5, df6]
labels = ['original', 'filtered', 'merged', 'propagated', 'dropped']

if not os.path.exists('out/'):
    os.mkdir('out/')

# Number of annotations
plt.bar(range(len(dfs)), [len(df) for df in dfs], width=0.5, tick_label=labels)
plt.ylabel('Number of annotations')
plt.savefig('out/bar_numannot-df.png')
plt.close()

# Number of annotations by aspect
counts = [df['aspect'].value_counts() for df in dfs]
bottoms = [0 for count in counts]
for aspect, label in [('P', 'Process'), ('F', 'Function'), ('C', 'Component')]:
    plt.bar(range(len(counts)), [count[aspect] for count in counts], bottom=bottoms, label=label, width=0.5, tick_label=labels)
    bottoms = [b + count[aspect] for b, count in zip(bottoms, counts)]
plt.legend()
plt.ylabel('Number of annotations')
plt.savefig('out/bar_numannot-df_aspect.png')
plt.close()

# Number of annotations by evidence code
counts = [df['evidence'].value_counts() for df in dfs]
codes = reduce(lambda x, y: x.combine(y, add, fill_value=0), counts).sort_values(ascending=False).index[:10]
bottoms = [0 for count in counts]
for code in codes:
    plt.bar(range(len(counts)), [count.get(code, 0) for count in counts], bottom=bottoms, label=code, width=0.5, tick_label=labels)
    bottoms = [b + count.get(code, 0) for b, count in zip(bottoms, counts)]
plt.legend()
plt.ylabel('Number of annotations')
plt.savefig('out/bar_numannot-df_evidence.png')
plt.close()

# Number of terms
plt.bar(range(len(dfs)), [df['GOid'].nunique() for df in dfs], width=0.5, tick_label=labels)
plt.ylabel('Number of unique terms')
plt.savefig('out/bar_numterms-df.png')
plt.close()

# Number of terms by aspect
counts = [df[['GOid', 'aspect']].drop_duplicates()['aspect'].value_counts() for df in dfs]
bottoms = [0 for count in counts]
for aspect, label in [('P', 'Process'), ('F', 'Function'), ('C', 'Component')]:
    plt.bar(range(len(counts)), [count[aspect] for count in counts], bottom=bottoms, label=label, width=0.5, tick_label=labels)
    bottoms = [b + count[aspect] for b, count in zip(bottoms, counts)]
plt.legend()
plt.ylabel('Number of unique terms')
plt.savefig('out/bar_numterms-df_aspect.png')
plt.close()

"""
DEPENDENCIES
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/precomputed_files/gene_association_v2.1.fb
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/precomputed_files/go-basic.obo
../../ortho_cluster3/clique4+_pcommunity/clique4+_pcommunity2.py
    ../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt
../../ortho_MSA/OG_filter/OG_filter.py
    ../../ortho_MSA/OG_filter/out/OG_filter.tsv
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
"""