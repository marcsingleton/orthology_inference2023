"""Filter gene association file."""

import os
from functools import reduce
from operator import add

import matplotlib.pyplot as plt
import pandas as pd


def get_ancestors(GO, GOid):
    ancestors = set()
    for parent in GO[GOid]['parents']:
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

# Read ontology
GO = {}
with open('../../../data/GO/go-basic.obo') as file:
    line = file.readline()
    while line:
        if line.startswith('[Term]'):
            parents = []
            alt_ids = []
            is_obsolete = False
            line = file.readline()
            while line and not line.startswith('['):
                if line.startswith('id:'):
                    GOid = line[4:-1]
                if line.startswith('name:'):
                    name = line[5:-1]
                if line.startswith('alt_id:'):
                    alt_ids.append(line[8:-1])
                if line.startswith('is_a:'):
                    parents.append(line.split('!')[0][6:-1])
                if line.startswith('is_obsolete'):
                    is_obsolete = True if line[13:-1] == 'true' else False
                line = file.readline()
            GO[GOid] = {'name': name, 'primary_id': GOid, 'alt_ids': alt_ids, 'parents': parents, 'is_obsolete': is_obsolete}
            for alt_id in alt_ids:
                GO[alt_id] = {'name': name, 'primary_id': GOid, 'alt_ids': [], 'parents': [], 'is_obsolete': is_obsolete}
        else:  # Else suite is necessary so new [Term] line isn't skipped exiting if suite
            line = file.readline()

# Find term ancestors
rows = []
for GOid in GO:
    rows.append({'GOid': GOid, 'ancestor': GOid})  # Include self as ancestor so merge keeps original term
    for ancestor in get_ancestors(GO, GOid):
        rows.append({'GOid': GOid, 'ancestor': ancestor})
ancestors = pd.DataFrame(rows)

# Read raw table and add term names
df1 = pd.read_table('../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/precomputed_files/gene_association_v2.1.fb',
                    skiprows=5,
                    usecols=list(range(15)),  # File contains two spare tabs at end
                    names=['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO ID',  # Official column labels
                           'DB:Reference', 'Evidence', 'With (or) From', 'Aspect', 'DB_Object_Name',
                           'DB_Object_Synonym', 'DB_Object_Type', 'taxon', 'Date', 'Assigned_by'])
df1['name'] = df1['GO ID'].apply(lambda x: GO[x]['name'])

# Drop unneeded columns and filter
mapper = {'DB_Object_ID': 'gnid', 'DB_Object_Symbol': 'symbol', 'Qualifier': 'qualifier', 'GO ID': 'GOid',
          'Evidence': 'evidence', 'Aspect': 'aspect', 'Date': 'date', 'Assigned_by': 'origin'}
df2 = df1[['DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO ID', 'Evidence', 'Aspect', 'taxon', 'Date', 'Assigned_by', 'name']].rename(columns=mapper)

bool1 = df2['qualifier'].isna()  # Remove qualifiers from terms
bool2 = df2['evidence'].isin(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP',  # Remove lower quality annotations
                              'HTP', 'HDA', 'HMP', 'HGI', 'HEP',
                              'TAS', 'IC'])
bool3 = df2['taxon'] == 'taxon:7227'  # Keep only dmel annotations
bool4 = ~df2['GOid'].apply(lambda x: GO[x]['is_obsolete'])  # Remove obsolete annotations
df3 = df2[bool1 & bool2 & bool3 & bool4].drop(['qualifier', 'taxon'], axis=1)

counts = df2.loc[~bool4, ['GOid', 'name']].value_counts()
string = counts.head(10).to_string()
length = len(string.split('\n')[1])
total = str(counts.sum())
print('TOP 10 OBSOLETE ANNOTATIONS (ORIGINAL)')
print(counts.head(10).to_string())
print()
print('TOTAL', (length-6-len(total))*' ' + total)
print()

counts = df2.loc[bool1 & bool2 & bool3 & ~bool4, ['GOid', 'name']].value_counts()
string = counts.head(10).to_string()
length = len(string.split('\n')[1])
total = str(counts.sum())
print('TOP 10 OBSOLETE ANNOTATIONS (FILTERED)')
print(counts.head(10).to_string())
print()
print('TOTAL', (length-6-len(total))*' ' + total)
print()

# Rename ids
df4 = df3.copy()
df4['GOid'] = df4['GOid'].apply(lambda x: GO[x]['primary_id'])

counts = df3.loc[df3['GOid'] != df4['GOid'], ['GOid', 'name']].value_counts()
string = counts.head(10).to_string()
length = len(string.split('\n')[1])
total = str(counts.sum())
print('TOP 10 RENAMED ANNOTATIONS')
print(counts.head(10).to_string())
print()
print('TOTAL', (length-6-len(total))*' ' + total)
print()

# Intersect with OGs
df5 = OGs.merge(df4, on='gnid')

# Propagate ancestors to table and drop poorly represented annotations
df6 = df5.merge(ancestors, on='GOid').drop('GOid', axis=1).rename(columns={'ancestor': 'GOid'}).drop_duplicates()
df7 = df6.groupby('GOid').filter(lambda x: len(x) >= 30)

# Make plots
dfs = [df2, df3, df4, df5, df6, df7]
labels = ['original', 'filter', 'merge', 'intersect', 'propagate', 'drop']

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

# Write dfs to file
for df, label in zip(dfs[2:], labels[2:]):
    df.to_csv(f'out/GAF_{label}.tsv', sep='\t', index=False)

"""
OUTPUT
TOP 10 OBSOLETE ANNOTATIONS (ORIGINAL)
GOid        name                                                                  
GO:0055114   obsolete oxidation-reduction process                                     416
GO:0016458   obsolete gene silencing                                                   39
GO:0005671   obsolete Ada2/Gcn5/Ada3 transcription activator complex                   31
GO:0070868   obsolete heterochromatin organization involved in chromatin silencing     28
GO:0031936   obsolete negative regulation of chromatin silencing                       15
GO:0031937   obsolete positive regulation of chromatin silencing                       13
GO:0031935   obsolete regulation of chromatin silencing                                12
GO:0060968   obsolete regulation of gene silencing                                      9
GO:0000187   obsolete activation of MAPK activity                                       8
GO:0072321   obsolete chaperone-mediated protein transport                              6

TOTAL                                                                                 629

TOP 10 OBSOLETE ANNOTATIONS (FILTERED)
GOid        name                                                                  
GO:0016458   obsolete gene silencing                                                  35
GO:0070868   obsolete heterochromatin organization involved in chromatin silencing    28
GO:0005671   obsolete Ada2/Gcn5/Ada3 transcription activator complex                  27
GO:0031937   obsolete positive regulation of chromatin silencing                      13
GO:0031935   obsolete regulation of chromatin silencing                               12
GO:0031936   obsolete negative regulation of chromatin silencing                      11
GO:0060968   obsolete regulation of gene silencing                                     9
GO:0000185   obsolete activation of MAPKKK activity                                    4
GO:0120081   obsolete positive regulation of microfilament motor activity              3
GO:0000186   obsolete activation of MAPKK activity                                     3

TOTAL                                                                                173

TOP 10 RENAMED ANNOTATIONS
GOid        name                                                                        
GO:0006342   heterochromatin assembly                                                       55
GO:0070491   DNA-binding transcription factor binding                                       21
GO:0048096   epigenetic maintenance of chromatin in transcription-competent conformation    21
GO:0048747   muscle cell development                                                        11
GO:0033613   DNA-binding transcription factor binding                                       10
GO:0030898   microfilament motor activity                                                    8
GO:0001085   RNA polymerase II-specific DNA-binding transcription factor binding             8
GO:0043044   chromatin remodeling                                                            7
GO:0070615   ATP-dependent chromatin remodeler activity                                      6
GO:0008274   gamma-tubulin large complex                                                     5

TOTAL                                                                                      193

DEPENDENCIES
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/precomputed_files/gene_association_v2.1.fb
../../../data/GO/go-basic.obo
../../ortho_cluster3/clique4+_pcommunity/clique4+_pcommunity2.py
    ../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt
../../ortho_MSA/OG_filter/OG_filter.py
    ../../ortho_MSA/OG_filter/out/OG_filter.tsv
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
"""