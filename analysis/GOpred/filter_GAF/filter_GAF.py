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


def write_table(counts, title):
    if counts.empty:  # Immediately return on empty table
        return
    if os.path.exists('out/output.txt'):
        mode = 'a'
        padding = '\n'
    else:
        padding = ''

    table_string = counts.head(10).to_string()
    length = len(table_string.split('\n')[1])
    total_string = str(counts.sum()).rjust(length - 6, ' ')
    output = f"""\
    {title}
    {table_string}

    TOTAL {total_string}
    """
    with open('out/output.txt', mode) as file:
        file.write(padding + output)


# Load sequence data
ppid2gnid = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppid2gnid[fields['ppid']] = fields['gnid']

# Load regions
rows = []
with open('../../brownian/aucpred_filter/out/regions_30.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        for ppid in fields['ppids'].split(','):
            rows.append({'OGid': fields['OGid'], 'start': int(fields['start']), 'stop': int(fields['stop']),
                         'disorder': fields['disorder'] == 'True', 'gnid': ppid2gnid[ppid]})
regions = pd.DataFrame(rows)

# Load ontology
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
    rows.append({'GOid': GOid, 'ancestor_id': GOid, 'ancestor_name': GO[GOid]['name']})  # Include self as ancestor so merge keeps original term
    for ancestor_id in get_ancestors(GO, GOid):
        rows.append({'GOid': GOid, 'ancestor_id': ancestor_id, 'ancestor_name': GO[ancestor_id]['name']})
ancestors = pd.DataFrame(rows)

# Load raw table and add term names
df1 = pd.read_table('../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.45_FB2022_02/precomputed_files/gene_association.fb',
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

bool1 = df2['qualifier'].isin(['enables', 'contributes_to', 'involved_in',  # Select appropriate qualifiers (FB defines additional qualifiers)
                               'located_in', 'part_of', 'is_active_in'])
bool2 = df2['evidence'].isin(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP',  # Remove lower quality annotations
                              'HTP', 'HDA', 'HMP', 'HGI', 'HEP',
                              'TAS', 'IC'])
bool3 = df2['taxon'] == 'taxon:7227'  # Keep only dmel annotations
bool4 = ~df2['GOid'].apply(lambda x: GO[x]['is_obsolete'])  # Remove obsolete annotations
df3 = df2[bool1 & bool2 & bool3 & bool4].drop(['qualifier', 'taxon'], axis=1)

if not os.path.exists('out/'):
    os.mkdir('out/')

# Identify IDs of obsolete and renamed annotations
counts = df2.loc[~bool4, ['GOid', 'name']].value_counts()
write_table(counts, 'TOP 10 OBSOLETE ANNOTATIONS (ORIGINAL)')

counts = df2.loc[bool1 & bool2 & bool3 & ~bool4, ['GOid', 'name']].value_counts()
write_table(counts, 'TOP 10 OBSOLETE ANNOTATIONS (FILTERED)')

# Update IDs
df4 = df3.copy()
df4['GOid'] = df4['GOid'].apply(lambda x: GO[x]['primary_id'])

counts = df3.loc[df3['GOid'] != df4['GOid'], ['GOid', 'name']].value_counts()
write_table(counts, 'TOP 10 RENAMED ANNOTATIONS')

# Join with regions
df5 = regions.merge(df4, on='gnid')

# Propagate ancestors to table and drop poorly represented annotations
df6 = df5.merge(ancestors, on='GOid').drop(['GOid', 'name'], axis=1).rename(columns={'ancestor_id': 'GOid', 'ancestor_name': 'name'}).drop_duplicates()
df7 = df6.groupby(['GOid', 'disorder']).filter(lambda x: x['gnid'].nunique() >= 50)

# Make plots
dfs = [df2, df3, df4, df5, df6, df7]
labels = ['original', 'filter', 'update', 'join', 'propagate', 'drop']

# Number of annotations
plt.bar(range(len(dfs)), [len(df) for df in dfs], width=0.5, tick_label=labels)
plt.xlabel('Cleaning step')
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
plt.xlabel('Cleaning step')
plt.ylabel('Number of annotations')
plt.savefig('out/bar_numannot-df_aspect.png')
plt.close()

# Number of annotations by evidence code
counts = [df['evidence'].value_counts() for df in dfs]
codes = reduce(lambda x, y: x.combine(y, add, fill_value=0), counts).sort_values(ascending=False)
top_codes = list(codes.index[:9])
other_codes = list(codes.index[9:])
merged_counts = []
for count in counts:
    merged_count = {code: count.get(code, 0) for code in top_codes}
    merged_count['other'] = sum([count.get(code, 0) for code in other_codes])
    merged_counts.append(merged_count)

counts = merged_counts
bottoms = [0 for _ in counts]
for code in (top_codes + ['other']):
    plt.bar(range(len(counts)), [count[code] for count in counts], bottom=bottoms, label=code, width=0.5, tick_label=labels)
    bottoms = [b + count[code] for b, count in zip(bottoms, counts)]
plt.legend()
plt.xlabel('Cleaning step')
plt.ylabel('Number of annotations')
plt.savefig('out/bar_numannot-df_evidence.png')
plt.close()

# Number of terms
plt.bar(range(len(dfs)), [df['GOid'].nunique() for df in dfs], width=0.5, tick_label=labels)
plt.xlabel('Cleaning step')
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
plt.xlabel('Cleaning step')
plt.ylabel('Number of unique terms')
plt.savefig('out/bar_numterms-df_aspect.png')
plt.close()

# Write dfs to file
for df, label in zip(dfs[2:], labels[2:]):
    df.to_csv(f'out/GAF_{label}.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.45_FB2022_02/precomputed_files/gene_association.fb
../../../data/GO/go-basic.obo
../../brownian/aucpred_filter/aucpred_filter.py
    ../../brownian/aucpred_filter/out/regions_30.tsv
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
"""