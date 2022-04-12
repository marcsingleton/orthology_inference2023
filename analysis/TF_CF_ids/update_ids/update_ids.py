"""Update IDs from lists of transcription factors and cofactors."""

import os
import pandas as pd

# Load raw ids and report some basic stats
df1 = pd.read_csv('../../../data/TF_CF_ids/nature15545-s1.csv', usecols=['TYPE', 'SYM', 'FBGN', 'CG', 'FBTR'])  # Stampfel et al.
df2 = pd.read_csv('../../../data/TF_CF_ids/nmeth.1763-S2.csv', usecols=['FBgn', 'Symbol', 'CG', 'Name'])  # Hens et al.

# Load Stampfel et al. types
# (Not using pandas to_dict since it creates a nested dictionary)
FBgn2type = {}
with open('../../../data/TF_CF_ids/nature15545-s1.csv') as file:
    file.readline()  # Skip header
    for line in file:
        fields = line.rstrip('\n').split(',')
        if fields[2] in FBgn2type:
            raise RuntimeError('Duplicate FBgns detected.')
        FBgn2type[fields[2]] = fields[0]

# Load Stampfel et al. map
TFs1, CFs, FBids = set(), set(), set()
with open('../../../data/TF_CF_ids/nature15545-s1_map.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        fields = line.rstrip('\n').split('\t')
        if len(fields) == 3 or fields[3] == 'True':
            FBgn1, FBgn2 = fields[0], fields[1]
            if FBgn1 in FBids:
                raise RuntimeError(f'Multiply-mapped keys detected: {FBgn1}')
            if FBgn2type[FBgn1] == 'TF':
                TFs1.add(FBgn2)
            else:
                CFs.add(FBgn2)
            FBids.add(FBgn1)

# Load Hens et al. map
TFs2, FBids = set(), set()
with open('../../../data/TF_CF_ids/nmeth.1763-S2_map.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        fields = line.rstrip('\n').split('\t')
        if len(fields) == 3 or fields[3] == 'True':
            FBgn1, FBgn2 = fields[0], fields[1]
            if FBgn1 in FBids:
                raise RuntimeError(f'Multiply-mapped keys detected: {FBgn1}')
            TFs2.add(FBgn2)
            FBids.add(FBgn1)

TFs = TFs1 | TFs2  # Merge TFs

if not os.path.exists('out/'):
    os.mkdir('out/')

# Report some basic stats
output = f"""\
STAMPFEL RAW IDS
Number of ids: {len(df1)}
Number of TFs: {len(df1[df1['TYPE'] == 'TF'])}
Number of CFs: {len(df1[df1['TYPE'] == 'COFACTOR'])}
Number of duplicates: {len(df1) - df1['FBGN'].nunique()}

HENS RAW IDS
Number of ids: {len(df2)}
Number of duplicates: {len(df2) - df2['FBgn'].nunique()}

STAMPFEL UPDATED IDS
Number of ids: {len(df2)}
Number of duplicates: {len(df2) - df2['FBgn'].nunique()}

HENS UPDATED IDS
Number of ids: {len(TFs2)}

MERGED IDS
Total TFs: {len(TFs)}
Unique Stampfel et al.: {len(TFs1 - TFs2)}
Unique Hens et al.: {len(TFs2 - TFs1)}
"""
with open('out/output.txt', 'w') as file:
    file.write(output)

# Write updated IDs to file
with open('out/TFs.txt', 'w') as file:
    for TF in TFs:
        file.write(f'{TF}\n')
with open('out/CFs.txt', 'w') as file:
    for CF in CFs:
        file.write(f'{CF}\n')

"""
DEPENDENCIES
../../../data/TF_CF_ids/nature15545-s1.csv
../../../data/TF_CF_ids/nature15545-s1_map.tsv
../../../data/TF_CF_ids/nmeth.1763-S2.csv
../../../data/TF_CF_ids/nmeth.1763-S2_map.tsv
"""