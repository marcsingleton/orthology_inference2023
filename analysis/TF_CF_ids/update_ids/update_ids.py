"""Update IDs from lists of transcription factors and cofactors."""

import os
import pandas as pd

# Load raw ids and report some basic stats
df1 = pd.read_csv('../../../data/TF_CF_ids/nature15545-s1.csv', usecols=['TYPE', 'SYM', 'FBGN', 'CG', 'FBTR'])  # Stampfel et al.
df2 = pd.read_csv('../../../data/TF_CF_ids/nmeth.1763-S2.csv', usecols=['FBgn', 'Symbol', 'CG', 'Name'])  # Hens et al.

print('STAMPFEL RAW IDS')
print('Number of ids:', len(df1))
print('Number of TFs:', len(df1[df1['TYPE'] == 'TF']))
print('Number of CFs:', len(df1[df1['TYPE'] == 'COFACTOR']))
print('Number of duplicates:', len(df1) - df1['FBGN'].nunique())
print()
print('HENS RAW IDS')
print('Number of ids:', len(df2))
print('Number of duplicates:', len(df2) - df2['FBgn'].nunique())
print()

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

# Report some basic stats
print('STAMPFEL UPDATED IDS')
print('Number of ids:', len(TFs1) + len(CFs))
print('Number of TFs:', len(TFs1))
print('Number of CFs:', len(CFs))
print()
print('HENS UPDATED IDS')
print('Number of ids:', len(TFs2))
print()

# Merge TFs
TFs = TFs1 | TFs2
print('MERGED IDS')
print('Total TFs:', len(TFs))
print('Unique Stampfel et al.:', len(TFs1 - TFs2))
print('Unique Hens et al.:', len(TFs2 - TFs1))

# Write updated IDs to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/TFs.txt', 'w') as file:
    for TF in TFs:
        file.write(f'{TF}\n')
with open('out/CFs.txt', 'w') as file:
    for CF in CFs:
        file.write(f'{CF}\n')

"""
OUTPUT
STAMPFEL RAW IDS
Number of ids: 812
Number of TFs: 474
Number of CFs: 338
Number of duplicates: 0

HENS RAW IDS
Number of ids: 755
Number of duplicates: 6

STAMPFEL UPDATED IDS
Number of ids: 812
Number of TFs: 474
Number of CFs: 338

HENS UPDATED IDS
Number of ids: 746

MERGED IDS
Total TFs: 747
Unique Stampfel et al.: 1
Unique Hens et al.: 273

DEPENDENCIES
../../../data/TF_CF_ids/nature15545-s1.csv
../../../data/TF_CF_ids/nature15545-s1_map.tsv
../../../data/TF_CF_ids/nmeth.1763-S2.csv
../../../data/TF_CF_ids/nmeth.1763-S2_map.tsv
"""