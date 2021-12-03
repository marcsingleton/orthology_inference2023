"""Update trimmed lists with correct CG and FBgn identifiers."""

import pandas as pd

# Raw CSVs of TFs and CFs
df1_raw = pd.read_csv('../trim/nature15545-s1_trim.csv')  # Stampfel et al. list
df2_raw = pd.read_csv('../trim/nmeth.1763-S2_trim.csv')  # Hens et al. list

# Updated CSVs of TF and CFs from flybase
df1_FB = pd.read_csv('nature15545-s1_FBfields.tsv', sep='\t')
df2_FB = pd.read_csv('nmeth.1763-S2_FBfields.tsv', sep='\t')

# Calculate number of updated CGs and update
# Sorting is unnecessary as the FlyBase tsv was generated with the records in the same order
num_cg1_update = sum(df1_raw['CG'] != df1_FB['ANNOTATION_SYMBOL'])
num_cg2_update = sum(df2_raw['CG'] != df2_FB['ANNOTATION_SYMBOL'])

df1_raw['CG'] = df1_FB['ANNOTATION_SYMBOL']
df2_raw['CG'] = df2_FB['ANNOTATION_SYMBOL']

# Calculate number of updated FBgns and update
num_fbgn1_update = sum(df1_raw['FBGN'] != df1_FB['FBID_KEY'])
num_fbgn2_update = sum(df2_raw['FBgn'] != df2_FB['FBID_KEY'])

df1_raw['FBGN'] = df1_FB['FBID_KEY']
df2_raw['FBgn'] = df2_FB['FBID_KEY']

# Calculate number of updated symbols and update
num_sym1_update = sum(df1_raw['SYM'] != df1_FB['SYMBOL'])
num_sym2_update = sum(df2_raw['Symbol'] != df2_FB['SYMBOL'])

df1_raw['SYM'] = df1_FB['SYMBOL']
df2_raw['Symbol'] = df2_FB['SYMBOL']

# Save updated CSVs
df1_raw.to_csv('nature15545-s1_trim_update.csv', index=False)
df2_raw.to_csv('nmeth.1763-S2_trim_update.csv', index=False)

df1_raw[df1_raw['TYPE'] == 'TF'].to_csv('nature15545-s1_trim_update_TF.csv', index=False)
df1_raw[df1_raw['TYPE'] == 'COFACTOR'].to_csv('nature15545-s1_trim_update_CF.csv', index=False)

# Print analysis
print('DF1_RAW UPDATES')
print('number of CGs updated:', num_cg1_update)
print('number of FBGNs updated:', num_fbgn1_update)
print('number of symbols updated:', num_sym1_update)
print()
print('DF2_RAW UPDATES')
print('number of CGs updated:', num_cg2_update)
print('number of FBGNs updated:', num_fbgn2_update)
print('number of symbols updated:', num_sym2_update)

"""
OUTPUT
DF1_RAW UPDATES
number of CGs updated: 23
number of FBGNs updated: 34
number of symbols updated: 108

DF2_RAW UPDATES
number of CGs updated: 42
number of FBGNs updated: 103
number of symbols updated: 130

DEPENDENCIES
../trim/
    ../trim/nature15545-s1_trim.csv
    ../trim/nmeth.1763-S2_trim.csv
"""