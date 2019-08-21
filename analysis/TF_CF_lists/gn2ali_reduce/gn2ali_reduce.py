"""
Filter mappings from Drosophila genes to alignments to only genes which have alignments.
Determine if every gene corresponds to only one alignment.
"""

import pandas as pd

df_all = pd.read_csv('../gn2ali/gn2ali.tsv', sep='\t')
df_reduce = df_all.dropna(subset=['EggNOGv5_ID'])
df_reduce_dedup = df_reduce.drop_duplicates(subset='FlyBase_FBgn')

df_reduce.to_csv('gn2ali_reduce.csv', index=False)

print('length of reduced gene to alignment mapping:', len(df_reduce))
print('length of unique reduced gene to alignment mapping:', len(df_reduce_dedup))

"""
OUTPUT
length of reduced gene to alignment mapping: 11485
length of unique reduced gene to alignment mapping: 11485

NOTES
Each gene is represented by only one alignment.
Because each polypeptide sequence is also represented by only alignment, the relationship between genes and polypeptide
sequences is 1-to-1.

DEPENDENCIES
../gn2ali/gn2ali.py
    ../gn2ali/gn2ali.tsv
"""