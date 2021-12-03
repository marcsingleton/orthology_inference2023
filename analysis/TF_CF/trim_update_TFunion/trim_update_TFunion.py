"""Merge updated trimmed lists of transcription factors."""

import pandas as pd

df1 = pd.read_csv('../trim_update/nature15545-s1_trim_update_TF.csv')
df2 = pd.read_csv('../trim_update/nmeth.1763-S2_trim_update.csv')

df1.drop(columns=['TYPE', 'FBTR', 'CLUSTER'], inplace=True)
df2.drop(columns=['Name', 'Splice variant', 'CDS length'], inplace=True)
df2.rename(columns={'FBgn': 'FBGN', 'Symbol': 'SYM'}, inplace=True)

df_union = pd.concat([df1, df2], ignore_index=True, sort=True).drop_duplicates()
df_union.to_csv('TF_union.csv', index=False)

"""
DEPENDENCIES
../trim_update/trim_update.py
    ../trim_update/nature15545-s1_trim_update_TF.csv
    ../trim_update/nmeth.1763-S2_trim_update.csv
"""
