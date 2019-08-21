"""Compare the two updated lists of TFs to determine the amount of overlap."""

import pandas as pd

TF1_df = pd.read_csv('../trim_update/nature15545-s1_trim_update_TF.csv')  # Stampfel et al. list
TF2_df = pd.read_csv('../trim_update/nmeth.1763-S2_trim_update.csv')  # Hens et al. list

TF1_set_cg = set(TF1_df['CG'])
TF2_set_cg = set(TF2_df['CG'])
TF_int_cg = TF1_set_cg.intersection(TF2_set_cg)
TF_union_cg = TF1_set_cg.union(TF2_set_cg)

TF1_set_fb = set(TF1_df['FBGN'])
TF2_set_fb = set(TF2_df['FBgn'])
TF_int_fb = TF1_set_fb.intersection(TF2_set_fb)
TF_union_fb = TF1_set_fb.union(TF2_set_fb)

fb2cg = set(TF1_df.loc[TF1_df['FBGN'].isin(TF_int_fb)]['CG'])  # Convert list of FBGNs to CGs
cg_uniint = TF_int_cg.difference(fb2cg)
fb_uniint = fb2cg.difference(TF_int_cg)

with open('CG_intersection_unique.txt', 'w') as file:
    for cg in cg_uniint:
        file.write(cg + '\n')
with open('FBGN_intersection_unique.txt', 'w') as file:
    for fb in fb_uniint:
        file.write(fb + '\n')

print('CG')
print('size of set 1:', len(TF1_set_cg))
print('size of set 2:', len(TF2_set_cg))
print('size of intersection:', len(TF_int_cg))
print('size of union:', len(TF_union_cg))
print('number of unique elements in TF1:', len(TF1_set_cg) - len(TF_int_cg))
print('number of unique elements in TF2:', len(TF2_set_cg) - len(TF_int_cg))
print()
print('FBGN')
print('size of set 1:', len(TF1_set_fb))
print('size of set 2:', len(TF2_set_fb))
print('size of intersection:', len(TF_int_fb))
print('size of union:', len(TF_union_fb))
print('number of unique elements in TF1:', len(TF1_set_fb) - len(TF_int_fb))
print('number of unique elements in TF2:', len(TF2_set_fb) - len(TF_int_fb))
print()
print('size of unique elements in CG intersection:', len(cg_uniint))
print('size of unique elements in FBGN intersection:', len(fb_uniint))

"""
OUTPUT
CG
size of set 1: 474
size of set 2: 747
size of intersection: 472
size of union: 749
number of unique elements in TF1: 2
number of unique elements in TF2: 275

FBGN
size of set 1: 474
size of set 2: 747
size of intersection: 472
size of union: 749
number of unique elements in TF1: 2
number of unique elements in TF2: 275

size of unique elements in CG intersection: 0
size of unique elements in FBGN intersection: 0

NOTES
Set 2 is has 755 records, meaning 8 records are duplicates. This is an increase of 2 duplicates from the original list.
All set 1 records are still unique.

DEPENDENCIES
../trim_update/
    ../trim_update/nature15545-s1_trim_update_TF.csv
    ../trim_update/nmeth.1763-S2_trim_update.csv
"""