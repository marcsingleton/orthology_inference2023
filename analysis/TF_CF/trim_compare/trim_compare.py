"""Compare the two lists of TFs to determine the amount of overlap."""

import pandas as pd

TF1_df = pd.read_csv('../trim/nature15545-s1_trim_TF.csv')  # Stampfel et al. list
TF2_df = pd.read_csv('../trim/nmeth.1763-S2_trim.csv')  # Hens et al. list

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
size of set 2: 749
size of intersection: 459
size of union: 764
number of unique elements in TF1: 15
number of unique elements in TF2: 290

FBGN
size of set 1: 474
size of set 2: 749
size of intersection: 440
size of union: 783
number of unique elements in TF1: 34
number of unique elements in TF2: 309

size of unique elements in CG intersection: 23
size of unique elements in FBGN intersection: 4

NOTES
Examining some unique elements in the CG intersection showed cases where the CGs were the same between records, but the
FBgns differed. Looking into these cases on FlyBase showed genes can have multiple CG and FBgn identifiers, though one
of each is considered the primary. FlyBase was then queried for the updated information via the batch download tool
(flybase.org/batchdownload) using the CGs as input. See ../trim_update for further information.

Set 2 has 755 records, meaning 6 records are duplicates. All set 1 records are unique.

DEPENDENCIES
../trim/
    ../trim/nature15545-s1_trim_TF.csv
    ../trim/nmeth.1763-S2_trim.csv
"""