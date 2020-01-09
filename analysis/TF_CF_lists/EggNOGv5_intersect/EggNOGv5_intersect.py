"""Create lists of TF and CF alignments from the filtered and re-aligned list."""

import os.path

paths_ali = ['../../../data/EggNOGv5/drosophilidae/7214_members.tsv',
             '../../EggNOGv5_validation/filter_unknown_realign/out/7214_noX_members.tsv',
             '../../EggNOGv5_validation/filter_count/out/7214_noX_members/equal_+5_members.tsv']
path_gn = '../gn2ali_reduce/gn2ali_reduce.csv'
path_TF = '../trim_update_TFunion/TF_union.csv'
path_CF = '../trim_update/nature15545-s1_trim_update_CF.csv'

# Dictionaries relating alignment ID to line in member files (for all alignments)
ali_dicts = []
for path_ali in paths_ali:
    ali_dict = {}
    with open(path_ali) as file_ali:
        for line in file_ali:
            fields = line.rstrip().split('\t')
            ali_dict[fields[1]] = line
    ali_dicts.append(ali_dict)

# Dictionary relating gn to alignment ID (for genes with alignments)
gn_dict = {}
with open(path_gn) as file_gn:
    next(file_gn)  # Skip header
    for line in file_gn:
        fields = line.rstrip().split(',')
        gn_dict[fields[0]] = fields[3]

# Write alignments to file for all and filtered member lists
for path_ali, ali_dict in zip(paths_ali, ali_dicts):
    root_out = os.path.splitext(os.path.basename(path_ali))[0]
    with open(root_out + '_TF.tsv', 'w') as member_TF:
        num_TF_all = 0
        num_TF_filt = 0
        with open(path_TF) as file_TF:
            next(file_TF)  # Skip header
            for line in file_TF:
                num_TF_all += 1
                fields = line.rstrip().split(',')
                gn = fields[1]

                # Retrieve member line for alignment ID
                try:
                    ali = gn_dict[gn]
                    ln_all = ali_dict[ali]
                    member_TF.write(ln_all)
                    num_TF_filt += 1
                except KeyError:
                    pass

    print(root_out.upper())
    print('number of TFs:', num_TF_all)
    print('number of TFs with alignments:', num_TF_filt)
    print()

# Write alignments to file for all and filtered member lists
for path_ali, ali_dict in zip(paths_ali, ali_dicts):
    root_out = os.path.splitext(os.path.basename(path_ali))[0]
    with open(root_out + '_CF.tsv', 'w') as member_CF:
        num_CF_all = 0
        num_CF_filt = 0
        with open(path_CF) as file_CF:
            next(file_CF)  # Skip header
            for line in file_CF:
                num_CF_all += 1
                fields = line.rstrip().split(',')
                gn = fields[2]

                # Retrieve member line for alignment ID
                try:
                    ali = gn_dict[gn]
                    ln_all = ali_dict[ali]
                    member_CF.write(ln_all)
                    num_CF_filt += 1
                except KeyError:
                    pass

    print(root_out.upper())
    print('number of CFs:', num_CF_all)
    print('number of CFs with alignments:', num_CF_filt)
    print()

"""
OUTPUT
7214_MEMBERS
number of TFs: 749
number of TFs with alignments: 681

7214_NOX_MEMBERS
number of TFs: 749
number of TFs with alignments: 681

EQUAL_+5_MEMBERS
number of TFs: 749
number of TFs with alignments: 484

7214_MEMBERS
number of CFs: 338
number of CFs with alignments: 316

7214_NOX_MEMBERS
number of CFs: 338
number of CFs with alignments: 316

EQUAL_+5_MEMBERS
number of CFs: 338
number of CFs with alignments: 197

DEPENDENCIES
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214_members.tsv
../../EggNOGv5_validation/filter_count/filter_count.py
    ../../EggNOGv5_validation/filter_count/out/7214_noX_members/equal_+5_members.tsv
../../EggNOGv5_validation/filter_unknown_realign/filter_unknown_realign.py
    ../../EggNOGv5_validation/filter_unknown_realign/out/7214_noX_members.tsv
../gn2ali_reduce/gn2ali_reduce.py
    ../gn2ali_reduce/gn2ali_reduce.csv
../trim_update_TFunion/trim_update_TFunion.py
    ../trim_update_TFunion/TF_union.csv
../trim_update/trim_update.py
    ../trim_update/nature15545-s1_trim_update_CF.csv
"""