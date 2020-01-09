"""Calculate distribution of alignment percentages."""

import gzip
import json
import os
import scipy.stats as stats
from Bio import AlignIO
from numpy import quantile

pathdir = [('../../../data/EggNOGv5/drosophilidae/7214_members.tsv', '../../../data/EggNOGv5/drosophilidae/7214/'),
           ('../filter_unknown_realign/out/7214_noX_members.tsv', '../filter_unknown_realign/out/align/'),
           ('../filter_count/out/7214_members/equal_+5_members.tsv', '../filter_unknown_realign/out/align/'),
           ('../../TF_CF_lists/EggNOGv5_intersect/equal_+5_members_CF.tsv', '../filter_unknown_realign/out/align/'),
           ('../../TF_CF_lists/EggNOGv5_intersect/7214_noX_members_CF.tsv', '../filter_unknown_realign/out/align/'),
           ('../../TF_CF_lists/EggNOGv5_intersect/equal_+5_members_TF.tsv', '../filter_unknown_realign/out/align/'),
           ('../../TF_CF_lists/EggNOGv5_intersect/7214_noX_members_TF.tsv', '../filter_unknown_realign/out/align/')]

def fraction_ungapped(MSA):
    fractions = []
    length = MSA.get_alignment_length()
    for record in MSA:
        fractions.append(1 - record.seq.count('-') / length)
    return fractions


for path_members, dir_msa in pathdir:
    # Create list of alignments
    members = []
    with open(path_members) as file_member:
        for line in file_member:
            fields = line.rstrip().split('\t')
            members.append(fields[1])

    # Open alignments and calculate statistics
    means = []
    stds = []
    iqrs = []
    for member in members:
        with gzip.open(dir_msa + member + '.raw_alg.faa.gz', 'rt') as file:
            MSA = AlignIO.read(file, 'fasta')
            fractions = fraction_ungapped(MSA)
            means.append(stats.tmean(fractions))
            stds.append(stats.tstd(fractions))
            iqrs.append(quantile(fractions, 0.75) - quantile(fractions, 0.25))

    # Save statistics to folder
    root, _ = os.path.splitext(os.path.basename(path_members))  # Get name of member file without extension
    if not os.path.exists('out/' + root):
        os.makedirs('out/' + root)  # Recursive folder creation

    with open('out/' + root + '/means.json', 'w') as file:
        json.dump(means, file)
    with open('out/' + root + '/stds.json', 'w') as file:
        json.dump(stds, file)
    with open('out/' + root + '/iqrs.json', 'w') as file:
        json.dump(iqrs, file)

"""
DEPENDENCIES
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214/
    ../../../data/EggNOGv5/drosophilidae/7214_members.tsv
../filter_unknown_realign/filter_unknown_realign.py
    ../filter_unknown_realign/out/7214_noX_members.tsv
    ../filter_unknown_realign/out/align/
../filter_count/filter_count.py
    ../filter_count/out/7214_members/equal_+5_members.tsv
../../TF_CF_lists/EggNOGv5_intersect/EggNOGv5_intersect.py
    ../../TF_CF_lists/EggNOGv5_intersect/equal_+5_members_CF.tsv
    ../../TF_CF_lists/EggNOGv5_intersect/7214_noX_members_CF.tsv
    ../../TF_CF_lists/EggNOGv5_intersect/equal_+5_members_TF.tsv
    ../../TF_CF_lists/EggNOGv5_intersect/7214_noX_members_TF.tsv
"""