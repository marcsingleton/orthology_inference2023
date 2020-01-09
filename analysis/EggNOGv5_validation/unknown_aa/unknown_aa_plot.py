"""Plot distribution of counts and percentages of unknown amino acids."""

import matplotlib.pyplot as plt
import os
import pandas as pd

for dir in os.listdir('out/'):
    os.chdir('out/' + dir)
    ali_df = pd.read_pickle('ali_df.pkl')
    seq_df = pd.read_pickle('seq_df.pkl')

    # Overall counts
    print(dir.upper())
    print('number of alignments:', len(ali_df))
    print('number of sequences:', len(seq_df))
    print()

    # Statistics of numbers of unknown amino acids
    ali_df0 = ali_df[ali_df['numseq_x'] > 0]
    seq_df0 = seq_df[seq_df['num_x'] > 0]

    print('number of alignments with unknowns amino acids:', len(ali_df0))
    print('fraction of sequences with unknown amino acids:', len(ali_df0) / len(ali_df))
    print()
    print('number of sequences with unknown amino acids:', len(seq_df0))
    print('fraction of sequences with unknown amino acids:', len(seq_df0) / len(seq_df))
    print()
    print('summary statistics of distribution of numbers of unknown amino acids')
    print(seq_df0['num_x'].describe())
    print()
    print('summary statistics of distribution of fractions of unknown amino acids')
    print(seq_df0['frac_x'].describe())
    print()

    # Distributions of number of unknown amino acids
    numx_hist = seq_df0['num_x']
    plt.figure()
    plt.hist(numx_hist, 50)
    plt.title('Distribution of Numbers of \nUnknown Amino Acids (Excluding 0)')
    plt.xlabel('Number of Unknown Amino Acids')
    plt.ylabel('Number of Sequences')
    plt.savefig('numx.png')

    plt.figure()
    numx_hist_filter = [n for n in numx_hist if n < 200]
    plt.hist(numx_hist_filter, 25)
    plt.title('Distribution of Numbers of \nUnknown Amino Acids Less than 200 (Excluding 0)')
    plt.xlabel('Number of Unknown Amino Acids')
    plt.ylabel('Number of Sequences')
    plt.savefig('numx_filter.png')

    # Distribution of fraction of unknown amino acids
    fracx_hist = seq_df0['frac_x']
    plt.figure()
    plt.hist(fracx_hist, 50)
    plt.title('Distribution of Fractions of \nUnknown Amino Acids (Excluding 0)')
    plt.xlabel('Fraction of Unknown Amino Acids')
    plt.ylabel('Number of Sequences')
    plt.savefig('fracx.png')

    # Distribution of numbers of contiguous regions of unknown amino acids
    numlensx_hist = {}
    for lens in seq_df0['lens_x']:
        numlensx_hist[len(lens)] = numlensx_hist.get(len(lens), 0) + 1
    x, h = zip(*numlensx_hist.items())
    plt.figure()
    plt.bar(x, h, width=1, align='edge')
    plt.title('Distribution of Numbers of \nContiguous Regions of Unknown Amino Acids (Excluding 0)')
    plt.xlabel('Number of Regions')
    plt.ylabel('Number of Sequences')
    plt.savefig('numlensx.png')

    params_corr = {'s': 20, 'edgecolors': 'none', 'alpha': '0.25'}
    # Correlation of number of contiguous regions and number of unknown amino acids
    x = [len(lenx) for lenx in seq_df0['lens_x']]
    y = list(seq_df0['num_x'])
    plt.figure()
    plt.scatter(x, y, **params_corr)
    plt.title('Correlation of Number of \nContiguous Regions and Number of Unknown Amino Acids')
    plt.xlabel('Number of Regions')
    plt.ylabel('Number of Unknown Amino Acids')
    plt.savefig('corr_numlensx_numx.png')

    # Correlation of number of contiguous regions and fraction of unknown amino acids
    x = [len(lenx) for lenx in seq_df0['lens_x']]
    y = list(seq_df0['frac_x'])
    plt.figure()
    plt.scatter(x, y, **params_corr)
    plt.title('Correlation of Number of \nContiguous Regions and Fraction of Unknown Amino Acids')
    plt.xlabel('Number of Regions')
    plt.ylabel('Fraction of Unknown Amino Acids')
    plt.savefig('corr_numlensx_fracx.png')

    # Correlation of maximum contiguous region length and number of unknown amino acids
    x = [max(lenx) for lenx in seq_df0['lens_x']]
    y = list(seq_df0['num_x'])
    plt.figure()
    plt.scatter(x, y, **params_corr)
    plt.title('Correlation of Maximum Contiguous Region\n Length and Number of Unknown Amino Acids')
    plt.xlabel('Maximum Contiguous Region Length')
    plt.ylabel('Number of Unknown Amino Acids')
    plt.savefig('corr_maxlensx_numx.png')

    # Correlation of maximum contiguous region fractional length and fraction of unknown amino acids
    x = [max(lenx) / seq_len for lenx, seq_len in zip(seq_df0['lens_x'], seq_df0['seq_len'])]
    y = list(seq_df0['frac_x'])
    plt.figure()
    plt.scatter(x, y, **params_corr)
    plt.title('Correlation of Maximum Contiguous Region\n Fractional Length and Fraction of Unknown Amino Acids')
    plt.xlabel('Maximum Contiguous Region Fractional Length')
    plt.ylabel('Fraction of Unknown Amino Acids')
    plt.savefig('corr_maxlensx_fracx.png')

    # Correlation of mean contiguous region length and number of unknown amino acids
    x = [sum(lenx) / len(lenx) for lenx in seq_df0['lens_x']]
    y = list(seq_df0['num_x'])
    plt.figure()
    plt.scatter(x, y, **params_corr)
    plt.title('Correlation of Mean Contiguous Region\n Fractional Length and Fraction of Unknown Amino Acids')
    plt.xlabel('Mean Contiguous Region Fractional Length')
    plt.ylabel('Number of Unknown Amino Acids')
    plt.savefig('corr_meanlensx_numx.png')

    # Correlation of mean contiguous region length and fraction of unknown amino acids
    x = [sum(lenx) / (seq_len * len(lenx)) for lenx, seq_len in zip(seq_df0['lens_x'], seq_df0['seq_len'])]
    y = list(seq_df0['frac_x'])
    plt.figure()
    plt.scatter(x, y, **params_corr)
    plt.title('Correlation of Mean Contiguous Region\n Fractional Length and Fraction of Unknown Amino Acids')
    plt.xlabel('Mean Contiguous Region Fractional Length')
    plt.ylabel('Fraction of Unknown Amino Acids')
    plt.savefig('corr_meanlensx_fracx.png')

    # Distribution of lengths of contiguous regions of unknown amino acids
    lensx_hist = []
    for lens in seq_df['lens_x']:
        lensx_hist.extend(lens)
    plt.figure()
    plt.hist(lensx_hist, 50)
    plt.title('Distribution of Lengths of \nContiguous Regions of Unknown Amino Acids')
    plt.xlabel('Length of Region')
    plt.ylabel('Count')
    plt.savefig('lensx.png')

    # Distribution of numbers of sequences with unknown amino acids in alignments
    numseqx_hist = {}
    for numseqx in ali_df0['numseq_x']:
        numseqx_hist[numseqx] = numseqx_hist.get(numseqx, 0) + 1
    x, h = zip(*numseqx_hist.items())
    plt.figure()
    plt.bar(x, h, width=0.25)
    plt.title('Distribution of Numbers of Sequences with \nUnknown Amino Acids in Alignments (Excluding 0)')
    plt.xlabel('Number of Sequences')
    plt.ylabel('Number of Alignments')
    plt.savefig('numseqx.png')

    plt.close('all')
    os.chdir('../..')

"""
OUTPUT
7214_MEMBERS
number of alignments: 13686
number of sequences: 121579

number of alignments with unknowns amino acids: 446
fraction of sequences with unknown amino acids: 0.03258804617857665

number of sequences with unknown amino acids: 467
fraction of sequences with unknown amino acids: 0.0038411238783013515

summary statistics of distribution of numbers of unknown amino acids
count    467.000000
mean      73.404711
std      119.019497
min        1.000000
25%        4.000000
50%       26.000000
75%       98.000000
max      869.000000
Name: num_x, dtype: float64

summary statistics of distribution of fractions of unknown amino acids
count    467.000000
mean       0.099130
std        0.127664
min        0.000297
25%        0.008365
50%        0.039326
75%        0.138619
max        0.492099
Name: frac_x, dtype: float64

EQUAL_+5_MEMBERS
number of alignments: 7800
number of sequences: 72429

number of alignments with unknowns amino acids: 185
fraction of sequences with unknown amino acids: 0.023717948717948717

number of sequences with unknown amino acids: 193
fraction of sequences with unknown amino acids: 0.002664678512750418

summary statistics of distribution of numbers of unknown amino acids
count    193.000000
mean      67.777202
std      112.619645
min        1.000000
25%        4.000000
50%       28.000000
75%       85.000000
max      869.000000
Name: num_x, dtype: float64

summary statistics of distribution of fractions of unknown amino acids
count    193.000000
mean       0.092991
std        0.121126
min        0.000297
25%        0.006993
50%        0.038309
75%        0.129363
max        0.461871
Name: frac_x, dtype: float64

DEPENDENCIES
./unknown_aa_calc.py
    ./out/*/*.pkl
"""