"""Plot correlations of mean, STD, and IQR as well as the distribution of means for the alignment percentage."""

import json
import matplotlib.pyplot as plt
import os

params = {'s': 5, 'color': 'black', 'edgecolors': 'none', 'alpha': '0.125'}
dirs = filter(os.path.isdir, os.listdir('out/'))

for dir in os.listdir('out/'):
    os.chdir('out/' + dir)

    with open('means.json') as file:
        means = json.load(file)
    with open('stds.json') as file:
        stds = json.load(file)
    with open('iqrs.json') as file:
        iqrs = json.load(file)

    plt.figure()
    plt.scatter(means, stds, **params)
    plt.xlabel('Mean of Fraction Ungapped')
    plt.ylabel('Standard Deviation of Fraction Ungapped')
    plt.savefig('corr_mean_std.png')

    plt.figure()
    plt.scatter(means, iqrs, **params)
    plt.xlabel('Mean of Fraction Ungapped')
    plt.ylabel('IQR of Fraction Ungapped')
    plt.savefig('corr_mean_iqr.png')

    plt.figure()
    plt.scatter(stds, iqrs, **params)
    plt.xlabel('Standard Deviation of Fraction Ungapped')
    plt.ylabel('IQR of Fraction Ungapped')
    plt.savefig('corr_std_iqr.png')

    plt.figure()
    plt.hist(means, bins=50)
    plt.xlabel('Fraction Ungapped')
    plt.ylabel('Count')
    plt.savefig('hist_mean.png')

    plt.close('all')
    os.chdir('../..')

"""
NOTES
Without plotting the scatters on the same graph, it is difficult to determine if the relationships between the means,
STDs, and IQRs significantly differ between the data sets. I should directly compare each data subset to its parent data
set.

The histograms of means more obviously show differences between CFs and TFs. CFs have an exponential distribution with
fewer alignments at lower alignment percentages. TFs have a bell-shaped distribution centered between 0.8 and 0.9.

DEPENDENCIES
./ungapped_percentage_dists_calc.py
    ./out/*/*.json
"""