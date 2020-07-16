"""Plot statistics of open reading frames."""

import matplotlib.pyplot as plt
import os
import pandas as pd

for path in os.listdir('out/'):
    max_orfs = pd.read_csv('out/' + path, sep='\t')

    plt.hist(max_orfs['len0'] / max_orfs['len'], bins=50)
    plt.title(f'{path[:4]}:\nFraction of Transcript of Longest ORF')
    plt.xlabel('Fraction of Transcript')
    plt.ylabel('Count')
    plt.savefig(f'out/{path[:4]}_maxhist.png')
    plt.close()

    plt.scatter(max_orfs['len0'], max_orfs['len1'], s=3, alpha=0.25, edgecolors='none')
    plt.title(f'{path[:4]}:\nLengths of Two Longest ORFs')
    plt.xlabel('Longest ORF (nt)')
    plt.ylabel('Second Longest ORF (nt)')
    plt.savefig(f'out/{path[:4]}_lenscatter.png')
    plt.close()

    plt.scatter(max_orfs['len0'] / max_orfs['len'], max_orfs['len1'] / max_orfs['len'], s=3, alpha=0.25, edgecolors='none')
    plt.title(f'{path[:4]}:\nFractions of Transcript of Two Longest ORFs')
    plt.xlabel('Longest ORF')
    plt.ylabel('Second Longest ORF')
    plt.savefig(f'out/{path[:4]}_fracscatter.png')
    plt.close()

"""
NOTES
dgri has an extremely high fraction of transcripts where the longest ORF is nearly the complete transcript
    dgri gene models were validated with long-read sequencing, but it is unclear if or why this is related
In all plots, there is a loose "cluster" between 0.4 and 1 along the x-axis and 0 to 0.2 along the y-axis
    This region likely corresponds to the cases where the longest ORF is likely the best choice
    However the edges are poorly defined, so this is not a robust metric
I will proceed to identify ORFs by homology
    Translate the entire transcript in all three reading frames, then BLAST to nearest species with annotated CDSs

DEPENDENCIES
./orf_stats_calc.py
    ./out/*_maxorfs.tsv
"""