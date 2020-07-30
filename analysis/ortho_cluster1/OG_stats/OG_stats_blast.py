"""Plot various statistics of OGs relating to their BLAST parameters."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from math import log10
from numpy import linspace


def hist1(df, bins, file_label, title_label, x_label, df_label, color, capital=True, wrap=False):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel((x_label[0].upper() + x_label[1:]) if capital else x_label)
    plt.ylabel(f'Number of {title_label}')
    plt.title(f'Distribution of {title_label} across' + ('\n' if wrap else ' ') + x_label)
    plt.legend()
    plt.savefig(f'out/blast/hist_{file_label}.png')
    plt.close()


def hist3(dfs, bins, file_label, title_label, x_label, df_labels, colors, capital=True, wrap=False):
    fig, axs = plt.subplots(3, 1, figsize=(4.8, 6), sharex=True)
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for ax, df, data_label, color in zip(axs, dfs, df_labels, colors):
        ax.hist(df, bins=b, label=data_label, color=color)
        ax.legend()
    axs[2].set_xlabel(x_label[0].upper() + x_label[1:] if capital else x_label)
    axs[1].set_ylabel(f'Number of {title_label}')
    fig.suptitle(f'Distribution of {title_label} across' + ('\n' if wrap else ' ') + x_label)
    fig.subplots_adjust(left=0.175)
    fig.savefig(f'out/blast/hist_{file_label}.png')
    plt.close()


def scatter1(x, y, file_label, xy_label, df_label, color):
    fig, ax = plt.subplots()
    ax.scatter(x, y, alpha=0.5, s=10, label=df_label, color=color, edgecolors='none')
    ax.set_xlabel(f'Mean {xy_label} in OG')
    ax.set_ylabel(f'Variance of {xy_label} in OG')
    leg = ax.legend(markerscale=2)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    fig.savefig(f'out/blast/scatter_{file_label}.png')
    plt.close()


def scatter2(x, y, file_label, y_label):
    plt.scatter(x, y, alpha=0.5, s=10, label='all', color='C0', edgecolors='none')
    plt.xlabel(f'Number of genes in OG')
    plt.ylabel(y_label)
    plt.savefig(f'out/blast/scatter_{file_label}.png')
    plt.close()


df0 = pd.read_csv('out/ggraph.tsv', sep='\t', dtype={col: 'string' for col in ['qppid', 'qgnid', 'qspid', 'sppid', 'sgnid', 'sspid', 'OGid']})
df0['fali'] = (df0['qend'] - df0['qstart'] + 1) / df0['qlen']

# Segment OGs
OGs = df0.groupby('OGid')
gns = df0.groupby('qgnid')

OG_spidnum = OGs['qspid'].nunique()
OG_gnidnum = OGs['qgnid'].nunique()
gn_OGidnum = gns['OGid'].nunique()

OGs_10sps = OG_spidnum[OG_spidnum == 10].index
OGs_10gns = OG_gnidnum[OG_gnidnum == 10].index
gns_gn1OG = gn_OGidnum[gn_OGidnum > 1].index
OGs_gn1OG = df0.loc[df0['qgnid'].isin(gns_gn1OG), 'OGid'].unique()

OGs1 = set(OGs_10sps) & set(OGs_10gns)
OGs2 = OGs1 - set(OGs_gn1OG)

df1 = df0[df0['OGid'].isin(OGs1)]
df2 = df0[df0['OGid'].isin(OGs2)]
dfs = [df0, df1, df2]
OGs = [df.groupby('OGid') for df in dfs]

# Make output directory
if not os.path.exists('out/blast/'):
    os.makedirs('out/blast/')  # Recursive folder creation

# Size of groups
xs = ['all', 'filter1', 'filter2']
ys = [df['OGid'].nunique() for df in dfs]
plt.bar(xs, ys, color=['C0', 'C1', 'C2'], width=0.25)
plt.xlim((-0.75, 2.75))
plt.ylabel('Number of OGs')
plt.savefig('out/blast/bar_OGnum-filter.png')
plt.close()

# EVALUE PLOTS
evalues = [df.loc[df0['evalue'] != 0, 'evalue'].apply(log10) for df in dfs]

# Stacked bar
xs = list(range(3))
ys_g0 = [len(evalue) / len(df) for evalue, df in zip(evalues, dfs)]  # Fraction of hits with evalue > 0
ys_e0 = [1 - y_g0 for y_g0 in ys_g0]  # Fraction of hits with evalue == 0
plt.bar(xs, ys_g0, label='non-zero', width=0.25)
plt.bar(xs, ys_e0, label='zero', width=0.25, bottom=ys_g0)
plt.xticks(xs, ['all', 'filter1', 'filter2'])
plt.xlim((-0.75, 2.75))
plt.ylabel('Fraction of total hits in OGs')
plt.title('Fraction of hits in OGs with zero and non-zero E-values')
plt.legend(bbox_to_anchor=(0.5, -0.1875), loc='lower center', ncol=2)
plt.subplots_adjust(bottom=0.15)
plt.savefig('out/blast/bar_evalue_zero.png')
plt.close()

# Histograms
hist3([evalue for evalue in evalues], 200, 'hitnum-evalue', 'hits in OGs', 'log10(E-value)',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'], capital=False)
hist1(evalues[0], 200, 'hitnum-evalue_all', 'hits in OGs', 'log10(E-value)', 'all', 'C0', capital=False)
hist1(evalues[1], 200, 'hitnum-evalue_filter1', 'hits in OGs', 'log10(E-value)', 'filter1', 'C1', capital=False)
hist1(evalues[2], 200, 'hitnum-evalue_filter2', 'hits in OGs', 'log10(E-value)', 'filter2', 'C2', capital=False)

# BITSCORE HISTOGRAMS
hist3([df['bitscore'] for df in dfs], 200, 'hitnum-bitscore', 'hits in OGs', 'bitscore',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'])
hist1(df0['bitscore'], 200, 'hitnum-bitscore_all', 'hits in OGs', 'bitscore', 'all', 'C0')
hist1(df1['bitscore'], 200, 'hitnum-bitscore_filter1', 'hits in OGs', 'bitscore', 'filter1', 'C1')
hist1(df2['bitscore'], 200, 'hitnum-bitscore_filter2', 'hits in OGs', 'bitscore', 'filter2', 'C2')

hist3([OG['bitscore'].mean() for OG in OGs], 200, 'OGnum-bitscoremean', 'OGs', 'mean bitscore of hits in OG',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'], wrap=True)
hist3([OG['bitscore'].var() for OG in OGs], 200, 'OGnum-bitscorevar', 'OGs', 'variance of bitscore of hits in OG',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'], wrap=True)
scatter1(OGs[0]['bitscore'].mean(), OGs[0]['bitscore'].var(),
        'bitscorevar-bitscoremean_all.png', 'bitscore', 'all', 'C0')
scatter1(OGs[1]['bitscore'].mean(), OGs[1]['bitscore'].var(),
        'bitscorevar-bitscoremean_filter1.png', 'bitscore', 'filter1', 'C1')
scatter1(OGs[2]['bitscore'].mean(), OGs[2]['bitscore'].var(),
        'bitscorevar-bitscoremean_filter2.png', 'bitscore', 'filter2', 'C2')

# PIDENT HISTOGRAMS
hist3([df['pident'] for df in dfs], 50, 'hitnum-pident', 'hits in OGs', 'percent identity',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'])
hist1(df0['pident'], 50, 'hitnum-pident_all', 'hits in OGs', 'percent identity', 'all', 'C0')
hist1(df1['pident'], 50, 'hitnum-pident_filter1', 'hits in OGs', 'percent identity', 'filter1', 'C1')
hist1(df2['pident'], 50, 'hitnum-pident_filter2', 'hits in OGs', 'percent identity', 'filter2', 'C2')

hist3([OG['pident'].mean() for OG in OGs], 50, 'OGnum-pidentmean', 'OGs', 'mean percent identity of hits in OG',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'], wrap=True)
hist3([OG['pident'].var() for OG in OGs], 50, 'OGnum-pidentvar', 'OGs', 'variance of percent identity of hits in OG',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'], wrap=True)
scatter1(OGs[0]['pident'].mean(), OGs[0]['pident'].var(),
        'pidentvar-pidentmean_all.png', 'percent identity', 'all', 'C0')
scatter1(OGs[1]['pident'].mean(), OGs[1]['pident'].var(),
        'pidentvar-pidentmean_filter1.png', 'percent identity', 'filter1', 'C1')
scatter1(OGs[2]['pident'].mean(), OGs[2]['pident'].var(),
        'pidentvar-pidentmean_filter2.png', 'percent identity', 'filter2', 'C2')

# FRACTION ALIGNED HISTOGRAMS
hist3([df['fali'] for df in dfs], 50, 'hitnum-fali', 'hits in OGs', 'fraction of query aligned',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'], wrap=True)
hist1(df0['fali'], 50, 'hitnum-fali_all', 'hits in OGs', 'fraction of query aligned', 'all', 'C0', wrap=True)
hist1(df1['fali'], 50, 'hitnum-fali_filter1', 'hits in OGs', 'fraction of query aligned', 'filter1', 'C1', wrap=True)
hist1(df2['fali'], 50, 'hitnum-fali_filter2', 'hits in OGs', 'fraction of query aligned', 'filter2', 'C2', wrap=True)

hist3([OG['fali'].mean() for OG in OGs], 50, 'OGnum-falimean', 'OGs', 'mean fraction of query aligned of hits in OG',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'], wrap=True)
hist3([OG['fali'].var() for OG in OGs], 50, 'OGnum-falivar', 'OGs', 'variance of fraction of query aligned of hits in OG',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'], wrap=True)
scatter1(OGs[0]['fali'].mean(), OGs[0]['fali'].var(),
        'falivar-falimean_all.png', 'fraction aligned', 'all', 'C0')
scatter1(OGs[1]['fali'].mean(), OGs[1]['fali'].var(),
        'falivar-falimean_filter1.png', 'fraction aligned', 'filter1', 'C1')
scatter1(OGs[2]['fali'].mean(), OGs[2]['fali'].var(),
        'falivar-falimean_filter2.png', 'fraction aligned', 'filter2', 'C2')

# EDGES
edgenums = [df[['qgnid', 'sgnid', 'OGid']].drop_duplicates().groupby('OGid').size() / 2 for df in dfs]
gnidnums = [OG['qgnid'].nunique() for OG in OGs]
edgefracs = [2*edgenum / (gnidnum*(gnidnum-1)) for edgenum, gnidnum in zip(edgenums, gnidnums)]

hist3(edgenums, 100, 'OGnum-edgenum', 'OGs', 'number of edges',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'])
hist1(edgenums[0], 100, 'OGnum-edgenum_all', 'OGs', 'number of edges', 'all', 'C0')
hist1(edgenums[1], 50, 'OGnum-edgenum_filter1', 'OGs', 'number of edges', 'filter1', 'C1')
hist1(edgenums[2], 50, 'OGnum-edgenum_filter2', 'OGs', 'number of edges', 'filter2', 'C2')

hist3(edgefracs, 50, 'OGnum-edgefrac', 'OGs', 'fraction of possible edges',
      ['all', 'filter1', 'filter2'], ['C0', 'C1', 'C2'])
hist1(edgefracs[0], 50, 'OGnum-edgefrac_all', 'OGs', 'fraction of possible edges', 'all', 'C0')
hist1(edgefracs[1], 50, 'OGnum-edgefrac_filter1', 'OGs', 'fraction of possible edges', 'filter1', 'C1')
hist1(edgefracs[2], 50, 'OGnum-edgefrac_filter2', 'OGs', 'fraction of possible edges', 'filter2', 'C2')

# CORRELATIONS
scatter2(gnidnums[0], OGs[0]['bitscore'].mean(),
         f'bitscore-OGgnnum_all', 'Mean bitscore of hits in OG')
scatter2(gnidnums[0], OGs[0]['pident'].mean(),
         f'pident-OGgnnum_all', 'Mean percent identity of hits in OG')
scatter2(gnidnums[0], OGs[0]['fali'].mean(),
         f'fali-OGgnnum_all', 'Mean fraction of query aligned of hits in OG')
scatter2(gnidnums[0], edgenums[0],
         f'edgenum-OGgnnum_all', 'Number of edges in OG')
scatter2(gnidnums[0], edgefracs[0],
         f'edgefrac-OGgnnum_all', 'Fraction of possible edges in OG')

# Remove His OGs
df3 = df0[~df0['OGid'].isin(['0122', '0044', '0059', '004b', '011e'])]
OG3 = df3.groupby('OGid')
edgenum = df3[['qgnid', 'sgnid', 'OGid']].drop_duplicates().groupby('OGid').size() / 2
gnidnum = OG3['qgnid'].nunique()
edgefrac = 2*edgenum / (gnidnum*(gnidnum-1))

scatter2(gnidnum, OG3['bitscore'].mean(),
         f'bitscore-OGgnnum_his', 'Mean bitscore of hits in OG')
scatter2(gnidnum, OG3['pident'].mean(),
         f'pident-OGgnnum_his', 'Mean percent identity of hits in OG')
scatter2(gnidnum, OG3['fali'].mean(),
         f'fali-OGgnnum_his', 'Mean fraction of query aligned of hits in OG')
scatter2(gnidnum, edgenum,
         f'edgenum-OGgnnum_his', 'Number of edges in OG')
scatter2(gnidnum, edgefrac,
         f'edgefrac-OGgnnum_his', 'Fraction of possible edges in OG')

"""
DEPENDENCIES
./OG_stats_extract.py
    ./out/ggraph.tsv
"""