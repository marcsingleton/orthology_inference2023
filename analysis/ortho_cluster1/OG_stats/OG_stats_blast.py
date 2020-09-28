"""Plot various statistics of OGs relating to their BLAST parameters."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from math import log10
from numpy import linspace


def hist1(df, bins, file_label, title_label, x_label, df_label, color, capital=True, wrap=False):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel(x_label[0].upper() + x_label[1:] if capital else x_label)
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
    fig.savefig(f'out/blast/hist3_{file_label}.png')
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
    plt.xlabel('Number of genes in OG')
    plt.ylabel(y_label)
    plt.savefig(f'out/blast/scatter_{file_label}.png')
    plt.close()


hsp_dtypes = {'qppid': 'string', 'qgnid': 'string', 'qspid': 'string',
              'sppid': 'string', 'sgnid': 'string', 'sspid': 'string',
              'length': int, 'nident': int,
              'qlen': int, 'qstart': int, 'qend': int,
              'evalue': float, 'bitscore': float}
id_dtypes = {'qgnid': 'string', 'sgnid': 'string',
             'CCid': 'string', 'OGid': 'string'}

# Load data
df = pd.read_csv('../blast2hsps/out/hsps.tsv', sep='\t', usecols=hsp_dtypes.keys(), dtype=hsp_dtypes)
r = pd.read_csv('../hsps2reciprocal/out/hsps.tsv', sep='\t', usecols=['reciprocal'], memory_map=True)
ids = pd.read_csv('../OGs2hsps/out/hsps.tsv', sep='\t', usecols=id_dtypes.keys(), dtype=id_dtypes)

dfr = df.join(r)
hsps0 = dfr[dfr['reciprocal']].merge(ids, how='left', on=['qgnid', 'sgnid']).dropna()
hsps0['pident'] = hsps0['nident'] / hsps0['length']
hsps0['nqa'] = hsps0['qend'] - hsps0['qstart'] + 1
hsps0['fqa'] = hsps0['nqa'] / hsps0['qlen']

# Segment OGs
OGs = hsps0.groupby('OGid')
gns = hsps0.groupby('qgnid')

OG_spidnum = OGs['qspid'].nunique()
OG_gnidnum = OGs['qgnid'].nunique()
gn_OGidnum = gns['OGid'].nunique()

OGs_10sps = OG_spidnum[OG_spidnum == 10].index
OGs_10gns = OG_gnidnum[OG_gnidnum == 10].index
gns_gn1OG = gn_OGidnum[gn_OGidnum > 1].index
OGs_gn1OG = hsps0.loc[hsps0['qgnid'].isin(gns_gn1OG), 'OGid'].unique()

OGs1 = set(OGs_10sps) & set(OGs_10gns)
OGs2 = OGs1 - set(OGs_gn1OG)

hsps1 = hsps0[hsps0['OGid'].isin(OGs1)]
hsps2 = hsps0[hsps0['OGid'].isin(OGs2)]
hsps = [hsps0, hsps1, hsps2]
OGs = [hsp.groupby('OGid') for hsp in hsps]

labels = ['all', 'filter1', 'filter2']
colors = ['C0', 'C1', 'C2']

# Make output directory
if not os.path.exists('out/blast/'):
    os.makedirs('out/blast/')  # Recursive folder creation

# 1 FILTER PLOT
ys = [hsp['OGid'].nunique() for hsp in hsps]
plt.bar(labels, ys, color=colors, width=0.25)
plt.xlim((-0.75, 2.75))
plt.ylabel('Number of OGs')
plt.savefig('out/blast/bar_OGnum-filter.png')
plt.close()

# 2 EVALUE PLOTS
evalues = [hsp.loc[hsps0['evalue'] != 0, 'evalue'].apply(log10) for hsp in hsps]

# 2.1 Stacked bars of zero and non-zero E-values in all and reciprocal data sets
xs = list(range(3))
ys_g0 = [len(evalue) / len(hsp) for evalue, hsp in zip(evalues, hsps)]  # Fraction of HSPs with evalue > 0
ys_e0 = [1 - y_g0 for y_g0 in ys_g0]  # Fraction of HSPs with evalue == 0
plt.bar(xs, ys_g0, label='non-zero', width=0.25)
plt.bar(xs, ys_e0, label='zero', width=0.25, bottom=ys_g0)
plt.xticks(xs, labels)
plt.xlim((-0.75, 2.75))
plt.ylabel('Fraction of total HSPs in OGs')
plt.title('Fraction of HSPs in OGs with zero and non-zero E-values')
plt.legend(bbox_to_anchor=(0.5, -0.1875), loc='lower center', ncol=2)
plt.subplots_adjust(bottom=0.15)
plt.savefig('out/blast/bar_evalue_zero.png')
plt.close()

# 2.2 Histograms of non-zero E-values
hist3([evalue for evalue in evalues], 200, 'hspnum-evalue', 'HSPs in OGs', 'log10(E-value)', labels, colors, capital=False)
hist1(evalues[0], 200, 'hspnum-evalue_all', 'HSPs in OGs', 'log10(E-value)', labels[0], colors[0], capital=False)
hist1(evalues[1], 200, 'hspnum-evalue_filter1', 'HSPs in OGs', 'log10(E-value)', labels[1], colors[1], capital=False)
hist1(evalues[2], 200, 'hspnum-evalue_filter2', 'HSPs in OGs', 'log10(E-value)', labels[2], colors[2], capital=False)

# 3 BITSCORE PLOTS
# 3.1 Hit histograms
hist3([hsp['bitscore'] for hsp in hsps], 200, 'hspnum-bitscore', 'HSPs in OGs', 'bitscore', labels, colors)
hist1(hsps0['bitscore'], 200, 'hspnum-bitscore_all', 'HSPs in OGs', 'bitscore', labels[0], colors[0])
hist1(hsps1['bitscore'], 200, 'hspnum-bitscore_filter1', 'HSPs in OGs', 'bitscore', labels[1], colors[1])
hist1(hsps2['bitscore'], 200, 'hspnum-bitscore_filter2', 'HSPs in OGs', 'bitscore', labels[2], colors[2])

# 3.2 OG histograms
hist3([OG['bitscore'].mean() for OG in OGs], 200, 'OGnum-bitscoremean', 'OGs', 'mean bitscore of HSPs in OG',
      labels, colors, wrap=True)
hist3([OG['bitscore'].var() for OG in OGs], 200, 'OGnum-bitscorevar', 'OGs', 'variance of bitscore of HSPs in OG',
      labels, colors, wrap=True)

# 3.3 OG scatters
scatter1(OGs[0]['bitscore'].mean(), OGs[0]['bitscore'].var(),
         'bitscorevar-bitscoremean_all.png', 'bitscore', labels[0], colors[0])
scatter1(OGs[1]['bitscore'].mean(), OGs[1]['bitscore'].var(),
         'bitscorevar-bitscoremean_filter1.png', 'bitscore', labels[1], colors[1])
scatter1(OGs[2]['bitscore'].mean(), OGs[2]['bitscore'].var(),
         'bitscorevar-bitscoremean_filter2.png', 'bitscore', labels[2], colors[2])

# 4 PIDENT PLOTS
# 4.1 Hit histograms
hist3([hsp['pident'] for hsp in hsps], 50, 'hspnum-pident', 'HSPs in OGs', 'percent identity', labels, colors)
hist1(hsps0['pident'], 50, 'hspnum-pident_all', 'HSPs in OGs', 'percent identity', labels[0], colors[0])
hist1(hsps1['pident'], 50, 'hspnum-pident_filter1', 'HSPs in OGs', 'percent identity', labels[1], colors[1])
hist1(hsps2['pident'], 50, 'hspnum-pident_filter2', 'HSPs in OGs', 'percent identity', labels[2], colors[2])

# 4.2 OG histograms
hist3([OG['pident'].mean() for OG in OGs], 50, 'OGnum-pidentmean', 'OGs', 'mean percent identity of HSPs in OG',
      labels, colors, wrap=True)
hist3([OG['pident'].var() for OG in OGs], 50, 'OGnum-pidentvar', 'OGs', 'variance of percent identity of HSPs in OG',
      labels, colors, wrap=True)

# 4.3 OG scatters
scatter1(OGs[0]['pident'].mean(), OGs[0]['pident'].var(),
         'pidentvar-pidentmean_all.png', 'percent identity', labels[0], colors[0])
scatter1(OGs[1]['pident'].mean(), OGs[1]['pident'].var(),
         'pidentvar-pidentmean_filter1.png', 'percent identity', labels[1], colors[1])
scatter1(OGs[2]['pident'].mean(), OGs[2]['pident'].var(),
         'pidentvar-pidentmean_filter2.png', 'percent identity', labels[2], colors[2])

# 5 NQA HISTOGRAMS
# 5.1 Hit histograms
hist3([hsp['nqa'] for hsp in hsps], 50, 'hspnum-nqa', 'HSPs in OGs', 'number of query aligned', labels, colors, wrap=True)
hist1(hsps0['nqa'], 50, 'hspnum-nqa_all', 'HSPs in OGs', 'number of query aligned', labels[0], colors[0], wrap=True)
hist1(hsps1['nqa'], 50, 'hspnum-nqa_filter1', 'HSPs in OGs', 'number of query aligned', labels[1], colors[1], wrap=True)
hist1(hsps2['nqa'], 50, 'hspnum-nqa_filter2', 'HSPs in OGs', 'number of query aligned', labels[2], colors[2], wrap=True)

# 5.2 OG histograms
hist3([OG['nqa'].mean() for OG in OGs], 50, 'OGnum-nqamean', 'OGs', 'mean number of query aligned of HSPs in OG',
      labels, colors, wrap=True)
hist3([OG['nqa'].var() for OG in OGs], 50, 'OGnum-nqavar', 'OGs', 'variance of number of query aligned of HSPs in OG',
      labels, colors, wrap=True)

# 5.3 OG scatters
scatter1(OGs[0]['nqa'].mean(), OGs[0]['nqa'].var(),
         'nqavar-nqamean_all.png', 'number of query aligned', labels[0], colors[0])
scatter1(OGs[1]['nqa'].mean(), OGs[1]['nqa'].var(),
         'nqavar-nqamean_filter1.png', 'number of query aligned', labels[1], colors[1])
scatter1(OGs[2]['nqa'].mean(), OGs[2]['nqa'].var(),
         'nqavar-nqamean_filter2.png', 'number of query aligned', labels[2], colors[2])

# 6 FQA HISTOGRAMS
# 6.1 Hit histograms
hist3([hsp['fqa'] for hsp in hsps], 50, 'hspnum-fqa', 'HSPs in OGs', 'fraction of query aligned', labels, colors, wrap=True)
hist1(hsps0['fqa'], 50, 'hspnum-fqa_all', 'HSPs in OGs', 'fraction of query aligned', labels[0], colors[0], wrap=True)
hist1(hsps1['fqa'], 50, 'hspnum-fqa_filter1', 'HSPs in OGs', 'fraction of query aligned', labels[1], colors[1], wrap=True)
hist1(hsps2['fqa'], 50, 'hspnum-fqa_filter2', 'HSPs in OGs', 'fraction of query aligned', labels[2], colors[2], wrap=True)

# 6.2 OG histograms
hist3([OG['fqa'].mean() for OG in OGs], 50, 'OGnum-fqamean', 'OGs', 'mean fraction of query aligned of HSPs in OG',
      labels, colors, wrap=True)
hist3([OG['fqa'].var() for OG in OGs], 50, 'OGnum-fqavar', 'OGs', 'variance of fraction of query aligned of HSPs in OG',
      labels, colors, wrap=True)

# 6.3 OG scatters
scatter1(OGs[0]['fqa'].mean(), OGs[0]['fqa'].var(),
         'fqavar-fqamean_all.png', 'fraction of query aligned', labels[0], colors[0])
scatter1(OGs[1]['fqa'].mean(), OGs[1]['fqa'].var(),
         'fqavar-fqamean_filter1.png', 'fraction of query aligned', labels[1], colors[1])
scatter1(OGs[2]['fqa'].mean(), OGs[2]['fqa'].var(),
         'fqavar-fqamean_filter2.png', 'fraction of query aligned', labels[2], colors[2])

# 7 EDGES
edgenums = [hsp[['qgnid', 'sgnid', 'OGid']].drop_duplicates().groupby('OGid').size() / 2 for hsp in hsps]
gnidnums = [OG['qgnid'].nunique() for OG in OGs]
edgefracs = [2*edgenum / (gnidnum*(gnidnum-1)) for edgenum, gnidnum in zip(edgenums, gnidnums)]

# 7.1 Edge number histograms
hist3(edgenums, 100, 'OGnum-edgenum', 'OGs', 'number of edges', labels, colors)
hist1(edgenums[0], 100, 'OGnum-edgenum_all', 'OGs', 'number of edges', labels[0], colors[0])
hist1(edgenums[1], 50, 'OGnum-edgenum_filter1', 'OGs', 'number of edges', labels[1], colors[1])
hist1(edgenums[2], 50, 'OGnum-edgenum_filter2', 'OGs', 'number of edges', labels[2], colors[2])

# 7.2 Edge fraction histograms
hist3(edgefracs, 50, 'OGnum-edgefrac', 'OGs', 'fraction of possible edges', labels, colors)
hist1(edgefracs[0], 50, 'OGnum-edgefrac_all', 'OGs', 'fraction of possible edges', labels[0], colors[0])
hist1(edgefracs[1], 50, 'OGnum-edgefrac_filter1', 'OGs', 'fraction of possible edges', labels[1], colors[1])
hist1(edgefracs[2], 50, 'OGnum-edgefrac_filter2', 'OGs', 'fraction of possible edges', labels[2], colors[2])

# 8 CORRELATIONS
# 8.1 Including His OGs
scatter2(gnidnums[0], OGs[0]['bitscore'].mean(), 'bitscore-OGgnnum_all', 'Mean bitscore of HSPs in OG')
scatter2(gnidnums[0], OGs[0]['pident'].mean(), 'pident-OGgnnum_all', 'Mean percent identity of HSPs in OG')
scatter2(gnidnums[0], OGs[0]['fqa'].mean(), 'fqa-OGgnnum_all', 'Mean fraction of query aligned of HSPs in OG')
scatter2(gnidnums[0], edgenums[0], 'edgenum-OGgnnum_all', 'Number of edges in OG')
scatter2(gnidnums[0], edgefracs[0], 'edgefrac-OGgnnum_all', 'Fraction of possible edges in OG')

# 11.2 Excluding His OGs
hsp3 = hsps0[~hsps0['OGid'].isin(['0122', '0044', '0059', '004b', '011e'])]
OG3 = hsp3.groupby('OGid')
edgenum = hsp3[['qgnid', 'sgnid', 'OGid']].drop_duplicates().groupby('OGid').size() / 2
gnidnum = OG3['qgnid'].nunique()
edgefrac = 2*edgenum / (gnidnum*(gnidnum-1))

scatter2(gnidnum, OG3['bitscore'].mean(), 'bitscore-OGgnnum_his', 'Mean bitscore of HSPs in OG')
scatter2(gnidnum, OG3['pident'].mean(), 'pident-OGgnnum_his', 'Mean percent identity of HSPs in OG')
scatter2(gnidnum, OG3['fqa'].mean(), 'fqa-OGgnnum_his', 'Mean fraction of query aligned of HSPs in OG')
scatter2(gnidnum, edgenum, 'edgenum-OGgnnum_his', 'Number of edges in OG')
scatter2(gnidnum, edgefrac, 'edgefrac-OGgnnum_his', 'Fraction of possible edges in OG')

"""
DEPENDENCIES
../blast2hsps/blast2hsps.py
    ../blast2hsps/out/hsps.tsv
../hits2reciprocal/hits2reciprocal.py
    ../hits2reciprocal/out/hsps.tsv
../OGs2hsps/OGs2hsps.py
    ../OGs2hsps/out/hsps.tsv
"""