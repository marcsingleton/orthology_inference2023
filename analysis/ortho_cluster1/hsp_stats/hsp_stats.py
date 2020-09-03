"""Plot various statistics of HSPs."""

import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
import os
import pandas as pd
from math import log10
from numpy import linspace


def hist1(df, bins, file_label, x_label, df_label, color, capital=True, wrap=False):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel(x_label[0].upper() + x_label[1:] if capital else x_label)
    plt.ylabel('Number of best HSPs')
    plt.title('Distribution of best HSPs across' + ('\n' if wrap else ' ') + x_label)
    plt.legend()
    plt.savefig(f'out/blast/hist_{file_label}.png')
    plt.close()


def hist2_1(dfs, bins, file_label, x_label, df_labels, colors, capital=True, wrap=False):
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for df, df_label, color in zip(dfs, df_labels, colors):
        plt.hist(df, bins=b, label=df_label, color=color)
    plt.xlabel(x_label[0].upper() + x_label[1:] if capital else x_label)
    plt.ylabel('Number of best HSPs')
    plt.title('Distribution of best HSPs across' + ('\n' if wrap else ' ') + x_label)
    plt.legend()
    plt.savefig(f'out/blast/hist2_1_{file_label}.png')
    plt.close()


def hist2_2(dfs, bins, file_label, x_label, df_labels, colors, capital=True, wrap=False):
    fig, axs = plt.subplots(2, 1, sharex=True)
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for ax, df, df_label, color in zip(axs, dfs, df_labels, colors):
        ax.hist(df, bins=b, label=df_label, color=color)
        ax.set_ylabel('Number of best HSPs')
        ax.legend()
    axs[1].set_xlabel(x_label[0].upper() + x_label[1:] if capital else x_label)
    fig.suptitle('Distribution of best HSPs across' + ('\n' if wrap else ' ') + x_label)
    fig.savefig(f'out/blast/hist2_2_{file_label}.png')
    plt.close()


def bar_hits(counts1, counts2, data_label, file_label, ax_label):
    plt.bar(counts1.keys(), counts1.values(),
            width=1, label='NCBI + FlyBase')
    plt.bar(counts2.keys(), counts2.values(), bottom=[counts1.get(key, 0) for key in counts2.keys()],
            width=1, label=r'Yang $et\ al.$')
    plt.title(f'Distribution of genes across number of {ax_label}hits')
    plt.xlabel(f'Number of {ax_label}hits to gene')
    plt.ylabel('Number of genes')
    plt.legend()
    plt.savefig(f'out/hits_{data_label}/hist_gnidnum-hitnum_{file_label}.png')
    plt.close()


dtypes = {'qppid': 'string', 'qgnid': 'string', 'qspid': 'string',
          'sppid': 'string', 'sgnid': 'string', 'sspid': 'string',
          'length': int, 'nident': int,
          'qlen': int, 'qstart': int, 'qend': int,
          'evalue': float, 'bitscore': float}
df = pd.read_csv('../blast2hsps/out/hsps.tsv', sep='\t', dtype=dtypes)
r = pd.read_csv('../hsps2reciprocal/out/hsps.tsv', sep='\t', usecols=['reciprocal'])

df['pident'] = df['nident'] / df['length']
df['nqa'] = df['qend'] - df['qstart'] + 1
df['fqa'] = (df['qend'] - df['qstart'] + 1) / df['qlen']
df['logevalue'] = df['evalue'].apply(lambda x: log10(x) if x > 0 else -180)
hsps = df.join(r)

# 1 BLAST METRICS
# Subset ggraph into various categories
df0 = hsps
df1 = hsps.query('reciprocal == True')
labels = ['all', 'reciprocal']
colors = ['C0', 'C1']

# Make plots output directory
if not os.path.exists('out/blast/'):
    os.makedirs('out/blast/')  # Recursive folder creation

# 1.1 FILTER PLOT
plt.bar(labels, [len(df0), len(df1)], color=colors, width=0.25)
plt.xlim((-0.75, 1.75))
plt.ylabel('Number of best HSPs')
plt.savefig('out/blast/bar_reciprocal_filter.png')
plt.close()
print('Fraction reciprocal:', len(df1) / len(df0))

# 1.2 QUALITY METRICS
# 1.2.1 E-value plots
evalue0 = df0.loc[df0['evalue'] != 0, 'logevalue']
evalue1 = df1.loc[df1['evalue'] != 0, 'logevalue']

# 1.2.1.1 Stacked bars of zero and non-zero E-values in all and reciprocal data sets
xs = list(range(2))
ys_g0 = [len(evalue0) / len(df0), len(evalue1) / len(df1)]  # Fraction of hits with evalue > 0
ys_e0 = [1 - ys_g0[0], 1 - ys_g0[1]]  # Fraction of hits with evalue == 0
plt.bar(xs, ys_g0, label='non-zero', width=0.25)
plt.bar(xs, ys_e0, label='zero', width=0.25, bottom=ys_g0)
plt.xticks(xs, labels)
plt.xlim((-0.75, 1.75))
plt.ylabel('Fraction of total best HSPs')
plt.title('Fraction of best HSPs with zero and non-zero E-values')
plt.legend(bbox_to_anchor=(0.5, -0.1875), loc='lower center', ncol=2)
plt.subplots_adjust(bottom=0.15)
plt.savefig('out/blast/bar_evalue_zero.png')
plt.close()

# 1.2.1.2 Histograms of non-zero E-values
hist2_1([evalue0, evalue1], 200, 'evalue', 'log10(E-value)', labels, colors, capital=False)
hist2_2([evalue0, evalue1], 200, 'evalue', 'log10(E-value)', labels, colors, capital=False)
hist1(evalue0, 200, 'evalue_all', 'log10(E-value)', labels[0], colors[0], capital=False)
hist1(evalue1, 200, 'evalue_reciprocal', 'log10(E-value)', labels[1], colors[1], capital=False)

# 1.2.1.3 Scatters of E-value with other metrics
plt.hist2d(df1['logevalue'], df1['fqa'], bins=50, norm=mpl_colors.PowerNorm(0.3))
plt.xlabel('log10(E-value)')
plt.ylabel('Fraction of query aligned')
plt.colorbar()
plt.savefig('out/blast/hist2d_fqa-evalue.png')
plt.close()

g0 = df0.groupby('logevalue')
x = g0['nqa'].min()
fig, ax = plt.subplots()
ax.scatter(x.index, x.values, label=labels[0], color=colors[0], alpha=0.5, s=10, edgecolors='none')
ax.set_xlabel('log10(E-value)')
ax.set_ylabel('Minimum number of residues aligned')
leg = ax.legend(markerscale=2)
for lh in leg.legendHandles:
    lh.set_alpha(1)
fig.savefig('out/blast/scatter_nqamin-evalue_all.png')
plt.close()

g1 = df1.groupby('logevalue')
x = g1['nqa'].min()
fig, ax = plt.subplots()
ax.scatter(x.index, x.values, label=labels[1], color=colors[1], alpha=0.5, s=10, edgecolors='none')
ax.set_xlabel('log10(E-value)')
ax.set_ylabel('Minimum number of residues aligned')
leg = ax.legend(markerscale=2)
for lh in leg.legendHandles:
    lh.set_alpha(1)
fig.savefig('out/blast/scatter_nqamin-evalue_reciprocal.png')
plt.close()

x = g0['bitscore'].min()
fig, ax = plt.subplots()
ax.scatter(x.index, x.values, label=labels[0], color=colors[0], alpha=0.5, s=10, edgecolors='none')
ax.set_xlabel('log10(E-value)')
ax.set_ylabel('Minimum bitscore')
leg = ax.legend(markerscale=2)
for lh in leg.legendHandles:
    lh.set_alpha(1)
fig.savefig('out/blast/scatter_bitscoremin-evalue_all.png')
plt.close()

x = g1['bitscore'].min()
fig, ax = plt.subplots()
ax.scatter(x.index, x.values, label=labels[1], color=colors[1], alpha=0.5, s=10, edgecolors='none')
ax.set_xlabel('log10(E-value)')
ax.set_ylabel('Minimum bitscore')
leg = ax.legend(markerscale=2)
for lh in leg.legendHandles:
    lh.set_alpha(1)
fig.savefig('out/blast/scatter_bitscoremin-evalue_reciprocal.png')
plt.close()

# 1.2.2 Bitscore histograms
hist2_1([df0['bitscore'], df1['bitscore']], 200, 'bitscore', 'bitscore', labels, colors)
hist2_2([df0['bitscore'], df1['bitscore']], 200, 'bitscore', 'bitscore', labels, colors)
hist1(df0['bitscore'], 200, 'bitscore_all', 'bitscore', labels[0], colors[0])
hist1(df1['bitscore'], 200, 'bitscore_reciprocal', 'bitscore', labels[1], colors[1])

# 1.2.3 Pident histograms
hist2_1([df0['pident'], df1['pident']], 50, 'pident', 'percent identity', labels, colors)
hist2_2([df0['pident'], df1['pident']], 50, 'pident', 'percent identity', labels, colors)
hist1(df0['pident'], 50, 'pident_all', 'percent identity', labels[0], colors[0])
hist1(df1['pident'], 50, 'pident_reciprocal', 'percent identity', labels[1], colors[1])

# 1.2.4 NQA histograms
hist2_1([df0['nqa'], df1['nqa']], 200, 'nqa', 'number of query residues aligned', labels, colors, wrap=True)
hist2_2([df0['nqa'], df1['nqa']], 200, 'nqa', 'number of query residues aligned', labels, colors, wrap=True)
hist1(df0['nqa'], 50, 'nqa_all', 'number of query residues aligned', labels[0], colors[0], wrap=True)
hist1(df1['nqa'], 50, 'nqa_reciprocal', 'number of query residues aligned', labels[1], colors[1], wrap=True)

# 1.2.5 FQA histograms
hist2_1([df0['fqa'], df1['fqa']], 50, 'fqa', 'fraction of query aligned', labels, colors, wrap=True)
hist2_2([df0['fqa'], df1['fqa']], 50, 'fqa', 'fraction of query aligned', labels, colors, wrap=True)
hist1(df0['fqa'], 50, 'fqa_all', 'fraction of query aligned', labels[0], colors[0], wrap=True)
hist1(df1['fqa'], 50, 'fqa_reciprocal', 'fraction of query aligned', labels[1], colors[1], wrap=True)

# 2 HIT METRICS
for data_label, df in [('all', df0), ('reciprocal', df1)]:
    # Make hits output directory
    if not os.path.exists(f'out/hits_{data_label}/'):
        os.mkdir(f'out/hits_{data_label}/')

    # 2.1 PPID RANKS
    ids = df.loc[:, ['sppid', 'sgnid', 'sspid']].drop_duplicates().set_index('sppid')
    YO = df[~df['qspid'].isin(['dpse', 'dyak'])]

    sppid_hitnum = df['sppid'].value_counts().rename('sppid_hitnum').sort_values(ascending=False).to_frame()
    sppid_hitnum = sppid_hitnum.join(ids)
    sppid_hitnum.to_csv(f'out/hits_{data_label}/sppids.tsv', sep='\t', index_label='sppid')

    sppid_hitnum_dmel = sppid_hitnum.loc[sppid_hitnum['sspid'] == 'dmel', :]
    sppid_hitnum_dmel.to_csv(f'out/hits_{data_label}/sppids_dmel.tsv', sep='\t', index_label='sppid')

    sppid_hitnum_YO = YO['sppid'].value_counts().rename('sppid_hitnum').sort_values(ascending=False).to_frame()
    sppid_hitnum_YO = sppid_hitnum_YO.join(ids)
    sppid_hitnum_YO.to_csv(f'out/hits_{data_label}/sppids_YO.tsv', sep='\t', index_label='sppid')

    # 2.2 GNID RANKS
    ids = df.loc[:, ['sgnid', 'sspid']].drop_duplicates().set_index('sgnid')

    sgnid_hitnum = df.groupby('sgnid')['qgnid'].nunique().rename('sgnid_hitnum').sort_values(ascending=False).to_frame()
    sgnid_hitnum = sgnid_hitnum.join(ids)
    sgnid_hitnum.to_csv(f'out/hits_{data_label}/sgnids.tsv', sep='\t', index_label='sgnid')

    sgnid_hitnum_dmel = sgnid_hitnum.loc[sgnid_hitnum['sspid'] == 'dmel', :]
    sgnid_hitnum_dmel.to_csv(f'out/hits_{data_label}/sgnids_dmel.tsv', sep='\t', index_label='sgnid')

    sgnid_hitnum_YO = YO.groupby('sgnid')['qgnid'].nunique().rename('sgnid_hitnum').sort_values(ascending=False).to_frame()
    sgnid_hitnum_YO = sgnid_hitnum_YO.join(ids)
    sgnid_hitnum_YO.to_csv(f'out/hits_{data_label}/sgnids_YO.tsv', sep='\t', index_label='sgnid')

    # 2.3 PLOTS
    ax_label = 'reciprocal ' if data_label == 'reciprocal' else ''

    # Correlation of gene hits with number of associated polypeptides
    gnid_nums = pd.read_csv('../genome_stats/out/gnid_nums.tsv', sep='\t',
                               index_col='gnid', dtype={'gnid': 'string'})
    corr = sgnid_hitnum.join(gnid_nums)
    yo_gns = corr['spid'].isin(['dpse', 'dyak'])

    plt.scatter(corr['ppidnum'], corr['sgnid_hitnum'],
                alpha=0.5, s=10, edgecolors='none')
    plt.xlabel('Number of polypeptides associated with gene')
    plt.ylabel(f'Number of {ax_label}hits to gene')
    plt.title(f'Correlation of number of {ax_label}hits to gene\nwith number of associated polypeptides')
    plt.savefig(f'out/hits_{data_label}/scatter_hitnum-ppidnum_{data_label}.png')
    plt.close()

    fig, ax = plt.subplots()
    ax.scatter(corr.loc[~yo_gns, 'ppidnum'], corr.loc[~yo_gns, 'sgnid_hitnum'],
               label='NCBI + FlyBase', alpha=0.5, s=10, edgecolors='none')
    ax.scatter(corr.loc[yo_gns, 'ppidnum'], corr.loc[yo_gns, 'sgnid_hitnum'],
               label=r'Yang $et\ al.$', alpha=0.5, s=10, edgecolors='none')
    ax.set_xlabel('Number of polypeptides associated with gene')
    ax.set_ylabel(f'Number of {ax_label}hits to gene')
    ax.set_title(f'Correlation of number of {ax_label}hits to gene\nwith number of associated polypeptides')
    leg = ax.legend(markerscale=2)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    fig.savefig(f'out/hits_{data_label}/scatter_hitnum-ppidnum_NCBI-YO_{data_label}.png')
    plt.close()

    # Histograms of genes by number of hits
    yo_gns = sgnid_hitnum['sspid'].isin(['dpse', 'dyak'])
    yo_counts = sgnid_hitnum.loc[yo_gns, 'sgnid_hitnum'].value_counts().to_dict()
    ncbifb_counts = sgnid_hitnum.loc[~yo_gns, 'sgnid_hitnum'].value_counts().to_dict()

    bar_hits(ncbifb_counts, yo_counts, data_label, 'all', ax_label)
    bar_hits({key: val for key, val in ncbifb_counts.items() if key > 10},
             {key: val for key, val in yo_counts.items() if key > 10}, data_label, '10+', ax_label)
    bar_hits({key: val for key, val in ncbifb_counts.items() if key <= 10},
             {key: val for key, val in yo_counts.items() if key <= 10}, data_label, '10-', ax_label)

"""
OUTPUT
Fraction reciprocal: 0.7628523944345292

DEPENDENCIES
../blast2hsps/blast2hsps.py
    ../blast2hsps/hsps.tsv
../genome_stats/genome_stats.py
    ../genome_stats/out/gnid_nums.tsv
../hsps2reciprocal/hsps2reciprocal.py
    ../hsps2reciprocal/out/hsps.tsv
"""