"""Plot various statistics of the BLAST results."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from math import log10
from numpy import linspace


def hist1(df, bins, file_label, x_label, df_label, color, capital=True):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel((x_label[0].upper() + x_label[1:]) if capital else x_label)
    plt.ylabel('Number of hits')
    plt.title('Distribution of hits across ' + x_label)
    plt.legend()
    plt.savefig(f'out/plots/hist_{file_label}.png')
    plt.close()


def hist2_1(dfs, bins, file_label, x_label, df_labels, colors, capital=True):
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for df, data_label, color in zip(dfs, df_labels, colors):
        plt.hist(df, bins=b, label=data_label, color=color)
    plt.xlabel((x_label[0].upper() + x_label[1:]) if capital else x_label)
    plt.ylabel('Number of hits')
    plt.title('Distribution of hits across ' + x_label)
    plt.legend()
    plt.savefig(f'out/plots/hist2_1_{file_label}.png')
    plt.close()


def hist2_2(dfs, bins, file_label, x_label, df_labels, colors, capital=True):
    fig, axs = plt.subplots(2, 1, sharex=True)
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for ax, df, data_label, color in zip(axs, dfs, df_labels, colors):
        ax.hist(df, bins=b, label=data_label, color=color)
        ax.set_ylabel('Number of hits')
        ax.legend()
    axs[1].set_xlabel((x_label[0].upper() + x_label[1:]) if capital else x_label)
    fig.suptitle('Distribution of hits across ' + x_label)
    fig.savefig(f'out/plots/hist2_2_{file_label}.png')
    plt.close()


def bar(counts1, counts2, file_label):
    plt.bar(counts1.keys(), counts1.values(), label='NCBI + FlyBase')
    plt.bar(counts2.keys(), counts2.values(), bottom=[counts1.get(key, 0) for key in counts2.keys()], label=r'Yang $et\ al.$')
    plt.title('Distribution of genes across\nnumber of reciprocal hits')
    plt.xlabel('Number of reciprocal hits to gene')
    plt.ylabel('Number of genes')
    plt.legend()
    plt.savefig(f'out/plots/bar_reciprocal_gnids_{file_label}.png')
    plt.close()


cols = ['qppid', 'qgnid', 'qspid', 'sppid', 'sgnid', 'sspid']
ggraph = pd.read_csv('out/ggraph.tsv', sep='\t', dtype={col: str for col in cols})

# Subset ggraph into various categories
df0 = ggraph
df1 = ggraph.query('reciprocal == True')
print('Fraction reciprocal:', len(df1) / len(df0))

# Make plots output directory
if not os.path.exists('out/plots/'):
    os.makedirs('out/plots/')  # Recursive folder creation

# EVALUE PLOTS
evalue0 = df0.loc[df0['evalue'] != 0, 'evalue'].apply(log10)
evalue1 = df1.loc[df1['evalue'] != 0, 'evalue'].apply(log10)

# Stacked bar
xs = list(range(2))
ys_g0 = [len(evalue0) / len(df0), len(evalue1) / len(df1)]  # Fraction of hits with evalue > 0
ys_e0 = [1 - ys_g0[0], 1 - ys_g0[1]]  # Fraction of hits with evalue == 0
plt.bar(xs, ys_g0, label='non-zero', width=0.25)
plt.bar(xs, ys_e0, label='zero', width=0.25, bottom=ys_g0)
plt.xticks(xs, ['all', 'reciprocal'])
plt.xlim((-0.75, 1.75))
plt.ylabel('Fraction of total hits')
plt.title('Fraction of hits with zero and non-zero E-values')
plt.legend(bbox_to_anchor=(0.5, -0.1875), loc='lower center', ncol=2)
plt.subplots_adjust(bottom=0.15)
plt.savefig('out/plots/fraction_zero.png')
plt.close()

# Histograms
hist2_1([evalue0, evalue1], 200, 'evalue', 'log10(E-value)',
        ['all', 'reciprocal'], ['C0', 'C1'], capital=False)
hist2_2([evalue0, evalue1], 200, 'evalue', 'log10(E-value)',
        ['all', 'reciprocal'], ['C0', 'C1'], capital=False)
hist1(evalue0, 200, 'evalue_all', 'log10(E-value)', 'all', 'C0', capital=False)
hist1(evalue1, 200, 'evalue_reciprocal', 'log10(E-value)', 'reciprocal', 'C1', capital=False)

# BITSCORE HISTOGRAMS
hist2_1([df0['bitscore'], df1['bitscore']], 200, 'bitscore', 'bitscore',
        ['all', 'reciprocal'], ['C0', 'C1'])
hist2_2([df0['bitscore'], df1['bitscore']], 200, 'bitscore', 'bitscore',
        ['all', 'reciprocal'], ['C0', 'C1'])
hist1(df0['bitscore'], 200, 'bitscore_all', 'bitscore', 'all', 'C0')
hist1(df1['bitscore'], 200, 'bitscore_reciprocal', 'bitscore', 'reciprocal', 'C1')

# PIDENT HISTOGRAMS
hist2_1([df0['pident'], df1['pident']], 50, 'pident', 'percent identity',
        ['all', 'reciprocal'], ['C0', 'C1'])
hist2_2([df0['pident'], df1['pident']], 50, 'pident', 'percent identity',
        ['all', 'reciprocal'], ['C0', 'C1'])
hist1(df0['pident'], 50, 'pident_all', 'percent identity', 'all', 'C0')
hist1(df1['pident'], 50, 'pident_reciprocal', 'percent identity', 'reciprocal', 'C1')

# FRACTION ALIGNED HISTOGRAMS
fali0 = (df0['qend'] - df0['qstart'] + 1) / df0['qlen']
fali1 = (df1['qend'] - df1['qstart'] + 1) / df1['qlen']

hist2_1([fali0, fali1], 50, 'fali', 'fraction of query aligned',
        ['all', 'reciprocal'], ['C0', 'C1'])
hist2_2([fali0, fali1], 50, 'fali', 'fraction of query aligned',
        ['all', 'reciprocal'], ['C0', 'C1'])
hist1(fali0, 50, 'fali_all', 'fraction of query aligned', 'all', 'C0')
hist1(fali1, 50, 'fali_reciprocal', 'fraction of query aligned', 'reciprocal', 'C1')

# TOP HITS
for data_label, df in [('all', df0), ('reciprocal', df1)]:
    # Make top hits output directory
    if not os.path.exists(f'out/hitnum_{data_label}/'):
        os.mkdir(f'out/hitnum_{data_label}/')

    # PPID
    ids = df.loc[:, ['sppid', 'sgnid', 'sspid']].drop_duplicates().set_index('sppid')
    YO = df[~df['qspid'].isin(['dpse', 'dyak'])]

    # All hits
    sppid_hitnum = df['sppid'].value_counts().rename('sppid_hitnum').sort_values(ascending=False).to_frame()
    sppid_hitnum = sppid_hitnum.join(ids)
    sppid_hitnum.to_csv(f'out/hitnum_{data_label}/sppids.tsv', sep='\t', index_label='sppid')

    sppid_hitnum_dmel = sppid_hitnum.loc[sppid_hitnum['sspid'] == 'dmel', :]
    sppid_hitnum_dmel.to_csv(f'out/hitnum_{data_label}/sppids_dmel.tsv', sep='\t', index_label='sppid')

    # Hits excluding YO annotations
    sppid_hitnum = YO['sppid'].value_counts().rename('sppid_hitnum').sort_values(ascending=False).to_frame()
    sppid_hitnum = sppid_hitnum.join(ids)
    sppid_hitnum.to_csv(f'out/hitnum_{data_label}/sppids_YO.tsv', sep='\t', index_label='sppid')

    # GNID
    ids = df.loc[:, ['sgnid', 'sspid']].drop_duplicates().set_index('sgnid')

    # All hits
    sgnid_hitnum = df.groupby('sgnid')['qgnid'].nunique().rename('sgnid_hitnum').sort_values(ascending=False).to_frame()
    sgnid_hitnum = sgnid_hitnum.join(ids)
    sgnid_hitnum.to_csv(f'out/hitnum_{data_label}/sgnids.tsv', sep='\t', index_label='sgnid')

    sgnid_hitnum_dmel = sgnid_hitnum.loc[sgnid_hitnum['sspid'] == 'dmel', :]
    sgnid_hitnum_dmel.to_csv(f'out/hitnum_{data_label}/sgnids_dmel.tsv', sep='\t', index_label='sgnid')

    # Hits from Dmel only
    sgnid_hitnum = YO.groupby('sgnid')['qgnid'].nunique().rename('sgnid_hitnum').sort_values(ascending=False).to_frame()
    sgnid_hitnum = sgnid_hitnum.join(ids)
    sgnid_hitnum.to_csv(f'out/hitnum_{data_label}/sgnids_YO.tsv', sep='\t', index_label='sgnid')

    # Correlation of gene hits with number of associated polypeptides
    gnid_ppidnum = pd.read_csv('../genome_stats/out/gnid_ppidnum.tsv', sep='\t', index_col='gnid')
    corr = sgnid_hitnum.join(gnid_ppidnum)
    yo_gns = (corr['spid'] == 'dyak') | (corr['spid'] == 'dpse')
    type = 'reciprocal ' if data_label == 'reciprocal' else ''

    plt.scatter(corr['ppidnum'], corr['sgnid_hitnum'],
                alpha=0.5, s=10, edgecolors='none')
    plt.xlabel('Number of polypeptides associated with gene')
    plt.ylabel(f'Number of {type}hits to gene')
    plt.title(f'Correlation of number of {type}hits to gene\nwith number of associated polypeptides')
    plt.savefig(f'out/plots/scatter_hitnum-ppidnum_{data_label}.png')
    plt.close()

    fig, ax = plt.subplots()
    ax.scatter(corr.loc[~yo_gns, 'ppidnum'], corr.loc[~yo_gns, 'sgnid_hitnum'],
                label='NCBI + FlyBase', alpha=0.5, s=10, edgecolors='none')
    ax.scatter(corr.loc[yo_gns, 'ppidnum'], corr.loc[yo_gns, 'sgnid_hitnum'],
                label=r'Yang $et\ al.$', alpha=0.5, s=10, edgecolors='none')
    ax.set_xlabel('Number of polypeptides associated with gene')
    ax.set_ylabel(f'Number of {type}hits to gene')
    ax.set_title(f'Correlation of number of {type}hits to gene\nwith number of associated polypeptides')
    leg = ax.legend(markerscale=2)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    fig.savefig(f'out/plots/scatter_hitnum-ppidnum_NCBI-YO_{data_label}.png')
    plt.close()

# Bar plot
top_sgnids = pd.read_csv('out/hitnum_reciprocal/sgnids.tsv', sep='\t')
yo_gns = top_sgnids['sspid'].isin(['dpse', 'dyak'])
yo_counts = top_sgnids.loc[yo_gns, 'sgnid_hitnum'].value_counts().to_dict()
ncbifb_counts = top_sgnids.loc[~yo_gns, 'sgnid_hitnum'].value_counts().to_dict()

bar(ncbifb_counts, yo_counts, 'all')
bar({key: val for key, val in ncbifb_counts.items() if key > 10},
    {key: val for key, val in yo_counts.items() if key > 10}, 'g10')
bar({key: val for key, val in ncbifb_counts.items() if key <= 10},
    {key: val for key, val in yo_counts.items() if key <= 10}, 'le10')

"""
OUTPUT
Fraction reciprocal: 0.7628523944345292

DEPENDENCIES
../genome_stats/genome_stats.py
    ../genome_stats/out/gnid_ppidnum.tsv
./blast_stats_extract.py
    ./out/ggraph.tsv
"""