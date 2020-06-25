"""Plot various statistics of the BLAST results."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from math import log10


def hist1(df, bins, x_label, df_label, color, file_label, capital=True):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel((x_label[0].upper() + x_label[1:]) if capital else x_label)
    plt.ylabel('Number of top hits')
    plt.title('Distribution of top hits across ' + x_label)
    plt.legend()
    plt.savefig(f'out/plots/{file_label}.png')
    plt.close()


def hist2(dfs, bins, x_label, df_labels, colors, file_label):
    fig, axs = plt.subplots(2, 1, sharex=True)
    for ax, df, data_label, color in zip(axs, dfs, df_labels, colors):
        ax.hist(df, bins=bins, label=data_label, color=color)
        ax.set_ylabel('Number of top hits')
        ax.legend()
    axs[1].set_xlabel(x_label[0].upper() + x_label[1:])
    fig.suptitle('Distribution of top hits across ' + x_label)

    plt.savefig(f'out/plots/{file_label}.png')
    plt.close()


def bar(counts1, counts2, file_label):
    plt.bar(counts1.keys(), counts1.values(), label='NCBI + FlyBase')
    plt.bar(counts2.keys(), counts2.values(), bottom=[counts1.get(key, 0) for key in counts2.keys()], label=r'Yang $et\ al.$')
    plt.title('Distribution of genes across\nnumber of reciprocal top hits')
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
if not os.path.exists('out/plots'):
    os.makedirs('out/plots')  # Recursive folder creation

# EVALUE PLOTS
eval0 = df0.loc[df0['evalue'] != 0, 'evalue'].apply(log10)
eval1 = df1.loc[df1['evalue'] != 0, 'evalue'].apply(log10)

# Stacked bar
xs = list(range(2))
ys0 = [len(eval0) / len(df0), len(eval1) / len(df1)]  # Fraction of hits with evalue > 0
ys1 = [1 - ys0[0], 1 - ys0[1]]  # Fraction of hits with evalue == 0
plt.bar(xs, ys0, label='non-zero', width=0.25)
plt.bar(xs, ys1, label='zero', width=0.25, bottom=ys0)
plt.xticks(xs, ['all', 'reciprocal'])
plt.xlim((-0.75, 1.75))
plt.ylabel('Fraction of total top hits')
plt.title('Fraction of top hits with zero and non-zero E-values')
plt.legend(bbox_to_anchor=(0.5, -0.1875), loc='lower center', ncol=2)
plt.subplots_adjust(bottom=0.15)
plt.savefig('out/plots/fraction_zero.png')
plt.close()

# Histograms
plt.hist(eval0, bins=200, label='all', color='C0')
plt.hist(eval1, bins=200, label='reciprocal', color='C1')
plt.xlabel('log10(E-value)')
plt.ylabel('Number of top hits')
plt.title('Distribution of top hits across log10(E-value)')
plt.legend()
plt.savefig('out/plots/hist_evalue.png')
plt.close()

hist1(eval0, 200, 'log10(E-value)', 'all', 'C0', 'hist_evalue_all', capital=False)
hist1(eval0, 200, 'log10(E-value)', 'reciprocal', 'C1', 'hist_evalue_reciprocal', capital=False)

# BITSCORE HISTOGRAMS
hist2([df0['bitscore'], df1['bitscore']], 200, 'bitscore', ['all', 'reciprocal'], ['C0', 'C1'], 'hist_bitscore')
hist1(df0['bitscore'], 200, 'bitscore', 'all', 'C0', 'hist_bitscore_all')
hist1(df1['bitscore'], 200, 'bitscore', 'reciprocal', 'C1', 'hist_bitscore_reciprocal')

# PIDENT HISTOGRAMS
hist2([df0['pident'], df1['pident']], 50, 'percent identity', ['all', 'reciprocal'], ['C0', 'C1'], 'hist_pident')
hist1(df0['pident'], 50, 'percent identity', 'all', 'C0', 'hist_pident_all')
hist1(df1['pident'], 50, 'percent identity', 'reciprocal', 'C1', 'hist_pident_reciprocal')

# FRACTION ALIGNED HISTOGRAMS
fali0 = (df0['qend'] - df0['qstart'] + 1) / df0['qlen']
fali1 = (df1['qend'] - df1['qstart'] + 1) / df1['qlen']

hist2([fali0, fali1], 50, 'fraction of query aligned', ['all', 'reciprocal'], ['C0', 'C1'], 'hist_fali')
hist1(fali0, 50, 'fraction of query aligned', 'all', 'C0', 'hist_fali_all')
hist1(fali1, 50, 'fraction of query aligned', 'reciprocal', 'C1', 'hist_fali_reciprocal')

# TOP HITS
for data_label, df in [('all', df0), ('reciprocal', df1)]:
    # Make top hits output directory
    if not os.path.exists(f'out/hitnum_{data_label}'):
        os.mkdir(f'out/hitnum_{data_label}')

    # PPID
    sppid_hitnum = df['sppid'].value_counts().rename('sppid_hitnum').to_frame()
    ids = df.loc[:, ['sppid', 'sgnid', 'sspid']].drop_duplicates().set_index('sppid')
    sppid_hitnum = sppid_hitnum.join(ids)

    sppid_hitnum.to_csv(f'out/hitnum_{data_label}/sppids.tsv', sep='\t', index_label='sppid')
    sppid_hitnum_dmel = sppid_hitnum.loc[sppid_hitnum['sspid'] == 'dmel', :]
    sppid_hitnum_dmel.to_csv(f'out/hitnum_{data_label}/sppids_dmel.tsv', sep='\t', index_label='sppid')

    # GNID
    gnid_pairs = df[['qgnid', 'sgnid']].drop_duplicates()
    sgnid_hitnum = gnid_pairs['sgnid'].value_counts().rename('sgnid_hitnum').to_frame()
    ids = df.loc[:, ['sgnid', 'sspid']].drop_duplicates().set_index('sgnid')
    sgnid_hitnum = sgnid_hitnum.join(ids)

    sgnid_hitnum.to_csv(f'out/hitnum_{data_label}/sgnids.tsv', sep='\t', index_label='sgnid')
    sgnid_hitnum_dmel = sgnid_hitnum.loc[sgnid_hitnum['sspid'] == 'dmel', :]
    sgnid_hitnum_dmel.to_csv(f'out/hitnum_{data_label}/sgnids_dmel.tsv', sep='\t', index_label='sgnid')

    # Correlation of gene hits with number of associated polypeptides
    gnid_ppidnum = pd.read_csv('../genome_stats/out/gnid_ppidnum.tsv', sep='\t', index_col='gnid')
    corr = sgnid_hitnum.join(gnid_ppidnum)
    yo_gns = (corr['spid'] == 'dyak') | (corr['spid'] == 'dpse')

    plt.scatter(corr['ppidnum'], corr['sgnid_hitnum'],
                alpha=0.5, s=10, edgecolors='none')
    plt.xlabel('Number of polypeptides associated with gene')
    plt.ylabel('Number of hits to gene')
    plt.title('Correlation of number of hits to gene\nwith number of associated polypeptides')
    plt.savefig(f'out/plots/scatter_hitnum-ppidnum_{data_label}.png')
    plt.close()

    fig, ax = plt.subplots()
    ax.scatter(corr.loc[~yo_gns, 'ppidnum'], corr.loc[~yo_gns, 'sgnid_hitnum'],
                label='NCBI + FlyBase', alpha=0.5, s=10, edgecolors='none')
    ax.scatter(corr.loc[yo_gns, 'ppidnum'], corr.loc[yo_gns, 'sgnid_hitnum'],
                label=r'Yang $et\ al.$', alpha=0.5, s=10, edgecolors='none')
    ax.set_xlabel('Number of polypeptides associated with gene')
    ax.set_ylabel('Number of hits to gene')
    ax.set_title('Correlation of number of hits to gene\nwith number of associated polypeptides')
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
DEPENDENCIES
../genome_stats/genome_stats.py
    ../genome_stats/out/gnid_ppidnum.tsv
./blast_stats_extract.py
    ./out/ggraph.tsv
"""