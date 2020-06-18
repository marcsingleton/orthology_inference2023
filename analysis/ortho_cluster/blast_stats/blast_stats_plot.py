"""Plot various statistics of the BLAST results."""

import matplotlib.pyplot as plt
import os
import pandas as pd
from math import log10


def hist1(df, bins, x_label, df_label, color, file_label):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel(x_label)
    plt.ylabel('Count')
    plt.title(x_label + ' of top BLAST hits')
    plt.legend()
    plt.savefig(f'out/plots/{file_label}.png')
    plt.close()


def hist2(dfs, bins, x_label, df_labels, colors, file_label):
    fig, axs = plt.subplots(2, 1, sharex=True)
    for ax, df, data_label, color in zip(axs, dfs, df_labels, colors):
        ax.hist(df, bins=bins, label=data_label, color=color)
        ax.set_ylabel('Count')
        ax.legend()
    axs[1].set_xlabel(x_label)
    fig.suptitle(x_label + ' of top BLAST hits')

    plt.savefig(f'out/plots/{file_label}.png')
    plt.close()


def bar(counts1, counts2, file_label):
    plt.bar(counts1.keys(), counts1.values(), label='NCBI')
    plt.bar(counts2.keys(), counts2.values(), bottom=[counts1.get(key, 0) for key in counts2.keys()], label=r'Yang $et\ al.$')
    plt.title('Distribution of number of reciprocal hits to genes')
    plt.xlabel('Number of reciprocal hits to gene')
    plt.ylabel('Count')
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
if not os.path.exists(f'out/plots'):
    os.makedirs(f'out/plots')  # Recursive folder creation

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
plt.ylabel('Fraction')
plt.title('Fraction of Hits with Zero and Non-zero E-values')
plt.legend(bbox_to_anchor=(0.5, -0.1875), loc='lower center', ncol=2)
plt.subplots_adjust(bottom=0.15)
plt.savefig('out/plots/fraction_reciprocal.png')
plt.close()

# Histograms
plt.hist(eval0, bins=200, label='all', color='C0')
plt.hist(eval1, bins=200, label='reciprocal', color='C1')
plt.xlabel('log10(E-value)')
plt.ylabel('Count')
plt.title('log10(E-value) of top BLAST hits')
plt.legend()
plt.savefig('out/plots/hist_evalue.png')
plt.close()

hist1(eval0, 200, 'log10(E-value)', 'all', 'C0', 'hist_evalue_all')
hist1(eval0, 200, 'log10(E-value)', 'reciprocal', 'C1', 'hist_evalue_reciprocal')

# BITSCORE HISTOGRAMS
hist2([df0['bitscore'], df1['bitscore']], 200, 'Bitscore', ['all', 'reciprocal'], ['C0', 'C1'], 'hist_bitscore')
hist1(df0['bitscore'], 200, 'Bitscore', 'all', 'C0', 'hist_bitscore_all')
hist1(df1['bitscore'], 200, 'Bitscore', 'reciprocal', 'C1', 'hist_bitscore_reciprocal')

# PIDENT HISTOGRAMS
hist2([df0['pident'], df1['pident']], 50, 'Percent Identity', ['all', 'reciprocal'], ['C0', 'C1'], 'hist_pident')
hist1(df0['pident'], 50, 'Percent Identity', 'all', 'C0', 'hist_pident_all')
hist1(df1['pident'], 50, 'Percent Identity', 'reciprocal', 'C1', 'hist_pident_reciprocal')

# FRACTION ALIGNED HISTOGRAMS
fali0 = (df0['qend'] - df0['qstart'] + 1) / df0['qlen']
fali1 = (df1['qend'] - df1['qstart'] + 1) / df1['qlen']

hist2([fali0, fali1], 50, 'Fraction of Query Aligned', ['all', 'reciprocal'], ['C0', 'C1'], 'hist_fali')
hist1(fali0, 50, 'Fraction of Query Aligned', 'all', 'C0', 'hist_fali_all')
hist1(fali1, 50, 'Fraction of Query Aligned', 'reciprocal', 'C1', 'hist_fali_reciprocal')

# TOP HITS
for data_label, df in [('all', df0), ('reciprocal', df1)]:
    # Make top hits output directory
    if not os.path.exists(f'out/tophits_{data_label}'):
        os.mkdir(f'out/tophits_{data_label}')

    # PPID
    counts = df['sppid'].value_counts().rename('count').to_frame()
    ids = df.loc[:, ['sppid', 'sgnid', 'sspid']].drop_duplicates().set_index('sppid')
    top_sppids = counts.join(ids)
    top_sppids.to_csv(f'out/tophits_{data_label}/sppids.tsv', sep='\t', index_label='sppid')

    counts = df.loc[df['sspid'] == 'dmel', 'sppid'].value_counts().rename('count').to_frame()
    ids = df.loc[:, ['sppid', 'sgnid']].drop_duplicates().set_index('sppid')
    top_sppids = counts.join(ids)
    top_sppids.to_csv(f'out/tophits_{data_label}/sppids_dmel.tsv', sep='\t', index_label='sppid')

    # GNID
    gnid_pairs = df[['qgnid', 'sgnid']].drop_duplicates()
    counts = gnid_pairs['sgnid'].value_counts().rename('count').to_frame()
    ids = df.loc[:, ['sgnid', 'sspid']].drop_duplicates().set_index('sgnid')
    top_sgnids = counts.join(ids)
    top_sgnids.to_csv(f'out/tophits_{data_label}/sgnids.tsv', sep='\t', index_label='sgnid')

    gnid_pairs = df.loc[df['sspid'] == 'dmel', ['qgnid', 'sgnid']].drop_duplicates()
    top_sgnids = gnid_pairs['sgnid'].value_counts().rename('count').to_frame()
    top_sgnids.to_csv(f'out/tophits_{data_label}/sgnids_dmel.tsv', sep='\t', index_label='sgnid')

# Bar plot
top_sgnids = pd.read_csv('out/tophits_reciprocal/sgnids.tsv', sep='\t')
yo_gns = (top_sgnids['sspid'] == 'dyak') | (top_sgnids['sspid'] == 'dpse')
yo_counts = top_sgnids.loc[yo_gns, 'count'].value_counts().to_dict()
ncbi_counts = top_sgnids.loc[~yo_gns, 'count'].value_counts().to_dict()

bar(ncbi_counts, yo_counts, 'all')
bar({key: val for key, val in ncbi_counts.items() if key > 10},
    {key: val for key, val in yo_counts.items() if key > 10}, 'g10')
bar({key: val for key, val in ncbi_counts.items() if key <= 10},
    {key: val for key, val in yo_counts.items() if key <= 10}, 'le10')

"""
DEPENDENCIES
./blast_stats_extract.py
    ./out/ggraph.tsv
"""