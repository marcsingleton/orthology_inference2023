"""Plot various statistics of HSPs."""

import multiprocessing as mp
import os
from itertools import permutations
from math import log10

import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
import pandas as pd
from numpy import linspace


# Load functions
def load_hsp(qspid, sspid):
    df = pd.read_table(f'../../ortho_search/blast2hsps/out/hsps/{qspid}/{sspid}.tsv',
                       usecols=dtypes.keys(), dtype=dtypes, memory_map=True)
    r = pd.read_table(f'../../ortho_search/hsps2reciprocal/out/{qspid}/{sspid}.tsv',
                      usecols=['reciprocal'], memory_map=True)
    df = df.join(r)

    df['pident'] = df['nident'] / df['length']
    df['nqa'] = df['qend'] - df['qstart'] + 1
    df['fqa'] = df['nqa'] / df['qlen']
    df['nsa'] = df['send'] - df['sstart'] + 1
    df['fsa'] = df['nsa'] / df['slen']
    df['logevalue'] = df['evalue'].apply(lambda x: log10(x) if x > 0 else -180)

    df['qspid'] = qspid
    df['qspid'] = df['qspid'].astype('category')
    df['sspid'] = sspid
    df['sspid'] = df['sspid'].astype('category')

    return df


# Plot functions
def hist1(df, bins, data_label, file_label, x_label, df_label, color):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel(x_label)
    plt.ylabel('Number of HSPs')
    plt.title(data_label)
    plt.legend()
    plt.savefig(f'out/blast_{data_label}/hist_{file_label}.png')
    plt.close()


def hist2_1(dfs, bins, data_label, file_label, x_label, df_labels, colors):
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for df, df_label, color in zip(dfs, df_labels, colors):
        plt.hist(df, bins=b, label=df_label, color=color)
    plt.xlabel(x_label)
    plt.ylabel('Number of HSPs')
    plt.title(data_label)
    plt.legend()
    plt.savefig(f'out/blast_{data_label}/hist_{file_label}_2-1.png')
    plt.close()


def hist2_2(dfs, bins, data_label, file_label, x_label, df_labels, colors):
    fig, axs = plt.subplots(2, 1, sharex=True)
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for ax, df, df_label, color in zip(axs, dfs, df_labels, colors):
        ax.hist(df, bins=b, label=df_label, color=color)
        ax.set_ylabel('Number of HSPs')
        ax.legend()
    axs[1].set_xlabel(x_label)
    fig.suptitle(data_label)
    fig.savefig(f'out/blast_{data_label}/hist_{file_label}_2-2.png')
    plt.close()


def bar(counts, data_label, file_label):
    plt.bar(counts.keys(), counts.values(), width=1)
    plt.xlabel('Number of HSPs to protein')
    plt.ylabel('Number of proteins')
    plt.title(f'{data_label} index HSPs')
    plt.savefig(f'out/hsps_{data_label}/hist_ppidnum-hspnum_{file_label}.png')
    plt.close()


def scatter(x, y, data_label, file_label, x_label, y_label, df_label, color):
    fig, ax = plt.subplots()
    ax.scatter(x, y, label=df_label, color=color, alpha=0.2, s=10, edgecolors='none')
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(data_label)
    leg = ax.legend(markerscale=2)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    fig.savefig(f'out/blast_{data_label}/scatter_{file_label}.png')
    plt.close()


dtypes = {'qppid': 'string', 'qgnid': 'string',
          'sppid': 'string', 'sgnid': 'string',
          'length': int, 'nident': int,
          'qlen': int, 'qstart': int, 'qend': int,
          'slen': int, 'sstart': int, 'send': int,
          'evalue': float, 'bitscore': float,
          'index_hsp': bool, 'disjoint': bool, 'compatible': bool}
num_processes = 4

if __name__ == '__main__':
    # Load genomes
    spids = []
    with open('../config/genomes.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            spids.append(line.rstrip('\n').split('\t')[0])

    # Load data
    with mp.Pool(processes=num_processes) as pool:
        hsps0 = pd.concat(pool.starmap(load_hsp, permutations(spids, 2)))
        hsps1 = hsps0[hsps0['compatible']]
        hsps2 = hsps0[hsps0['disjoint']]
        hsps3 = hsps0[hsps0['index_hsp']]

    if not os.path.exists('out/'):
        os.mkdir('out/')

    # 0 UNFILTERED, COMPATIBLE, DISJOINT, AND INDEX SIZES
    plt.bar(['unfiltered', 'compatible', 'disjoint', 'index'], [len(hsps0), len(hsps1), len(hsps2), len(hsps3)], color=('C7', 'C6', 'C3', 'C0'), width=0.5)
    plt.xlabel('HSP subset')
    plt.ylabel('Number of HSPs')
    plt.savefig('out/bar_hspnum-filter.png')
    plt.close()

    # 1 BLAST METRICS
    plots = [('unfiltered', hsps0, ('C7', 'C4')),
             ('compatible', hsps1, ('C6', 'C5')),
             ('disjoint', hsps2, ('C3', 'C2')),
             ('index', hsps3, ('C0', 'C1'))]
    for data_label, hsps, colors in plots:
        if not os.path.exists(f'out/blast_{data_label}/'):
            os.mkdir(f'out/blast_{data_label}/')

        # Subset HSPs into non-reciprocal and reciprocal sets
        df0 = hsps
        df1 = hsps[hsps['reciprocal']]
        labels = ['all', 'reciprocal']

        # 1.1 FILTER PLOT
        plt.bar(labels, [len(df0), len(df1)], color=colors, width=0.25)
        plt.xlim((-0.75, 1.75))
        plt.ylabel('Number of HSPs')
        plt.title(data_label)
        plt.savefig(f'out/blast_{data_label}/bar_reciprocal_filter.png')
        plt.close()
        print(f'Fraction of {data_label} HSPs reciprocal:', round(len(df1) / len(df0), 3))

        # 1.2 QUALITY METRICS
        evalue0 = df0.loc[df0['evalue'] != 0, 'logevalue']
        evalue1 = df1.loc[df1['evalue'] != 0, 'logevalue']

        # 1.2.1 Stacked bars of zero and non-zero E-values in all and reciprocal data sets
        xs = [0, 1]
        ys1 = [len(evalue0) / len(df0), len(evalue1) / len(df1)]  # Fraction of HSPs with evalue > 0
        ys2 = [1 - ys1[0], 1 - ys1[1]]  # Fraction of HSPs with evalue == 0
        plt.bar(xs, ys1, label='non-zero', width=0.25)
        plt.bar(xs, ys2, label='zero', width=0.25, bottom=ys1)
        plt.xticks(xs, labels)
        plt.xlim((-0.75, 1.75))
        plt.ylabel('Fraction of HSPs')
        plt.title(f'Fraction of {data_label} HSPs with zero and non-zero E-values')
        plt.legend(bbox_to_anchor=(0.5, -0.1875), loc='lower center', ncol=2)
        plt.subplots_adjust(bottom=0.15)
        plt.savefig(f'out/blast_{data_label}/bar_evalue_zero.png')
        plt.close()

        # 1.2.2 Histograms of non-zero E-values
        hist2_1([evalue0, evalue1], 200, data_label, 'evalue', 'log10(E-value)', labels, colors)
        hist2_2([evalue0, evalue1], 200, data_label, 'evalue', 'log10(E-value)', labels, colors)
        hist1(evalue0, 200, data_label, 'evalue_all', 'log10(E-value)', labels[0], colors[0])
        hist1(evalue1, 200, data_label, 'evalue_reciprocal', 'log10(E-value)', labels[1], colors[1])

        # 1.2.3 Scatters of E-value with other metrics
        # 1.2.3.1 FQA-E-value (heatmap)
        plt.hist2d(df0['logevalue'], df0['fqa'], bins=50, norm=mpl_colors.LogNorm())
        plt.xlabel('log10(E-value)')
        plt.ylabel('Fraction of query aligned')
        plt.colorbar()
        plt.savefig(f'out/blast_{data_label}/hist2d_fqa-evalue_all.png')
        plt.close()

        plt.hist2d(df1['logevalue'], df1['fqa'], bins=50, norm=mpl_colors.LogNorm())
        plt.xlabel('log10(E-value)')
        plt.ylabel('Fraction of query aligned')
        plt.colorbar()
        plt.savefig(f'out/blast_{data_label}/hist2d_fqa-evalue_reciprocal.png')
        plt.close()

        # 1.2.3.2 NQAmin-E-value
        groups0 = df0.groupby('logevalue')
        x = groups0['nqa'].min()
        scatter(x.index, x.values, data_label, 'nqamin-evalue_all', 'log10(E-value)', 'Minimum number of residues aligned', labels[0], colors[0])

        groups1 = df1.groupby('logevalue')
        x = groups1['nqa'].min()
        fig, ax = plt.subplots()
        scatter(x.index, x.values, data_label, 'nqamin-evalue_reciprocal', 'log10(E-value)', 'Minimum number of residues aligned', labels[1], colors[1])

        # 1.2.3.3 NQAmax-E-value
        groups0 = df0.groupby('logevalue')
        x = groups0['nqa'].max()
        scatter(x.index, x.values, data_label, 'nqamax-evalue_all', 'log10(E-value)', 'Maximum number of residues aligned', labels[0], colors[0])

        groups1 = df1.groupby('logevalue')
        x = groups1['nqa'].max()
        scatter(x.index, x.values, data_label, 'nqamax-evalue_reciprocal', 'log10(E-value)', 'Maximum number of residues aligned', labels[1], colors[1])

        # 1.2.3.4 Bitscoremin-E-value
        x = groups0['bitscore'].min()
        scatter(x.index, x.values, data_label, 'bitscoremin-evalue_all', 'log10(E-value)', 'Minimum bitscore', labels[0], colors[0])

        x = groups1['bitscore'].min()
        scatter(x.index, x.values, data_label, 'bitscoremin-evalue_reciprocal', 'log10(E-value)', 'Minimum bitscore', labels[1], colors[1])

        # 1.2.4 Bitscore histograms
        hist2_1([df0['bitscore'], df1['bitscore']], 200, data_label, 'bitscore', 'Bitscore', labels, colors)
        hist2_2([df0['bitscore'], df1['bitscore']], 200, data_label, 'bitscore', 'Bitscore', labels, colors)
        hist1(df0['bitscore'], 200, data_label, 'bitscore_all', 'Bitscore', labels[0], colors[0])
        hist1(df1['bitscore'], 200, data_label, 'bitscore_reciprocal', 'Bitscore', labels[1], colors[1])

        # 1.2.5 Pident histograms
        hist2_1([df0['pident'], df1['pident']], 50, data_label, 'pident', 'Percent identity', labels, colors)
        hist2_2([df0['pident'], df1['pident']], 50, data_label, 'pident', 'Percent identity', labels, colors)
        hist1(df0['pident'], 50, data_label, 'pident_all', 'Percent identity', labels[0], colors[0])
        hist1(df1['pident'], 50, data_label, 'pident_reciprocal', 'Percent identity', labels[1], colors[1])

        # 1.2.6 NQA histograms
        hist2_1([df0['nqa'], df1['nqa']], 200, data_label, 'nqa', 'Number of query residues aligned', labels, colors)
        hist2_2([df0['nqa'], df1['nqa']], 200, data_label, 'nqa', 'Number of query residues aligned', labels, colors)
        hist1(df0['nqa'], 200, data_label, 'nqa_all', 'Number of query residues aligned', labels[0], colors[0])
        hist1(df1['nqa'], 200, data_label, 'nqa_reciprocal', 'Number of query residues aligned', labels[1], colors[1])

        # 1.2.5 FQA histograms
        hist2_1([df0['fqa'], df1['fqa']], 50, data_label, 'fqa', 'Fraction of query aligned', labels, colors)
        hist2_2([df0['fqa'], df1['fqa']], 50, data_label, 'fqa', 'Fraction of query aligned', labels, colors)
        hist1(df0['fqa'], 50, data_label, 'fqa_all', 'Fraction of query aligned', labels[0], colors[0])
        hist1(df1['fqa'], 50, data_label, 'fqa_reciprocal', 'Fraction of query aligned', labels[1], colors[1])

        # 1.2.7 FQA-FSA scatters
        plt.hist2d(df0['fqa'], df0['fsa'], bins=50, norm=mpl_colors.LogNorm())
        plt.xlabel('Fraction of query aligned')
        plt.ylabel('Fraction of subject aligned')
        plt.title(data_label)
        plt.colorbar()
        plt.savefig(f'out/blast_{data_label}/hist2d_fsa-fqa_all.png')
        plt.close()

        plt.hist2d(df1['fqa'], df1['fsa'], bins=50, norm=mpl_colors.LogNorm())
        plt.xlabel('Fraction of query aligned')
        plt.ylabel('Fraction of subject aligned')
        plt.title(data_label)
        plt.colorbar()
        plt.savefig(f'out/blast_{data_label}/hist2d_fsa-fqa_reciprocal.png')
        plt.close()

    # 2 HSP METRICS
    for data_label, hsps in [('all', hsps3), ('reciprocal', hsps3[hsps3['reciprocal']])]:
        if not os.path.exists(f'out/hsps_{data_label}/'):
            os.mkdir(f'out/hsps_{data_label}/')

        # 2.1 PPID RANKS
        ids = hsps.loc[:, ['sppid', 'sgnid', 'sspid']].drop_duplicates().set_index('sppid')

        sppid_hspnum = hsps['sppid'].value_counts().rename('sppid_hspnum').sort_values(ascending=False).to_frame()
        sppid_hspnum = sppid_hspnum.join(ids)
        sppid_hspnum.to_csv(f'out/hsps_{data_label}/sppids.tsv', sep='\t', index_label='sppid')

        sppid_hspnum_dmel = sppid_hspnum.loc[sppid_hspnum['sspid'] == 'dmel', :]
        sppid_hspnum_dmel.to_csv(f'out/hsps_{data_label}/sppids_dmel.tsv', sep='\t', index_label='sppid')

        # 2.2 GNID RANKS
        ids = hsps.loc[:, ['sgnid', 'sspid']].drop_duplicates().set_index('sgnid')

        sgnid_hspnum = hsps.groupby('sgnid')['qgnid'].nunique().rename('sgnid_hspnum').sort_values(ascending=False).to_frame()
        sgnid_hspnum = sgnid_hspnum.join(ids)
        sgnid_hspnum.to_csv(f'out/hsps_{data_label}/sgnids.tsv', sep='\t', index_label='sgnid')

        sgnid_hspnum_dmel = sgnid_hspnum.loc[sgnid_hspnum['sspid'] == 'dmel', :]
        sgnid_hspnum_dmel.to_csv(f'out/hsps_{data_label}/sgnids_dmel.tsv', sep='\t', index_label='sgnid')

        # 2.3 PLOTS
        # 2.3.1 Correlation of gene HSPs with number of associated proteins
        gnid_nums = pd.read_table('../genome_stats/out/gnidnums.tsv', index_col='gnid', dtype={'gnid': 'string'})
        corr = sgnid_hspnum.join(gnid_nums)

        plt.scatter(corr['ppidnum'], corr['sgnid_hspnum'], alpha=0.5, s=10, edgecolors='none')
        plt.xlabel('Number of proteins associated with gene')
        plt.ylabel('Number of HSPs to gene')
        plt.title(f'{data_label} index HSPs')
        plt.savefig(f'out/hsps_{data_label}/scatter_hspnum-ppidnum.png')
        plt.close()

        # 2.3.2 Histograms of genes by number of HSPs
        counts = sppid_hspnum['sppid_hspnum'].value_counts().to_dict()

        bar(counts, data_label, 'all')
        bar({key: value for key, value in counts.items() if key > len(spids)}, data_label, f'{len(spids)}+')
        bar({key: value for key, value in counts.items() if key <= len(spids)}, data_label, f'{len(spids)}-')

"""
OUTPUT

DEPENDENCIES
../../ortho_search/blast2hsps/blast2hsps.py
    ../../ortho_search/blast2hsps/out/hsps/*/*.tsv
../../ortho_search/hsps2reciprocal/hsps2reciprocal.py
    ../../ortho_search/hsps2reciprocal/out/*/*.tsv
../config/genomes.tsv
../genome_stats/genome_stats.py
    ../genome_stats/out/gnidnums.tsv
"""