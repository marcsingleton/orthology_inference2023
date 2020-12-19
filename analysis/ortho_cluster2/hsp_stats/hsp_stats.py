"""Plot various statistics of HSPs."""

import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import pandas as pd
from itertools import permutations
from math import log10
from numpy import linspace


# Load functions
def load_hsp(qspid, sspid):
    df = pd.read_csv(f'../../ortho_search/blast2hsps/out/hsps/{qspid}/{sspid}.tsv', sep='\t',
                     usecols=dtypes.keys(), dtype=dtypes, memory_map=True)
    r = pd.read_csv(f'../../ortho_search/hsps2reciprocal/out/{qspid}/{sspid}.tsv', sep='\t',
                    usecols=['reciprocal'], memory_map=True)

    df['pident'] = df['nident'] / df['length']
    df['nqa'] = df['qend'] - df['qstart'] + 1
    df['fqa'] = df['nqa'] / df['qlen']
    df['nsa'] = df['send'] - df['sstart'] + 1
    df['fsa'] = df['nsa'] / df['slen']
    df['logevalue'] = df['evalue'].apply(lambda x: log10(x) if x > 0 else -180)
    return df[df['disjoint']].join(r)


# Plot functions
def hist1(df, bins, data_label, file_label, x_label, df_label, color, capital=True, wrap=False):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel(x_label[0].upper() + x_label[1:] if capital else x_label)
    plt.ylabel(f'Number of {data_label} HSPs')
    plt.title(f'Distribution of {data_label} HSPs across' + ('\n' if wrap else ' ') + x_label)
    plt.legend()
    plt.savefig(f'out/blast_{data_label}/hist_{file_label}.png')
    plt.close()


def hist2_1(dfs, bins, data_label, file_label, x_label, df_labels, colors, capital=True, wrap=False):
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for df, df_label, color in zip(dfs, df_labels, colors):
        plt.hist(df, bins=b, label=df_label, color=color)
    plt.xlabel(x_label[0].upper() + x_label[1:] if capital else x_label)
    plt.ylabel(f'Number of {data_label} HSPs')
    plt.title(f'Distribution of {data_label} HSPs across' + ('\n' if wrap else ' ') + x_label)
    plt.legend()
    plt.savefig(f'out/blast_{data_label}/hist2_1_{file_label}.png')
    plt.close()


def hist2_2(dfs, bins, data_label, file_label, x_label, df_labels, colors, capital=True, wrap=False):
    fig, axs = plt.subplots(2, 1, sharex=True)
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for ax, df, df_label, color in zip(axs, dfs, df_labels, colors):
        ax.hist(df, bins=b, label=df_label, color=color)
        ax.set_ylabel(f'Number of {data_label} HSPs')
        ax.legend()
    axs[1].set_xlabel(x_label[0].upper() + x_label[1:] if capital else x_label)
    fig.suptitle(f'Distribution of {data_label} HSPs across' + ('\n' if wrap else ' ') + x_label)
    fig.savefig(f'out/blast_{data_label}/hist2_2_{file_label}.png')
    plt.close()


def bar_hits(counts, data_label, file_label, ax_label):
    plt.bar(counts.keys(), counts.values(), width=1)
    plt.title(f'Distribution of genes across number of {ax_label}hits')
    plt.xlabel(f'Number of {ax_label}hits to gene')
    plt.ylabel('Number of genes')
    plt.savefig(f'out/hits_{data_label}/hist_gnidnum-hitnum_{file_label}.png')
    plt.close()


dtypes = {'qppid': 'string', 'qgnid': 'string', 'qspid': 'string',
          'sppid': 'string', 'sgnid': 'string', 'sspid': 'string',
          'length': int, 'nident': int,
          'qlen': int, 'qstart': int, 'qend': int,
          'slen': int, 'sstart': int, 'send': int,
          'evalue': float, 'bitscore': float,
          'index_hsp': bool, 'disjoint': bool}
num_processes = 2

if __name__ == '__main__':
    # Parse parameters
    spids = []
    with open('params.tsv') as file:
        fields = file.readline().split()  # Skip header
        for line in file:
            spids.append(line.split()[0])

    # Load data
    with mp.Pool(processes=num_processes) as pool:
        hsps0 = pd.concat(pool.starmap(load_hsp, permutations(spids, 2)))
        hsps1 = hsps0[hsps0['index_hsp']]

    # 0 BEST AND DISJOINT SIZES
    # Make output directory
    if not os.path.exists('out/'):
        os.mkdir('out/')

    plt.bar(['best', 'disjoint'], [len(hsps1), len(hsps0)], color=['C0', 'C3'], width=0.25)
    plt.xlim((-0.75, 1.75))
    plt.ylabel('Number of HSPs')
    plt.savefig('out/bar_best-disjoint.png')
    plt.close()

    # 1 BLAST METRICS
    for data_label, hsps, colors in [('best', hsps1, ['C0', 'C1']), ('disjoint', hsps0, ['C3', 'C6'])]:
        # Subset HSPs into non-reciprocal and reciprocal sets
        df0 = hsps
        df1 = hsps.query('reciprocal == True')
        labels = ['all', 'reciprocal']

        # Make BLAST output directory
        if not os.path.exists(f'out/blast_{data_label}/'):
            os.mkdir(f'out/blast_{data_label}/')

        # 1.1 FILTER PLOT
        plt.bar(labels, [len(df0), len(df1)], color=colors, width=0.25)
        plt.xlim((-0.75, 1.75))
        plt.ylabel(f'Number of {data_label} HSPs')
        plt.savefig(f'out/blast_{data_label}/bar_reciprocal_filter.png')
        plt.close()
        print(f'Fraction of {data_label} HSPs reciprocal:', len(df1) / len(df0))

        # 1.2 QUALITY METRICS
        # 1.2.1 E-value plots
        evalue0 = df0.loc[df0['evalue'] != 0, 'logevalue']
        evalue1 = df1.loc[df1['evalue'] != 0, 'logevalue']

        # 1.2.1.1 Stacked bars of zero and non-zero E-values in all and reciprocal data sets
        xs = list(range(2))
        ys_g0 = [len(evalue0) / len(df0), len(evalue1) / len(df1)]  # Fraction of HSPs with evalue > 0
        ys_e0 = [1 - ys_g0[0], 1 - ys_g0[1]]  # Fraction of HSPs with evalue == 0
        plt.bar(xs, ys_g0, label='non-zero', width=0.25)
        plt.bar(xs, ys_e0, label='zero', width=0.25, bottom=ys_g0)
        plt.xticks(xs, labels)
        plt.xlim((-0.75, 1.75))
        plt.ylabel(f'Fraction of total {data_label} HSPs')
        plt.title(f'Fraction of {data_label} HSPs with zero and non-zero E-values')
        plt.legend(bbox_to_anchor=(0.5, -0.1875), loc='lower center', ncol=2)
        plt.subplots_adjust(bottom=0.15)
        plt.savefig(f'out/blast_{data_label}/bar_evalue_zero.png')
        plt.close()

        # 1.2.1.2 Histograms of non-zero E-values
        hist2_1([evalue0, evalue1], 200, data_label, 'evalue', 'log10(E-value)', labels, colors, capital=False)
        hist2_2([evalue0, evalue1], 200, data_label, 'evalue', 'log10(E-value)', labels, colors, capital=False)
        hist1(evalue0, 200, data_label, 'evalue_all', 'log10(E-value)', labels[0], colors[0], capital=False)
        hist1(evalue1, 200, data_label, 'evalue_reciprocal', 'log10(E-value)', labels[1], colors[1], capital=False)

        # 1.2.1.3 Scatters of E-value with other metrics
        plt.hist2d(df0['logevalue'], df0['fqa'], bins=50, norm=mpl_colors.PowerNorm(0.3))
        plt.xlabel('log10(E-value)')
        plt.ylabel('Fraction of query aligned')
        plt.colorbar()
        plt.savefig(f'out/blast_{data_label}/hist2d_fqa-evalue_all.png')
        plt.close()

        plt.hist2d(df1['logevalue'], df1['fqa'], bins=50, norm=mpl_colors.PowerNorm(0.3))
        plt.xlabel('log10(E-value)')
        plt.ylabel('Fraction of query aligned')
        plt.colorbar()
        plt.savefig(f'out/blast_{data_label}/hist2d_fqa-evalue_reciprocal.png')
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
        fig.savefig(f'out/blast_{data_label}/scatter_nqamin-evalue_all.png')
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
        fig.savefig(f'out/blast_{data_label}/scatter_nqamin-evalue_reciprocal.png')
        plt.close()

        x = g0['bitscore'].min()
        fig, ax = plt.subplots()
        ax.scatter(x.index, x.values, label=labels[0], color=colors[0], alpha=0.5, s=10, edgecolors='none')
        ax.set_xlabel('log10(E-value)')
        ax.set_ylabel('Minimum bitscore')
        leg = ax.legend(markerscale=2)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        fig.savefig(f'out/blast_{data_label}/scatter_bitscoremin-evalue_all.png')
        plt.close()

        x = g1['bitscore'].min()
        fig, ax = plt.subplots()
        ax.scatter(x.index, x.values, label=labels[1], color=colors[1], alpha=0.5, s=10, edgecolors='none')
        ax.set_xlabel('log10(E-value)')
        ax.set_ylabel('Minimum bitscore')
        leg = ax.legend(markerscale=2)
        for lh in leg.legendHandles:
            lh.set_alpha(1)
        fig.savefig(f'out/blast_{data_label}/scatter_bitscoremin-evalue_reciprocal.png')
        plt.close()

        # 1.2.2 Bitscore histograms
        hist2_1([df0['bitscore'], df1['bitscore']], 200, data_label, 'bitscore', 'bitscore', labels, colors)
        hist2_2([df0['bitscore'], df1['bitscore']], 200, data_label, 'bitscore', 'bitscore', labels, colors)
        hist1(df0['bitscore'], 200, data_label, 'bitscore_all', 'bitscore', labels[0], colors[0])
        hist1(df1['bitscore'], 200, data_label, 'bitscore_reciprocal', 'bitscore', labels[1], colors[1])

        # 1.2.3 Pident histograms
        hist2_1([df0['pident'], df1['pident']], 50, data_label, 'pident', 'percent identity', labels, colors)
        hist2_2([df0['pident'], df1['pident']], 50, data_label, 'pident', 'percent identity', labels, colors)
        hist1(df0['pident'], 50, data_label, 'pident_all', 'percent identity', labels[0], colors[0])
        hist1(df1['pident'], 50, data_label, 'pident_reciprocal', 'percent identity', labels[1], colors[1])

        # 1.2.4 NQA histograms
        hist2_1([df0['nqa'], df1['nqa']], 200, data_label, 'nqa', 'number of query residues aligned',
                labels, colors, wrap=True)
        hist2_2([df0['nqa'], df1['nqa']], 200, data_label, 'nqa', 'number of query residues aligned',
                labels, colors, wrap=True)
        hist1(df0['nqa'], 200, data_label, 'nqa_all', 'number of query residues aligned',
              labels[0], colors[0], wrap=True)
        hist1(df1['nqa'], 200, data_label, 'nqa_reciprocal', 'number of query residues aligned',
              labels[1], colors[1], wrap=True)

        # 1.2.5 FQA histograms
        hist2_1([df0['fqa'], df1['fqa']], 50, data_label, 'fqa', 'fraction of query aligned',
                labels, colors, wrap=True)
        hist2_2([df0['fqa'], df1['fqa']], 50, data_label, 'fqa', 'fraction of query aligned',
                labels, colors, wrap=True)
        hist1(df0['fqa'], 50, data_label, 'fqa_all', 'fraction of query aligned',
              labels[0], colors[0], wrap=True)
        hist1(df1['fqa'], 50, data_label, 'fqa_reciprocal', 'fraction of query aligned',
              labels[1], colors[1], wrap=True)

        # 1.2.6 FQA-FSA scatters
        plt.hist2d(df0['fqa'], df0['fsa'], bins=50, norm=mpl_colors.PowerNorm(0.3))
        plt.xlabel('Fraction of query aligned')
        plt.ylabel('Fraction of subject aligned')
        plt.colorbar()
        plt.savefig(f'out/blast_{data_label}/hist2d_fsa-fqa_all.png')
        plt.close()

        plt.hist2d(df1['fqa'], df1['fsa'], bins=50, norm=mpl_colors.PowerNorm(0.3))
        plt.xlabel('Fraction of query aligned')
        plt.ylabel('Fraction of subject aligned')
        plt.colorbar()
        plt.savefig(f'out/blast_{data_label}/hist2d_fsa-fqa_reciprocal.png')
        plt.close()

    # 2 HIT METRICS
    for data_label, hsps in [('all', hsps1), ('reciprocal', hsps1[hsps1['reciprocal']])]:
        # Make hits output directory
        if not os.path.exists(f'out/hits_{data_label}/'):
            os.mkdir(f'out/hits_{data_label}/')

        # 2.1 PPID RANKS
        ids = hsps.loc[:, ['sppid', 'sgnid', 'sspid']].drop_duplicates().set_index('sppid')

        sppid_hitnum = hsps['sppid'].value_counts().rename('sppid_hitnum').sort_values(ascending=False).to_frame()
        sppid_hitnum = sppid_hitnum.join(ids)
        sppid_hitnum.to_csv(f'out/hits_{data_label}/sppids.tsv', sep='\t', index_label='sppid')

        sppid_hitnum_dmel = sppid_hitnum.loc[sppid_hitnum['sspid'] == 'dmel', :]
        sppid_hitnum_dmel.to_csv(f'out/hits_{data_label}/sppids_dmel.tsv', sep='\t', index_label='sppid')

        # 2.2 GNID RANKS
        ids = hsps.loc[:, ['sgnid', 'sspid']].drop_duplicates().set_index('sgnid')

        sgnid_hitnum = hsps.groupby('sgnid')['qgnid'].nunique().rename('sgnid_hitnum').sort_values(ascending=False).to_frame()
        sgnid_hitnum = sgnid_hitnum.join(ids)
        sgnid_hitnum.to_csv(f'out/hits_{data_label}/sgnids.tsv', sep='\t', index_label='sgnid')

        sgnid_hitnum_dmel = sgnid_hitnum.loc[sgnid_hitnum['sspid'] == 'dmel', :]
        sgnid_hitnum_dmel.to_csv(f'out/hits_{data_label}/sgnids_dmel.tsv', sep='\t', index_label='sgnid')

        # 2.3 PLOTS
        ax_label = 'reciprocal ' if data_label == 'reciprocal' else ''

        # 2.3.1 Correlation of gene hits with number of associated polypeptides
        gnid_nums = pd.read_csv('../genome_stats/out/gnid_nums.tsv', sep='\t',
                                index_col='gnid', dtype={'gnid': 'string'})
        corr = sgnid_hitnum.join(gnid_nums)

        plt.scatter(corr['ppidnum'], corr['sgnid_hitnum'],
                    alpha=0.5, s=10, edgecolors='none')
        plt.xlabel('Number of polypeptides associated with gene')
        plt.ylabel(f'Number of {ax_label}hits to gene')
        plt.title(f'Correlation of number of {ax_label}hits to gene\nwith number of associated polypeptides')
        plt.savefig(f'out/hits_{data_label}/scatter_hitnum-ppidnum.png')
        plt.close()

        # 2.3.2 Histograms of genes by number of hits
        counts = sgnid_hitnum['sgnid_hitnum'].value_counts().to_dict()

        bar_hits(counts, data_label, 'all', ax_label)
        bar_hits({key: val for key, val in counts.items() if key > 31}, data_label, '31+', ax_label)
        bar_hits({key: val for key, val in counts.items() if key <= 31}, data_label, '31-', ax_label)

"""
OUTPUT
Fraction of best HSPs reciprocal: 0.9252354199762229
Fraction of disjoint HSPs reciprocal: 0.9273319495519593

DEPENDENCIES
../../ortho_search/blast2hsps/blast2hsps.py
    ../../ortho_search/blast2hsps/out/hsps/*/*.tsv
../../ortho_search/hsps2reciprocal/hsps2reciprocal.py
    ../../ortho_search/hsps2reciprocal/out/*/*.tsv
../genome_stats/genome_stats.py
    ../genome_stats/out/gnid_nums.tsv
./params.tsv
"""