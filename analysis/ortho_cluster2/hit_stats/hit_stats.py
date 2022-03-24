"""Plot various statistics of hits."""

import multiprocessing as mp
import os
from itertools import permutations

import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
import pandas as pd
from numpy import linspace


# Load functions
def load_hsp(qspid, sspid):
    df = pd.read_csv(f'../../ortho_search/blast2hsps/out/hsps/{qspid}/{sspid}.tsv', sep='\t',
                     usecols=hsp_dtypes.keys(), dtype=hsp_dtypes, memory_map=True)

    return df[df['disjoint']]


def load_hit(qspid, sspid):
    df = pd.read_csv(f'../../ortho_search/hsps2hits/out/{qspid}/{sspid}.tsv', sep='\t',
                     usecols=hit_dtypes.keys(), dtype=hit_dtypes, memory_map=True)
    r = pd.read_csv(f'../../ortho_search/hits2reciprocal/out/{qspid}/{sspid}.tsv', sep='\t',
                    usecols=['reciprocal'], memory_map=True)
    max_bitscore = df.groupby('qppid')['bitscore'].max().rename('bitscore_max')
    dfr = df.join(r).join(max_bitscore, on='qppid')

    return dfr


# Plot functions
def hist1(df, bins, file_label, x_label, df_label, color, wrap=False):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel(x_label[0].upper() + x_label[1:])
    plt.ylabel('Number of hits')
    plt.title('Distribution of hits across' + ('\n' if wrap else ' ') + x_label)
    plt.legend()
    plt.savefig(f'out/blast/hist_{file_label}.png')
    plt.close()


def hist3(dfs, bins, file_label, x_label, df_labels, colors, wrap=False):
    fig, axs = plt.subplots(3, 1, figsize=(4.8, 6), sharex=True)
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for ax, df, df_label, color in zip(axs, dfs, df_labels, colors):
        ax.hist(df, bins=b, label=df_label, color=color)
        ax.legend()
    axs[2].set_xlabel(x_label[0].upper() + x_label[1:])
    axs[1].set_ylabel('Number of hits')
    fig.suptitle('Distribution of hits across' + ('\n' if wrap else ' ') + x_label)
    fig.subplots_adjust(left=0.175)
    fig.savefig(f'out/blast/hist_{file_label}_3.png')
    plt.close()


def bar1(df, file_label, x_label, df_label, color, wrap=False):
    plt.bar(df.index, df.values, label=df_label, color=color, width=1)
    plt.xlabel(x_label[0].upper() + x_label[1:])
    plt.ylabel('Number of hits')
    plt.title('Distribution of hits across' + ('\n' if wrap else ' ') + x_label)
    plt.legend()
    plt.savefig(f'out/blast/hist_{file_label}.png')
    plt.close()


def bar3(dfs, file_label, x_label, df_labels, colors, wrap=False):
    fig, axs = plt.subplots(3, 1, figsize=(4.8, 6), sharex=True)
    for ax, df, df_label, color in zip(axs, dfs, df_labels, colors):
        ax.bar(df.index, df.values, label=df_label, color=color, width=1)
        ax.legend()
    axs[2].set_xlabel(x_label[0].upper() + x_label[1:])
    axs[1].set_ylabel('Number of hits')
    fig.suptitle('Distribution of hits across' + ('\n' if wrap else ' ') + x_label)
    fig.subplots_adjust(left=0.175)
    fig.savefig(f'out/blast/hist3_{file_label}.png')
    plt.close()


def bar_hits(counts, file_label):
    plt.bar(counts.keys(), counts.values(), width=1)
    plt.title('Distribution of polypeptides across number of reciprocal hits')
    plt.xlabel('Number of reciprocal hits to polypeptides')
    plt.ylabel('Number of polypeptides')
    plt.savefig(f'out/hits/hist_ppidnum-hitnum_{file_label}.png')
    plt.close()


hsp_dtypes = {'qppid': 'string', 'qgnid': 'string',
              'sppid': 'string', 'sgnid': 'string',
              'bitscore': float,
              'disjoint': bool}
hit_dtypes = {'qppid': 'string', 'qgnid': 'string', 'qspid': 'string',
              'sppid': 'string', 'sgnid': 'string', 'sspid': 'string',
              'hspnum': int, 'chspnum': int,
              'qlen': int, 'nqa': int, 'cnqa': int,
              'slen': int, 'nsa': int, 'cnsa': int,
              'bitscore': float}
num_processes = 2

if __name__ == '__main__':
    # Parse genomes
    spids = []
    with open('../config/genomes.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            spids.append(line.split()[0])

    # Load data
    with mp.Pool(processes=num_processes) as pool:
        hsps0 = pd.concat(pool.starmap(load_hsp, permutations(spids, 2)))
        hsps1 = hsps0[hsps0['bitscore'] >= 50]

        hits0 = pd.concat(pool.starmap(load_hit, permutations(spids, 2)))
        hits0['fqa'] = hits0['nqa'] / hits0['qlen']
        hits0['fsa'] = hits0['nsa'] / hits0['slen']
        hits0['cfqa'] = hits0['cnqa'] / hits0['qlen']
        hits0['xfqa'] = hits0['cfqa'] - hits0['fqa']
        hits0['xhspnum'] = hits0['chspnum'] - hits0['hspnum']

        hits1 = hits0[hits0['bitscore'] == hits0['max_bitscore']]
        hits2 = hits1[hits1['cfqa'] >= 0.5]
        hits3 = hits2[hits2['reciprocal']]
        hits = [hits0, hits1, hits2, hits3]

        labels = ['bitscore ≥ 50', 'best', 'cFQA ≥ 0.5', 'reciprocal']
        colors = ['C1', 'C2', 'C3', 'C4']

    # 1 BLAST METRICS
    # Make BLAST output directory
    if not os.path.exists('out/blast/'):
        os.makedirs('out/blast/')  # Recursive folder creation

    # 1 FILTER PLOTS
    ys1 = [len(hsps0), len(hsps1),
           hits1['hspnum'].sum(), hits2['hspnum'].sum(), hits3['hspnum'].sum()]
    ys2 = [len(hsps0[['qppid', 'sppid']].drop_duplicates()), len(hsps1[['qppid', 'sppid']].drop_duplicates()),
           len(hits1), len(hits2), len(hits3)]
    xs = ['all'] + labels
    cs = ['C0'] + colors

    # 1.1 Number of HSPs in each filter condition
    plt.bar(xs, ys1, color=cs, width=0.25)
    plt.xlim((-0.75, 4.75))
    plt.ylabel('Number of disjoint HSPs')
    plt.savefig('out/blast/bar_hspnum-filter.png')
    plt.close()

    # 1.2 Number of hits in each filter condition
    plt.bar(xs, ys2, color=cs, width=0.25)
    plt.xlim((-0.75, 4.75))
    plt.ylabel('Number of hits')
    plt.savefig('out/blast/bar_hitnum-filter.png')
    plt.close()

    # 1.3 QUALITY METRICS
    # 1.3.1 Bitscore histograms
    hist3([hit['bitscore'] for hit in hits[1:]], 200, 'bitscore', 'bitscore', labels[1:], colors[1:])
    hist1(hits0['bitscore'], 200, 'bitscore_bs50', 'bitscore', labels[0], colors[0])
    hist1(hits1['bitscore'], 200, 'bitscore_best', 'bitscore', labels[1], colors[1])
    hist1(hits2['bitscore'], 200, 'bitscore_cfqa50', 'bitscore', labels[2], colors[2])
    hist1(hits3['bitscore'], 200, 'bitscore_reciprocal', 'bitscore', labels[3], colors[3])

    # 1.3.2 HSPnum histograms
    counts = [hit['hspnum'].value_counts() for hit in hits]
    bar3(counts, 'hspnum', 'number of disjoint HSPs in hit', labels[1:], colors[1:], wrap=True)
    bar1(counts[0], 'hspnum_bs50', 'number of disjoint HSPs in hit', labels[0], colors[0], wrap=True)
    bar1(counts[1], 'hspnum_best', 'number of disjoint HSPs in hit', labels[1], colors[1], wrap=True)
    bar1(counts[2], 'hspnum_cfqa50', 'number of disjoint HSPs in hit', labels[2], colors[2], wrap=True)
    bar1(counts[3], 'hspnum_reciprocal', 'number of disjoint HSPs in hit', labels[3], colors[3], wrap=True)

    # 1.3.3 cHSPnum histograms
    counts = [hit['chspnum'].value_counts() for hit in hits]
    bar3(counts, 'chspnum', 'number of compatible HSPs in hit', labels[1:], colors[1:], wrap=True)
    bar1(counts[0], 'chspnum_bs50', 'number of compatible HSPs in hit', labels[0], colors[0], wrap=True)
    bar1(counts[1], 'chspnum_best', 'number of compatible HSPs in hit', labels[1], colors[1], wrap=True)
    bar1(counts[2], 'chspnum_cfqa50', 'number of compatible HSPs in hit', labels[2], colors[2], wrap=True)
    bar1(counts[3], 'chspnum_reciprocal', 'number of compatible HSPs in hit', labels[3], colors[3], wrap=True)

    # 1.3.4 xHSPnum histograms
    counts = [hit['xhspnum'].value_counts() for hit in hits]
    bar3(counts, 'xhspnum', 'excess number of HSPs in hit', labels[1:], colors[1:], wrap=True)
    bar1(counts[0], 'xhspnum_bs50', 'excess number of HSPs in hit', labels[0], colors[0], wrap=True)
    bar1(counts[1], 'xhspnum_best', 'excess number of HSPs in hit', labels[1], colors[1], wrap=True)
    bar1(counts[2], 'xhspnum_cfqa50', 'excess number of HSPs in hit', labels[2], colors[2], wrap=True)
    bar1(counts[3], 'xhspnum_reciprocal', 'excess number of HSPs in hit', labels[3], colors[3], wrap=True)

    # 1.3.5 NQA histograms
    hist3([hit['nqa'] for hit in hits[1:]], 200, 'nqa', 'number of query aligned', labels[1:], colors[1:], wrap=True)
    hist1(hits0['nqa'], 200, 'nqa_bs50', 'number of query aligned', labels[0], colors[0], wrap=True)
    hist1(hits1['nqa'], 200, 'nqa_best', 'number of query aligned', labels[1], colors[1], wrap=True)
    hist1(hits2['nqa'], 200, 'nqa_cfqa50', 'number of query aligned', labels[2], colors[2], wrap=True)
    hist1(hits3['nqa'], 200, 'nqa_reciprocal', 'number of query aligned', labels[3], colors[3], wrap=True)

    # 1.3.6 cNQA histograms
    hist3([hit['cnqa'] for hit in hits[1:]], 200, 'cnqa', 'compatible number of query aligned', labels[1:], colors[1:], wrap=True)
    hist1(hits0['cnqa'], 200, 'cnqa_bs50', 'compatible number of query aligned', labels[0], colors[0], wrap=True)
    hist1(hits1['cnqa'], 200, 'cnqa_best', 'compatible number of query aligned', labels[1], colors[1], wrap=True)
    hist1(hits2['cnqa'], 200, 'cnqa_cfqa50', 'compatible number of query aligned', labels[2], colors[2], wrap=True)
    hist1(hits3['cnqa'], 200, 'cnqa_reciprocal', 'compatible number of query aligned', labels[3], colors[3], wrap=True)

    # 1.3.7 FQA histograms
    hist3([hit['fqa'] for hit in hits[1:]], 50, 'fqa', 'fraction of query aligned', labels[1:], colors[1:], wrap=True)
    hist1(hits0['fqa'], 50, 'fqa_bs50', 'fraction of query aligned', labels[0], colors[0], wrap=True)
    hist1(hits1['fqa'], 50, 'fqa_best', 'fraction of query aligned', labels[1], colors[1], wrap=True)
    hist1(hits2['fqa'], 50, 'fqa_cfqa50', 'fraction of query aligned', labels[2], colors[2], wrap=True)
    hist1(hits3['fqa'], 50, 'fqa_reciprocal', 'fraction of query aligned', labels[3], colors[3], wrap=True)

    # 1.3.8 cFQA histograms
    hist3([hit['cfqa'] for hit in hits[1:]], 50, 'cfqa', 'compatible fraction of query aligned', labels[1:], colors[1:], wrap=True)
    hist1(hits0['cfqa'], 50, 'cfqa_bs50', 'compatible fraction of query aligned', labels[0], colors[0], wrap=True)
    hist1(hits1['cfqa'], 50, 'cfqa_best', 'compatible fraction of query aligned', labels[1], colors[1], wrap=True)
    hist1(hits2['cfqa'], 50, 'cfqa_cfqa50', 'compatible fraction of query aligned', labels[2], colors[2], wrap=True)
    hist1(hits3['cfqa'], 50, 'cfqa_reciprocal', 'compatible fraction of query aligned', labels[3], colors[3], wrap=True)

    # 1.3.9 xFQA histograms
    hist3([hit['xfqa'] for hit in hits[1:]], 50, 'xfqa', 'excess fraction of query aligned', labels[1:], colors[1:], wrap=True)
    hist1(hits0['xfqa'], 50, 'xfqa_bs50', 'excess fraction of query aligned', labels[0], colors[0], wrap=True)
    hist1(hits1['xfqa'], 50, 'xfqa_best', 'excess fraction of query aligned', labels[1], colors[1], wrap=True)
    hist1(hits2['xfqa'], 50, 'xfqa_cfqa50', 'excess fraction of query aligned', labels[2], colors[2], wrap=True)
    hist1(hits3['xfqa'], 50, 'xfqa_reciprocal', 'excess fraction of query aligned', labels[3], colors[3], wrap=True)

    # 1.3.10 FQA-FSA scatters
    for label, hit in zip(['best', 'bs50', 'cfqa50', 'reciprocal'], hits):
        plt.hist2d(hit['fqa'], hit['fsa'], bins=50, norm=mpl_colors.PowerNorm(0.3))
        plt.xlabel('Fraction of query aligned')
        plt.ylabel('Fraction of subject aligned')
        plt.colorbar()
        plt.savefig(f'out/blast/hist2d_fsa-fqa_{label}.png')
        plt.close()

    # 2 HIT METRICS
    # Make hits output directory
    if not os.path.exists('out/hits/'):
        os.mkdir('out/hits/')

    # 2.1 GNID RANKS
    ids = hits2.loc[:, ['sgnid', 'sspid']].drop_duplicates().set_index('sgnid')

    sgnid_hitnum = hits2.groupby('sgnid')['qgnid'].nunique().rename('sgnid_hitnum').sort_values(ascending=False).to_frame()
    sgnid_hitnum = sgnid_hitnum.join(ids)
    sgnid_hitnum.to_csv('out/hits/sgnids.tsv', sep='\t', index_label='sgnid')

    sgnid_hitnum_dmel = sgnid_hitnum.loc[sgnid_hitnum['sspid'] == 'dmel', :]
    sgnid_hitnum_dmel.to_csv('out/hits/sgnids_dmel.tsv', sep='\t', index_label='sgnid')

    # 2.2 PPID RANKS
    ids = hits2.loc[:, ['sppid', 'sspid']].drop_duplicates().set_index('sppid')

    sppid_hitnum = hits2.groupby('sppid')['qppid'].nunique().rename('sppid_hitnum').sort_values(ascending=False).to_frame()
    sppid_hitnum = sppid_hitnum.join(ids)
    sppid_hitnum.to_csv('out/hits/sppids.tsv', sep='\t', index_label='sppid')

    sppid_hitnum_dmel = sppid_hitnum.loc[sppid_hitnum['sspid'] == 'dmel', :]
    sppid_hitnum_dmel.to_csv('out/hits/sppids_dmel.tsv', sep='\t', index_label='sppid')

    # 2.3 PLOTS
    # 2.3.1 Correlation of gene hits with number of associated polypeptides
    gnid_ppidnum = pd.read_csv('../genome_stats/out/gnid_nums.tsv', sep='\t', index_col='gnid', dtype={'gnid': 'string'})
    corr = sgnid_hitnum.join(gnid_ppidnum)

    plt.scatter(corr['ppidnum'], corr['sgnid_hitnum'], alpha=0.5, s=10, edgecolors='none')
    plt.xlabel('Number of polypeptides associated with gene')
    plt.ylabel('Number of reciprocal hits to gene')
    plt.title('Correlation of number of reciprocal hits to gene\nwith number of associated polypeptides')
    plt.savefig('out/hits/scatter_hitnum-ppidnum.png')
    plt.close()

    # 2.3.2 Histograms of polypeptides by number of hits
    counts = sppid_hitnum['sppid_hitnum'].value_counts().to_dict()

    bar_hits(counts, 'all')
    bar_hits({key: value for key, value in counts.items() if key > len(spids)}, f'{len(spids)}+')
    bar_hits({key: value for key, value in counts.items() if key <= len(spids)}, f'{len(spids)}-')

"""
DEPENDENCIES
../../ortho_search/blast2hsps/blast2hsps.py
    ../../ortho_search/blast2hsps/out/hsps/*/*.tsv
../../ortho_search/hsps2hits/hsps2hits.py
    ../../ortho_search/hsps2hits/out/*/*.tsv
../../ortho_search/hits2reciprocal/hits2reciprocal.py
    ../../ortho_search/hits2reciprocal/out/*/*.tsv
../config/genomes.tsv
../genome_stats/genome_stats.py
    ../genome_stats/out/gnid_nums.tsv
"""