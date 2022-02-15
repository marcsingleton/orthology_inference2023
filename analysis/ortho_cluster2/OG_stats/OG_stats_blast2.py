"""Plot various statistics of OGs relating to their BLAST parameters."""

import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import pandas as pd
from itertools import permutations
from numpy import linspace


# Load functions
def load_hit(qspid, sspid):
    df = pd.read_csv(f'../../ortho_search/hsps2hits/out/{qspid}/{sspid}.tsv', sep='\t',
                     usecols=hit_dtypes.keys(), dtype=hit_dtypes, memory_map=True)
    r = pd.read_csv(f'../../ortho_search/hits2reciprocal/out/{qspid}/{sspid}.tsv', sep='\t',
                    usecols=['reciprocal2'], memory_map=True)

    return df[r['reciprocal2']]


# Plot functions
def hist1(df, bins, file_label, title_label, x_label, df_label, color, wrap=False):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel(x_label[0].upper() + x_label[1:])
    plt.ylabel(f'Number of {title_label}')
    plt.title(f'Distribution of {title_label} across' + ('\n' if wrap else ' ') + x_label)
    plt.legend()
    plt.savefig(f'out/pgraph2/blast/hist_{file_label}.png')
    plt.close()


def hist3(dfs, bins, file_label, title_label, x_label, df_labels, colors, wrap=False):
    fig, axs = plt.subplots(3, 1, figsize=(4.8, 6), sharex=True)
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for ax, df, data_label, color in zip(axs, dfs, df_labels, colors):
        ax.hist(df, bins=b, label=data_label, color=color)
        ax.legend()
    axs[2].set_xlabel(x_label[0].upper() + x_label[1:])
    axs[1].set_ylabel(f'Number of {title_label}')
    fig.suptitle(f'Distribution of {title_label} across' + ('\n' if wrap else ' ') + x_label)
    fig.subplots_adjust(left=0.175)
    fig.savefig(f'out/pgraph2/blast/hist3_{file_label}.png')
    plt.close()


def bar1(df, file_label, x_label, df_label, color, wrap=False):
    plt.bar(df.index, df.values, label=df_label, color=color, width=1)
    plt.xlabel(x_label[0].upper() + x_label[1:])
    plt.ylabel('Number of hits')
    plt.title('Distribution of hits across' + ('\n' if wrap else ' ') + x_label)
    plt.legend()
    plt.savefig(f'out/pgraph2/blast/hist_{file_label}.png')
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
    fig.savefig(f'out/pgraph2/blast/hist_{file_label}.png')
    plt.close()


def scatter1(x, y, file_label, xy_label, df_label, color):
    fig, ax = plt.subplots()
    ax.scatter(x, y, alpha=0.5, s=10, label=df_label, color=color, edgecolors='none')
    ax.set_xlabel(f'Mean {xy_label} in OG')
    ax.set_ylabel(f'Variance of {xy_label} in OG')
    leg = ax.legend(markerscale=2)
    for lh in leg.legendHandles:
        lh.set_alpha(1)
    fig.savefig(f'out/pgraph2/blast/scatter_{file_label}.png')
    plt.close()


def scatter2(x, y, file_label, y_label):
    plt.scatter(x, y, alpha=0.5, s=10, label='all', color='C0', edgecolors='none')
    plt.xlabel('Number of polypeptides in OG')
    plt.ylabel(y_label)
    plt.savefig(f'out/pgraph2/blast/scatter_{file_label}.png')
    plt.close()


hit_dtypes = {'qppid': 'string', 'qgnid': 'string', 'qspid': 'string',
              'sppid': 'string', 'sgnid': 'string', 'sspid': 'string',
              'hspnum': int, 'chspnum': int,
              'qlen': int, 'nqa': int, 'cnqa': int,
              'slen': int, 'nsa': int, 'cnsa': int,
              'bitscore': float}
num_processes = 4

if __name__ == '__main__':
    # Parse genomes
    spids = []
    with open('../config/genomes.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            spids.append(line.split()[0])

    # Load data
    rows = []
    with open('../subcluster_pgraph/out/pgraph2/pclusters.txt') as file:
        for line in file:
            CCid, OGid, edges = line.rstrip().split(':')
            for edge in edges.split('\t'):
                node1, node2 = edge.split(',')
                rows.append({'CCid': CCid, 'OGid': OGid, 'qppid': node1, 'sppid': node2})
                rows.append({'CCid': CCid, 'OGid': OGid, 'qppid': node2, 'sppid': node1})
    edges = pd.DataFrame(rows)

    with mp.Pool(processes=num_processes) as pool:
        hits0 = pd.concat(pool.starmap(load_hit, permutations(spids, 2)))
        hits0 = edges.merge(hits0, how='left', on=['qppid', 'sppid'])
        hits0['fqa'] = hits0['nqa'] / hits0['qlen']
        hits0['fsa'] = hits0['nsa'] / hits0['slen']
        hits0['cfqa'] = hits0['cnqa'] / hits0['qlen']
        hits0['xfqa'] = hits0['cfqa'] - hits0['fqa']
        hits0['xhspnum'] = hits0['chspnum'] - hits0['hspnum']

    # Segment OGs
    OGs = hits0.groupby('OGid')
    gns = hits0.groupby('qgnid')

    OG_spidnum = OGs['qspid'].nunique()
    OG_gnidnum = OGs['qgnid'].nunique()
    gn_OGidnum = gns['OGid'].nunique()

    OGs_sps = OG_spidnum[OG_spidnum == len(spids)].index
    OGs_gns = OG_gnidnum[OG_gnidnum == len(spids)].index
    gns_gn1OG = gn_OGidnum[gn_OGidnum > 1].index
    OGs_gn1OG = hits0.loc[hits0['qgnid'].isin(gns_gn1OG), 'OGid'].unique()

    OGs1 = set(OGs_sps) & set(OGs_gns)
    OGs2 = OGs1 - set(OGs_gn1OG)

    hits1 = hits0[hits0['OGid'].isin(OGs1)]
    hits2 = hits0[hits0['OGid'].isin(OGs2)]
    hits = [hits0, hits1, hits2]
    OGs = [hit.groupby('OGid') for hit in hits]

    labels = ['all', 'filter1', 'filter2']
    colors = ['C0', 'C1', 'C2']

    # Make output directory
    if not os.path.exists('out/pgraph2/blast/'):
        os.makedirs('out/pgraph2/blast/')

    # 1 FILTER PLOT
    ys = [hit['OGid'].nunique() for hit in hits]
    plt.bar(labels, ys, color=colors, width=0.25)
    plt.xlim((-0.75, 2.75))
    plt.ylabel('Number of OGs')
    plt.savefig('out/pgraph2/blast/bar_OGnum-filter.png')
    plt.close()

    # 2 BITSCORE PLOTS
    # 2.1 Hit histograms
    hist3([hit['bitscore'] for hit in hits], 200, 'hitnum-bitscore', 'hits in OGs', 'bitscore', labels, colors)
    hist1(hits0['bitscore'], 200, 'hitnum-bitscore_all', 'hits in OGs', 'bitscore', labels[0], colors[0])
    hist1(hits1['bitscore'], 200, 'hitnum-bitscore_filter1', 'hits in OGs', 'bitscore', labels[1], colors[1])
    hist1(hits2['bitscore'], 200, 'hitnum-bitscore_filter2', 'hits in OGs', 'bitscore', labels[2], colors[2])

    # 2.2 OG histograms
    hist3([OG['bitscore'].mean() for OG in OGs], 200, 'OGnum-bitscoremean', 'OGs',
          'mean bitscore of hits in OG', labels, colors, wrap=True)
    hist3([OG['bitscore'].var() for OG in OGs], 200, 'OGnum-bitscorevar', 'OGs',
          'variance of bitscore of hits in OG', labels, colors, wrap=True)

    # 2.3 OG scatters
    scatter1(OGs[0]['bitscore'].mean(), OGs[0]['bitscore'].var(),
             'bitscorevar-bitscoremean_all.png', 'bitscore', labels[0], colors[0])
    scatter1(OGs[1]['bitscore'].mean(), OGs[1]['bitscore'].var(),
             'bitscorevar-bitscoremean_filter1.png', 'bitscore', labels[1], colors[1])
    scatter1(OGs[2]['bitscore'].mean(), OGs[2]['bitscore'].var(),
             'bitscorevar-bitscoremean_filter2.png', 'bitscore', labels[2], colors[2])

    # 3 HSPNUM HISTOGRAMS
    counts = [hit['hspnum'].value_counts() for hit in hits]
    bar3(counts, 'hitnum-hspnum', 'number of disjoint HSPs in hit', labels, colors, wrap=True)
    bar1(counts[0], 'hitnum-hspnum_all', 'number of disjoint HSPs in hit', labels[0], colors[0], wrap=True)
    bar1(counts[1], 'hitnum-hspnum_filter1', 'number of disjoint HSPs in hit', labels[1], colors[1], wrap=True)
    bar1(counts[2], 'hitnum-hspnum_filter2', 'number of disjoint HSPs in hit', labels[2], colors[2], wrap=True)

    # 4 cHSPNUM HISTOGRAMS
    counts = [hit['chspnum'].value_counts() for hit in hits]
    bar3(counts, 'hitnum-chspnum', 'number of compatible HSPs in hit', labels, colors, wrap=True)
    bar1(counts[0], 'hitnum-chspnum_all', 'number of compatible HSPs in hit', labels[0], colors[0], wrap=True)
    bar1(counts[1], 'hitnum-chspnum_filter1', 'number of compatible HSPs in hit', labels[1], colors[1], wrap=True)
    bar1(counts[2], 'hitnum-chspnum_filter2', 'number of compatible HSPs in hit', labels[2], colors[2], wrap=True)

    # 5 xHSPNUM HISTOGRAMS
    counts = [hit['xhspnum'].value_counts() for hit in hits]
    bar3(counts, 'hitnum-xhspnum', 'excess number of HSPs in hit', labels, colors, wrap=True)
    bar1(counts[0], 'hitnum-xhspnum_all', 'excess number of HSPs in hit', labels[0], colors[0], wrap=True)
    bar1(counts[1], 'hitnum-xhspnum_filter1', 'excess number of HSPs in hit', labels[1], colors[1], wrap=True)
    bar1(counts[2], 'hitnum-xhspnum_filter2', 'excess number of HSPs in hit', labels[2], colors[2], wrap=True)

    # 6 NQA PLOTS
    # 6.1 Hit histograms
    hist3([hit['nqa'] for hit in hits], 200, 'hitnum-nqa', 'hits in OGs', 'number of query aligned',
          labels, colors, wrap=True)
    hist1(hits0['nqa'], 200, 'hitnum-nqa_all', 'hits in OGs', 'number of query aligned',
          labels[0], colors[0], wrap=True)
    hist1(hits1['nqa'], 200, 'hitnum-nqa_filter1', 'hits in OGs', 'number of query aligned',
          labels[1], colors[1], wrap=True)
    hist1(hits2['nqa'], 200, 'hitnum-nqa_filter2', 'hits in OGs', 'number of query aligned',
          labels[2], colors[2], wrap=True)

    # 6.2 OG histograms
    hist3([OG['nqa'].mean() for OG in OGs], 50, 'OGnum-nqamean', 'OGs', 'mean NQA of hits in OG',
          labels, colors, wrap=True)
    hist3([OG['nqa'].var() for OG in OGs], 50, 'OGnum-nqavar', 'OGs', 'variance of NQA of hits in OG',
          labels, colors, wrap=True)

    # 6.3 OG scatters
    scatter1(OGs[0]['nqa'].mean(), OGs[0]['nqa'].var(), 'nqavar-nqamean_all.png', 'NQA', labels[0], colors[0])
    scatter1(OGs[1]['nqa'].mean(), OGs[1]['nqa'].var(), 'nqavar-nqamean_filter1.png', 'NQA', labels[1], colors[1])
    scatter1(OGs[2]['nqa'].mean(), OGs[2]['nqa'].var(), 'nqavar-nqamean_filter2.png', 'NQA', labels[2], colors[2])

    # 7 cNQA PLOTS
    # 7.1 Hit histograms
    hist3([hit['cnqa'] for hit in hits], 200, 'hitnum-cnqa', 'hits in OGs', 'compatible number of query aligned',
          labels, colors, wrap=True)
    hist1(hits0['cnqa'], 200, 'hitnum-cnqa_all', 'hits in OGs', 'compatible number of query aligned',
          labels[0], colors[0], wrap=True)
    hist1(hits1['cnqa'], 200, 'hitnum-cnqa_filter1', 'hits in OGs', 'compatible number of query aligned',
          labels[1], colors[1], wrap=True)
    hist1(hits2['cnqa'], 200, 'hitnum-cnqa_filter2', 'hits in OGs', 'compatible number of query aligned',
          labels[2], colors[2], wrap=True)

    # 7.2 OG histograms
    hist3([OG['cnqa'].mean() for OG in OGs], 50, 'OGnum-cnqamean', 'OGs', 'mean cNQA of hits in OG',
          labels, colors, wrap=True)
    hist3([OG['cnqa'].var() for OG in OGs], 50, 'OGnum-cnqavar', 'OGs', 'variance of cNQA of hits in OG',
          labels, colors, wrap=True)

    # 7.3 OG scatters
    scatter1(OGs[0]['cnqa'].mean(), OGs[0]['cnqa'].var(), 'cnqavar-cnqamean_all.png', 'cNQA', labels[0], colors[0])
    scatter1(OGs[1]['cnqa'].mean(), OGs[1]['cnqa'].var(), 'cnqavar-cnqamean_filter1.png', 'cNQA', labels[1], colors[1])
    scatter1(OGs[2]['cnqa'].mean(), OGs[2]['cnqa'].var(), 'cnqavar-cnqamean_filter2.png', 'cNQA', labels[2], colors[2])

    # 8 FQA PLOTS
    # 8.1 Hit histograms
    hist3([hit['fqa'] for hit in hits], 50, 'hitnum-fqa', 'hits in OGs', 'fraction of query aligned',
          labels, colors, wrap=True)
    hist1(hits0['fqa'], 50, 'hitnum-fqa_all', 'hits in OGs', 'fraction of query aligned',
          labels[0], colors[0], wrap=True)
    hist1(hits1['fqa'], 50, 'hitnum-fqa_filter1', 'hits in OGs', 'fraction of query aligned',
          labels[1], colors[1], wrap=True)
    hist1(hits2['fqa'], 50, 'hitnum-fqa_filter2', 'hits in OGs', 'fraction of query aligned',
          labels[2], colors[2], wrap=True)

    # 8.2 OG histograms
    hist3([OG['fqa'].mean() for OG in OGs], 50, 'OGnum-fqamean', 'OGs',
          'mean FQA of hits in OG', labels, colors, wrap=True)
    hist3([OG['fqa'].var() for OG in OGs], 50, 'OGnum-fqavar', 'OGs',
          'variance of FQA of hits in OG', labels, colors, wrap=True)

    # 8.3 OG scatters
    scatter1(OGs[0]['fqa'].mean(), OGs[0]['fqa'].var(), 'fqavar-fqamean_all.png', 'FQA', labels[0], colors[0])
    scatter1(OGs[1]['fqa'].mean(), OGs[1]['fqa'].var(), 'fqavar-fqamean_filter1.png', 'FQA', labels[1], colors[1])
    scatter1(OGs[2]['fqa'].mean(), OGs[2]['fqa'].var(), 'fqavar-fqamean_filter2.png', 'FQA', labels[2], colors[2])

    # 9 cFQA PLOTS
    # 9.1 Hit histograms
    hist3([hit['cfqa'] for hit in hits], 50, 'hitnum-cfqa', 'hits in OGs', 'compatible fraction of query aligned',
          labels, colors, wrap=True)
    hist1(hits0['cfqa'], 50, 'hitnum-cfqa_all', 'hits in OGs', 'compatible fraction of query aligned',
          labels[0], colors[0], wrap=True)
    hist1(hits1['cfqa'], 50, 'hitnum-cfqa_filter1', 'hits in OGs', 'compatible fraction of query aligned',
          labels[1], colors[1], wrap=True)
    hist1(hits2['cfqa'], 50, 'hitnum-cfqa_filter2', 'hits in OGs', 'compatible fraction of query aligned',
          labels[2], colors[2], wrap=True)

    # 9.2 OG histograms
    hist3([OG['cfqa'].mean() for OG in OGs], 50, 'OGnum-cfqamean', 'OGs',
          'mean cFQA of hits in OG', labels, colors, wrap=True)
    hist3([OG['cfqa'].var() for OG in OGs], 50, 'OGnum-cfqavar', 'OGs',
          'variance of cFQA of hits in OG', labels, colors, wrap=True)

    # 9.3 OG scatters
    scatter1(OGs[0]['cfqa'].mean(), OGs[0]['cfqa'].var(), 'cfqavar-cfqamean_all.png', 'cFQA', labels[0], colors[0])
    scatter1(OGs[1]['cfqa'].mean(), OGs[1]['cfqa'].var(), 'cfqavar-cfqamean_filter1.png', 'cFQA', labels[1], colors[1])
    scatter1(OGs[2]['cfqa'].mean(), OGs[2]['cfqa'].var(), 'cfqavar-cfqamean_filter2.png', 'cFQA', labels[2], colors[2])

    # 10 xFQA PLOTS
    # 10.1 Hit histograms
    hist3([hit['xfqa'] for hit in hits], 50, 'hitnum-xfqa', 'hits in OGs', 'excess fraction of query aligned',
          labels, colors, wrap=True)
    hist1(hits0['xfqa'], 50, 'hitnum-xfqa_all', 'hits in OGs', 'excess fraction of query aligned',
          labels[0], colors[0], wrap=True)
    hist1(hits1['xfqa'], 50, 'hitnum-xfqa_filter1', 'hits in OGs', 'excess fraction of query aligned',
          labels[1], colors[1], wrap=True)
    hist1(hits2['xfqa'], 50, 'hitnum-xfqa_filter2', 'hits in OGs', 'excess fraction of query aligned',
          labels[2], colors[2], wrap=True)

    # 10.2 OG histograms
    hist3([OG['xfqa'].mean() for OG in OGs], 50, 'OGnum-xfqamean', 'OGs',
          'mean xFQA of hits in OG', labels, colors, wrap=True)
    hist3([OG['xfqa'].var() for OG in OGs], 50, 'OGnum-xfqavar', 'OGs',
          'variance of xFQA of hits in OG', labels, colors, wrap=True)

    # 10.3 OG scatters
    scatter1(OGs[0]['xfqa'].mean(), OGs[0]['xfqa'].var(), 'xfqavar-xfqamean_all.png', 'xFQA', labels[0], colors[0])
    scatter1(OGs[1]['xfqa'].mean(), OGs[1]['xfqa'].var(), 'xfqavar-xfqamean_filter1.png', 'xFQA', labels[1], colors[1])
    scatter1(OGs[2]['xfqa'].mean(), OGs[2]['xfqa'].var(), 'xfqavar-xfqamean_filter2.png', 'xFQA', labels[2], colors[2])

    # 11 FQA-FSA SCATTERS
    for label, hit in zip(labels, hits):
        plt.hist2d(hit['fqa'], hit['fsa'], bins=50, norm=mpl_colors.PowerNorm(0.3))
        plt.xlabel('Fraction of query aligned')
        plt.ylabel('Fraction of subject aligned')
        plt.colorbar()
        plt.savefig(f'out/pgraph2/blast/hist2d_fsa-fqa_{label}.png')
        plt.close()

    # 12 EDGES
    edgenums = [hit[['qppid', 'sppid', 'OGid']].drop_duplicates().groupby('OGid').size() / 2 for hit in hits]
    ppidnums = [OG['qppid'].nunique() for OG in OGs]
    edgefracs = [2 * edgenum / (gnidnum*(gnidnum-1)) for edgenum, gnidnum in zip(edgenums, ppidnums)]

    # 12.1 Edge number histograms
    hist3(edgenums, 100, 'OGnum-edgenum', 'OGs', 'number of edges', labels, colors)
    hist1(edgenums[0], 100, 'OGnum-edgenum_all', 'OGs', 'number of edges', labels[0], colors[0])
    hist1(edgenums[1], 50, 'OGnum-edgenum_filter1', 'OGs', 'number of edges', labels[1], colors[1])
    hist1(edgenums[2], 50, 'OGnum-edgenum_filter2', 'OGs', 'number of edges', labels[2], colors[2])

    # 12.2 Edge fraction histograms
    hist3(edgefracs, 50, 'OGnum-edgefrac', 'OGs', 'fraction of possible edges', labels, colors)
    hist1(edgefracs[0], 50, 'OGnum-edgefrac_all', 'OGs', 'fraction of possible edges', labels[0], colors[0])
    hist1(edgefracs[1], 50, 'OGnum-edgefrac_filter1', 'OGs', 'fraction of possible edges', labels[1], colors[1])
    hist1(edgefracs[2], 50, 'OGnum-edgefrac_filter2', 'OGs', 'fraction of possible edges', labels[2], colors[2])

    # 13 CORRELATIONS
    scatter2(ppidnums[0], OGs[0]['bitscore'].mean(), 'bitscore-OGppnum', 'Mean bitscore of hits in OG')
    scatter2(ppidnums[0], OGs[0]['fqa'].mean(), 'fqa-OGppnum', 'Mean fraction of query aligned of hits in OG')
    scatter2(ppidnums[0], edgenums[0], 'edgenum-OGppnum', 'Number of edges in OG')
    scatter2(ppidnums[0], edgefracs[0], 'edgefrac-OGppnum', 'Fraction of possible edges in OG')

"""
DEPENDENCIES
../../ortho_search/hsps2hits/hsps2hits.py
    ../../ortho_search/hsps2hits/out/*/*.tsv
../../ortho_search/hits2reciprocal/hits2reciprocal.py
    ../../ortho_search/hits2reciprocal/out/*/*.tsv
../config/genomes.tsv
../subcluster_pgraph/subcluster_pgraph2.py
    ../subcluster_pgraph/out/pgraph2/pclusters.txt
"""