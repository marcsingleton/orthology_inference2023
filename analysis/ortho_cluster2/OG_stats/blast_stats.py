"""Plot various statistics of OGs relating to their BLAST output."""

import multiprocessing as mp
import os
from itertools import permutations

import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt
import pandas as pd
from numpy import linspace


# Load functions
def load_hit(qspid, sspid):
    df = pd.read_csv(f'../../ortho_search/hsps2hits/out/{qspid}/{sspid}.tsv', sep='\t',
                     usecols=hit_dtypes.keys(), dtype=hit_dtypes, memory_map=True)
    r = pd.read_csv(f'../../ortho_search/hits2reciprocal/out/{qspid}/{sspid}.tsv', sep='\t',
                    usecols=['reciprocal'], memory_map=True)

    return df[r['reciprocal']]


# Plot functions
def hist1(df, bins, file_label, x_label, y_label, df_label, color):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel(x_label[0].upper() + x_label[1:])
    plt.ylabel(f'Number of {y_label}')
    plt.legend()
    plt.savefig(f'out/blast/hist_{file_label}.png')
    plt.close()


def hist3(dfs, bins, file_label, x_label, y_label, df_labels, colors):
    fig, axs = plt.subplots(3, 1, figsize=(4.8, 6), sharex=True)
    b = linspace(min([df.min() for df in dfs]), max([df.max() for df in dfs]), bins, endpoint=True)
    for ax, df, data_label, color in zip(axs, dfs, df_labels, colors):
        ax.hist(df, bins=b, label=data_label, color=color)
        ax.legend()
    axs[2].set_xlabel(x_label[0].upper() + x_label[1:])
    axs[1].set_ylabel(f'Number of {y_label}')
    fig.subplots_adjust(left=0.175)
    fig.savefig(f'out/blast/hist_{file_label}_3.png')
    plt.close()


def bar1(df, file_label, x_label, y_label, df_label, color):
    plt.bar(df.index, df.values, label=df_label, color=color, width=1)
    plt.xlabel(x_label[0].upper() + x_label[1:])
    plt.ylabel(f'Number of {y_label}')
    plt.legend()
    plt.savefig(f'out/blast/hist_{file_label}.png')
    plt.close()


def bar3(dfs, file_label, x_label, y_label, df_labels, colors):
    fig, axs = plt.subplots(3, 1, figsize=(4.8, 6), sharex=True)
    for ax, df, df_label, color in zip(axs, dfs, df_labels, colors):
        ax.bar(df.index, df.values, label=df_label, color=color, width=1)
        ax.legend()
    axs[2].set_xlabel(x_label[0].upper() + x_label[1:])
    axs[1].set_ylabel(f'Number of {y_label}')
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
    plt.xlabel('Number of proteins in OG')
    plt.ylabel(y_label)
    plt.savefig(f'out/blast/scatter_{file_label}.png')
    plt.close()


hit_dtypes = {'qppid': 'string', 'qgnid': 'string', 'qspid': 'string',
              'sppid': 'string', 'sgnid': 'string', 'sspid': 'string',
              'hspnum': int, 'chspnum': int,
              'qlen': int, 'nqa': int, 'cnqa': int,
              'slen': int, 'nsa': int, 'cnsa': int,
              'bitscore': float}
num_processes = 4

if __name__ == '__main__':
    # Load genomes
    spids = []
    with open('../config/genomes.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            spid, _, _, _ = line.split()
            spids.append(spid)

    # Load data
    rows = []
    with open('../cluster4+_graph/out/4clique/clusters.tsv') as file:
        file.readline()  # Skip header
        for line in file:
            component_id, OGid, _, edges = line.rstrip().split('\t')
            for edge in edges.split(','):
                node1, node2 = edge.split(':')
                rows.append({'component_id': component_id, 'OGid': OGid, 'qppid': node1, 'sppid': node2})
                rows.append({'component_id': component_id, 'OGid': OGid, 'qppid': node2, 'sppid': node1})
    edges = pd.DataFrame(rows)

    with mp.Pool(processes=num_processes) as pool:
        hits0 = pd.concat(pool.starmap(load_hit, permutations(spids[:3], 2)))
        hits0 = edges.merge(hits0, how='left', on=['qppid', 'sppid'])
        hits0['fqa'] = hits0['nqa'] / hits0['qlen']
        hits0['fsa'] = hits0['nsa'] / hits0['slen']
        hits0['cfqa'] = hits0['cnqa'] / hits0['qlen']
        hits0['xfqa'] = hits0['cfqa'] - hits0['fqa']
        hits0['xhspnum'] = hits0['chspnum'] - hits0['hspnum']

    # Segment OGs
    OGid_groups = hits0.groupby('OGid')
    gnid_groups = hits0.groupby('qgnid')
    spidnum = OGid_groups['qspid'].nunique()
    gnidnum = OGid_groups['qgnid'].nunique()

    species_filter = spidnum[spidnum == len(spids)].index
    gene_filter = gnidnum[gnidnum == len(spids)].index
    bijective_filter = gnid_groups.filter(lambda x: len(x) == 1)['OGid']

    OGs1 = set(species_filter) & set(gene_filter)
    OGs2 = OGs1 & set(bijective_filter)

    hits1 = hits0[hits0['OGid'].isin(OGs1)]
    hits2 = hits0[hits0['OGid'].isin(OGs2)]
    hits = [hits0, hits1, hits2]
    OGs = [hit.groupby('OGid') for hit in hits]

    labels = ['all', 'filter1', 'filter2']
    colors = ['C0', 'C1', 'C2']

    if not os.path.exists('out/blast/'):
        os.makedirs('out/blast/')

    # 1 FILTER PLOT
    ys = [hit['OGid'].nunique() for hit in hits]
    plt.bar(labels, ys, color=colors, width=0.25)
    plt.xlim((-0.75, 2.75))
    plt.ylabel('Number of OGs')
    plt.savefig('out/blast/bar_OGidnum-filter.png')
    plt.close()

    # 2 BITSCORE PLOTS
    # 2.1 Hit histograms
    hist3([hit['bitscore'] for hit in hits], 200, 'hitnum-bitscore', 'bitscore', 'hits in OGs', labels, colors)
    hist1(hits0['bitscore'], 200, 'hitnum-bitscore_all', 'bitscore', 'hits in OGs', labels[0], colors[0])
    hist1(hits1['bitscore'], 200, 'hitnum-bitscore_filter1', 'bitscore', 'hits in OGs', labels[1], colors[1])
    hist1(hits2['bitscore'], 200, 'hitnum-bitscore_filter2', 'bitscore', 'hits in OGs', labels[2], colors[2])

    # 2.2 OG histograms
    hist3([OG['bitscore'].mean() for OG in OGs], 200, 'OGidnum-bitscoremean', 'mean bitscore of hits in OG', 'OGs', labels, colors)
    hist3([OG['bitscore'].var() for OG in OGs], 200, 'OGidnum-bitscorevar', 'variance of bitscore of hits in OG', 'OGs', labels, colors)

    # 2.3 OG scatters
    scatter1(OGs[0]['bitscore'].mean(), OGs[0]['bitscore'].var(), 'bitscorevar-bitscoremean_all.png', 'bitscore', labels[0], colors[0])
    scatter1(OGs[1]['bitscore'].mean(), OGs[1]['bitscore'].var(), 'bitscorevar-bitscoremean_filter1.png', 'bitscore', labels[1], colors[1])
    scatter1(OGs[2]['bitscore'].mean(), OGs[2]['bitscore'].var(), 'bitscorevar-bitscoremean_filter2.png', 'bitscore', labels[2], colors[2])

    # 3 HSPNUM HISTOGRAMS
    counts = [hit['hspnum'].value_counts() for hit in hits]
    bar3(counts, 'hitnum-hspnum', 'number of disjoint HSPs in hit', 'hits', labels, colors)
    bar1(counts[0], 'hitnum-hspnum_all', 'number of disjoint HSPs in hit', 'hits', labels[0], colors[0])
    bar1(counts[1], 'hitnum-hspnum_filter1', 'number of disjoint HSPs in hit', 'hits', labels[1], colors[1])
    bar1(counts[2], 'hitnum-hspnum_filter2', 'number of disjoint HSPs in hit', 'hits', labels[2], colors[2])

    # 4 cHSPNUM HISTOGRAMS
    counts = [hit['chspnum'].value_counts() for hit in hits]
    bar3(counts, 'hitnum-chspnum', 'number of compatible HSPs in hit', 'hits', labels, colors)
    bar1(counts[0], 'hitnum-chspnum_all', 'number of compatible HSPs in hit', 'hits', labels[0], colors[0])
    bar1(counts[1], 'hitnum-chspnum_filter1', 'number of compatible HSPs in hit', 'hits', labels[1], colors[1])
    bar1(counts[2], 'hitnum-chspnum_filter2', 'number of compatible HSPs in hit', 'hits', labels[2], colors[2])

    # 5 xHSPNUM HISTOGRAMS
    counts = [hit['xhspnum'].value_counts() for hit in hits]
    bar3(counts, 'hitnum-xhspnum', 'excess number of HSPs in hit', 'hits', labels, colors)
    bar1(counts[0], 'hitnum-xhspnum_all', 'excess number of HSPs in hit', 'hits', labels[0], colors[0])
    bar1(counts[1], 'hitnum-xhspnum_filter1', 'excess number of HSPs in hit', 'hits', labels[1], colors[1])
    bar1(counts[2], 'hitnum-xhspnum_filter2', 'excess number of HSPs in hit', 'hits', labels[2], colors[2])

    # 6 NQA PLOTS
    # 6.1 Hit histograms
    hist3([hit['nqa'] for hit in hits], 200, 'hitnum-nqa', 'number of query aligned', 'hits in OGs', labels, colors)
    hist1(hits0['nqa'], 200, 'hitnum-nqa_all', 'number of query aligned', 'hits in OGs', labels[0], colors[0])
    hist1(hits1['nqa'], 200, 'hitnum-nqa_filter1', 'number of query aligned', 'hits in OGs', labels[1], colors[1])
    hist1(hits2['nqa'], 200, 'hitnum-nqa_filter2', 'number of query aligned', 'hits in OGs', labels[2], colors[2])

    # 6.2 OG histograms
    hist3([OG['nqa'].mean() for OG in OGs], 50, 'OGidnum-nqamean', 'mean NQA of hits in OG', 'OGs', labels, colors)
    hist3([OG['nqa'].var() for OG in OGs], 50, 'OGidnum-nqavar', 'variance of NQA of hits in OG', 'OGs', labels, colors)

    # 6.3 OG scatters
    scatter1(OGs[0]['nqa'].mean(), OGs[0]['nqa'].var(), 'nqavar-nqamean_all.png', 'NQA', labels[0], colors[0])
    scatter1(OGs[1]['nqa'].mean(), OGs[1]['nqa'].var(), 'nqavar-nqamean_filter1.png', 'NQA', labels[1], colors[1])
    scatter1(OGs[2]['nqa'].mean(), OGs[2]['nqa'].var(), 'nqavar-nqamean_filter2.png', 'NQA', labels[2], colors[2])

    # 7 cNQA PLOTS
    # 7.1 Hit histograms
    bar3([hit['cnqa'] for hit in hits], 'hitnum-cnqa', 'compatible number of query aligned', 'hits in OGs', labels, colors)
    bar1(hits0['cnqa'], 'hitnum-cnqa_all', 'compatible number of query aligned', 'hits in OGs', labels[0], colors[0])
    bar1(hits1['cnqa'], 'hitnum-cnqa_filter1', 'compatible number of query aligned', 'hits in OGs', labels[1], colors[1])
    bar1(hits2['cnqa'], 'hitnum-cnqa_filter2', 'compatible number of query aligned', 'hits in OGs', labels[2], colors[2])

    # 7.2 OG histograms
    hist3([OG['cnqa'].mean() for OG in OGs], 50, 'OGidnum-cnqamean', 'mean cNQA of hits in OG', 'OGs', labels, colors)
    hist3([OG['cnqa'].var() for OG in OGs], 50, 'OGidnum-cnqavar', 'variance of cNQA of hits in OG', 'OGs', labels, colors)

    # 7.3 OG scatters
    scatter1(OGs[0]['cnqa'].mean(), OGs[0]['cnqa'].var(), 'cnqavar-cnqamean_all.png', 'cNQA', labels[0], colors[0])
    scatter1(OGs[1]['cnqa'].mean(), OGs[1]['cnqa'].var(), 'cnqavar-cnqamean_filter1.png', 'cNQA', labels[1], colors[1])
    scatter1(OGs[2]['cnqa'].mean(), OGs[2]['cnqa'].var(), 'cnqavar-cnqamean_filter2.png', 'cNQA', labels[2], colors[2])

    # 8 FQA PLOTS
    # 8.1 Hit histograms
    hist3([hit['fqa'] for hit in hits], 50, 'hitnum-fqa', 'fraction of query aligned', 'hits in OGs', labels, colors)
    hist1(hits0['fqa'], 50, 'hitnum-fqa_all', 'fraction of query aligned', 'hits in OGs', labels[0], colors[0])
    hist1(hits1['fqa'], 50, 'hitnum-fqa_filter1', 'fraction of query aligned', 'hits in OGs', labels[1], colors[1])
    hist1(hits2['fqa'], 50, 'hitnum-fqa_filter2', 'fraction of query aligned', 'hits in OGs', labels[2], colors[2])

    # 8.2 OG histograms
    hist3([OG['fqa'].mean() for OG in OGs], 50, 'OGidnum-fqamean', 'mean FQA of hits in OG', 'OGs', labels, colors)
    hist3([OG['fqa'].var() for OG in OGs], 50, 'OGidnum-fqavar', 'variance of FQA of hits in OG', 'OGs', labels, colors)

    # 8.3 OG scatters
    scatter1(OGs[0]['fqa'].mean(), OGs[0]['fqa'].var(), 'fqavar-fqamean_all.png', 'FQA', labels[0], colors[0])
    scatter1(OGs[1]['fqa'].mean(), OGs[1]['fqa'].var(), 'fqavar-fqamean_filter1.png', 'FQA', labels[1], colors[1])
    scatter1(OGs[2]['fqa'].mean(), OGs[2]['fqa'].var(), 'fqavar-fqamean_filter2.png', 'FQA', labels[2], colors[2])

    # 9 cFQA PLOTS
    # 9.1 Hit histograms
    hist3([hit['cfqa'] for hit in hits], 50, 'hitnum-cfqa', 'compatible fraction of query aligned', 'hits in OGs', labels, colors)
    hist1(hits0['cfqa'], 50, 'hitnum-cfqa_all', 'compatible fraction of query aligned', 'hits in OGs', labels[0], colors[0])
    hist1(hits1['cfqa'], 50, 'hitnum-cfqa_filter1', 'compatible fraction of query aligned', 'hits in OGs', labels[1], colors[1])
    hist1(hits2['cfqa'], 50, 'hitnum-cfqa_filter2', 'compatible fraction of query aligned', 'hits in OGs', labels[2], colors[2])

    # 9.2 OG histograms
    hist3([OG['cfqa'].mean() for OG in OGs], 50, 'OGidnum-cfqamean', 'mean cFQA of hits in OG', 'OGs', labels, colors)
    hist3([OG['cfqa'].var() for OG in OGs], 50, 'OGidnum-cfqavar', 'variance of cFQA of hits in OG', 'OGs', labels, colors)

    # 9.3 OG scatters
    scatter1(OGs[0]['cfqa'].mean(), OGs[0]['cfqa'].var(), 'cfqavar-cfqamean_all.png', 'cFQA', labels[0], colors[0])
    scatter1(OGs[1]['cfqa'].mean(), OGs[1]['cfqa'].var(), 'cfqavar-cfqamean_filter1.png', 'cFQA', labels[1], colors[1])
    scatter1(OGs[2]['cfqa'].mean(), OGs[2]['cfqa'].var(), 'cfqavar-cfqamean_filter2.png', 'cFQA', labels[2], colors[2])

    # 10 xFQA PLOTS
    # 10.1 Hit histograms
    hist3([hit['xfqa'] for hit in hits], 50, 'hitnum-xfqa', 'excess fraction of query aligned', 'hits in OGs', labels, colors)
    hist1(hits0['xfqa'], 50, 'hitnum-xfqa_all', 'excess fraction of query aligned', 'hits in OGs', labels[0], colors[0])
    hist1(hits1['xfqa'], 50, 'hitnum-xfqa_filter1', 'excess fraction of query aligned', 'hits in OGs', labels[1], colors[1])
    hist1(hits2['xfqa'], 50, 'hitnum-xfqa_filter2', 'excess fraction of query aligned', 'hits in OGs', labels[2], colors[2])

    # 10.2 OG histograms
    hist3([OG['xfqa'].mean() for OG in OGs], 50, 'OGidnum-xfqamean', 'mean xFQA of hits in OG', 'OGs', labels, colors)
    hist3([OG['xfqa'].var() for OG in OGs], 50, 'OGidnum-xfqavar', 'variance of xFQA of hits in OG', 'OGs', labels, colors)

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
        plt.savefig(f'out/blast/hist2d_fsa-fqa_{label}.png')
        plt.close()

    # 12 EDGES
    edgenums = [hit[['qppid', 'sppid', 'OGid']].drop_duplicates().groupby('OGid').size() / 2 for hit in hits]
    ppidnums = [OG['qppid'].nunique() for OG in OGs]
    edgefracs = [2 * edgenum / (gnidnum*(gnidnum-1)) for edgenum, gnidnum in zip(edgenums, ppidnums)]

    # 12.1 Edge number histograms
    bar3(edgenums, 'OGidnum-edgenum', 'number of edges', 'OGs', labels, colors)
    bar1(edgenums[0], 'OGidnum-edgenum_all', 'number of edges', 'OGs', labels[0], colors[0])
    bar1(edgenums[1], 'OGidnum-edgenum_filter1', 'number of edges', 'OGs', labels[1], colors[1])
    bar1(edgenums[2], 'OGidnum-edgenum_filter2', 'number of edges', 'OGs', labels[2], colors[2])

    # 12.2 Edge fraction histograms
    hist3(edgefracs, 50, 'OGidnum-edgefrac', 'fraction of possible edges', 'OGs', labels, colors)
    hist1(edgefracs[0], 50, 'OGidnum-edgefrac_all', 'fraction of possible edges', 'OGs', labels[0], colors[0])
    hist1(edgefracs[1], 50, 'OGidnum-edgefrac_filter1', 'fraction of possible edges', 'OGs', labels[1], colors[1])
    hist1(edgefracs[2], 50, 'OGidnum-edgefrac_filter2', 'fraction of possible edges', 'OGs', labels[2], colors[2])

    # 13 CORRELATIONS
    scatter2(ppidnums[0], OGs[0]['bitscore'].mean(), 'bitscore-ppidnum', 'Mean bitscore of hits in OG')
    scatter2(ppidnums[0], OGs[0]['fqa'].mean(), 'fqa-ppidnum', 'Mean fraction of query aligned of hits in OG')
    scatter2(ppidnums[0], edgenums[0], 'edgenum-ppidnum', 'Number of edges in OG')
    scatter2(ppidnums[0], edgefracs[0], 'edgefrac-ppidnum', 'Fraction of possible edges in OG')

"""
DEPENDENCIES
../../ortho_search/hits2reciprocal/hits2reciprocal.py
    ../../ortho_search/hits2reciprocal/out/*/*.tsv
../../ortho_search/hsps2hits/hsps2hits.py
    ../../ortho_search/hsps2hits/out/*/*.tsv
../config/genomes.tsv
../cluster4+_graph/cluster4+_graph.py
    ../cluster4+_graph/out/4clique/clusters.tsv
"""