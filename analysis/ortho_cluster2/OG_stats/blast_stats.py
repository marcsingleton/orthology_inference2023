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
    df = pd.read_table(f'../../ortho_search/hsps2hits/out/{qspid}/{sspid}.tsv',
                       usecols=hit_dtypes.keys(), dtype=hit_dtypes, memory_map=True)
    r = pd.read_table(f'../../ortho_search/hits2reciprocal/out/{qspid}/{sspid}.tsv',
                      usecols=['reciprocal'], memory_map=True)
    df = df[r['reciprocal']]

    df['qspid'] = qspid
    df['qspid'] = df['qspid'].astype('category')
    df['sspid'] = sspid
    df['sspid'] = df['sspid'].astype('category')

    return df


# Plot functions
def hist1(df, bins, file_label, x_label, y_label, df_label, color):
    plt.hist(df, bins=bins, label=df_label, color=color)
    plt.xlabel(x_label)
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
    axs[2].set_xlabel(x_label)
    axs[1].set_ylabel(f'Number of {y_label}')
    fig.subplots_adjust(left=0.175)
    fig.savefig(f'out/blast/hist_{file_label}_3.png')
    plt.close()


def bar1(df, file_label, x_label, y_label, df_label, color):
    plt.bar(df.index, df.values, label=df_label, color=color, width=1)
    plt.xlabel(x_label)
    plt.ylabel(f'Number of {y_label}')
    plt.legend()
    plt.savefig(f'out/blast/hist_{file_label}.png')
    plt.close()


def bar3(dfs, file_label, x_label, y_label, df_labels, colors):
    fig, axs = plt.subplots(3, 1, figsize=(4.8, 6), sharex=True)
    for ax, df, df_label, color in zip(axs, dfs, df_labels, colors):
        ax.bar(df.index, df.values, label=df_label, color=color, width=1)
        ax.legend()
    axs[2].set_xlabel(x_label)
    axs[1].set_ylabel(f'Number of {y_label}')
    fig.subplots_adjust(left=0.175)
    fig.savefig(f'out/blast/hist_{file_label}_3.png')
    plt.close()


def hexbin1(x, y, file_label, xy_label, df_label):
    fig, ax = plt.subplots()
    hb = ax.hexbin(x, y, bins='log', gridsize=100, mincnt=1, linewidth=0)
    ax.set_xlabel(f'Mean {xy_label} in OG')
    ax.set_ylabel(f'Variance of {xy_label} in OG')
    ax.set_title(df_label)
    fig.colorbar(hb)
    fig.savefig(f'out/blast/hexbin_{file_label}.png')
    plt.close()


def hexbin2(x, y, file_label, y_label, df_label):
    fig, ax = plt.subplots()
    hb = ax.hexbin(x, y, bins='log', gridsize=100, mincnt=1, linewidth=0)
    ax.set_xlabel('Number of proteins in OG')
    ax.set_ylabel(y_label)
    ax.set_title(df_label)
    fig.colorbar(hb)
    plt.savefig(f'out/blast/hexbin_{file_label}.png')
    plt.close()


hit_dtypes = {'qppid': 'string', 'qgnid': 'string',
              'sppid': 'string', 'sgnid': 'string',
              'hspnum': int, 'chspnum': int,
              'qlen': int, 'nqa': int, 'cnqa': int,
              'slen': int, 'nsa': int, 'cnsa': int,
              'bitscore': float}
num_processes = 4

if __name__ == '__main__':
    # Load genomes
    spids = []
    with open('../config/genomes.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            spids.append(fields['spid'])

    # Load data
    rows = []
    with open('../cluster4+_graph/out/4clique/clusters.tsv') as file:
        field_names = file.readline().rstrip('\n').split('\t')
        for line in file:
            fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
            component_id, OGid = fields['component_id'], fields['OGid']
            for edge in fields['edges'].split(','):
                node1, node2 = edge.split(':')
                rows.append({'component_id': component_id, 'OGid': OGid, 'qppid': node1, 'sppid': node2})
                rows.append({'component_id': component_id, 'OGid': OGid, 'qppid': node2, 'sppid': node1})
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
    OGid_groups = hits0.groupby('OGid')
    ppid_groups = hits0[['qppid', 'OGid']].drop_duplicates().groupby('qppid')
    spidnum = OGid_groups['qspid'].nunique()
    gnidnum = OGid_groups['qgnid'].nunique()

    spid_filter = spidnum[spidnum == len(spids)].index  # All species represented
    gnid_filter = gnidnum[gnidnum == len(spids)].index  # One gene per species
    ppid_filter = ppid_groups.filter(lambda x: len(x) == 1)['OGid']  # One OG per protein

    OGs1 = set(spid_filter) & set(gnid_filter)
    OGs2 = OGs1 & set(ppid_filter)

    hits1 = hits0[hits0['OGid'].isin(OGs1)]
    hits2 = hits0[hits0['OGid'].isin(OGs2)]
    hits = [hits0, hits1, hits2]
    OGs = [hit.groupby('OGid') for hit in hits]

    file_labels = ['all', 'filter1', 'filter2']
    labels = ['all', 'filter1', 'filter2']
    colors = ['C0', 'C1', 'C2']

    if not os.path.exists('out/blast/'):
        os.makedirs('out/blast/')

    # 1 FILTER PLOT
    ys = [hit['OGid'].nunique() for hit in hits]
    plt.bar(labels, ys, color=colors, width=0.25)
    plt.xlim((-0.75, 2.75))
    plt.xlabel('OG subset')
    plt.ylabel('Number of OGs')
    plt.savefig('out/blast/bar_OGidnum-filter.png')
    plt.close()

    # 2 BITSCORE PLOTS
    # 2.1 Hit histograms
    hist3([hit['bitscore'] for hit in hits], 200, 'hitnum-bitscore', 'Bitscore', 'hits in OGs', labels, colors)
    for hit, file_label, label, color in zip(hits, file_labels, labels, colors):
        hist1(hit['bitscore'], 200, f'hitnum-bitscore_{file_label}', 'Bitscore', 'hits in OGs', label, color)

    # 2.2 OG histograms
    hist3([OG['bitscore'].mean() for OG in OGs], 200, 'OGidnum-bitscoremean', 'Mean bitscore of hits in OG', 'OGs', labels, colors)
    hist3([OG['bitscore'].var() for OG in OGs], 200, 'OGidnum-bitscorevar', 'Variance of bitscore of hits in OG', 'OGs', labels, colors)

    # 2.3 OG scatters
    for OG, file_label, label, color in zip(OGs, file_labels, labels, colors):
        hexbin1(OG['bitscore'].mean(), OG['bitscore'].var(), f'bitscorevar-bitscoremean_{file_label}', 'bitscore', label)

    # 3 HSPNUM HISTOGRAMS
    counts = [hit['hspnum'].value_counts() for hit in hits]
    bar3(counts, 'hitnum-hspnum', 'Number of disjoint HSPs in hit', 'hits', labels, colors)
    for count, file_label, label, color in zip(counts, file_labels, labels, colors):
        bar1(count, f'hitnum-hspnum_{file_label}', 'Number of disjoint HSPs in hit', 'hits', label, color)

    # 4 cHSPNUM HISTOGRAMS
    counts = [hit['chspnum'].value_counts() for hit in hits]
    bar3(counts, 'hitnum-chspnum', 'Number of compatible HSPs in hit', 'hits', labels, colors)
    for count, file_label, label, color in zip(counts, file_labels, labels, colors):
        bar1(count, f'hitnum-chspnum_{file_label}', 'Number of compatible HSPs in hit', 'hits', label, color)

    # 5 xHSPNUM HISTOGRAMS
    counts = [hit['xhspnum'].value_counts() for hit in hits]
    bar3(counts, 'hitnum-xhspnum', 'Excess number of HSPs in hit', 'hits', labels, colors)
    for count, file_label, label, color in zip(counts, file_labels, labels, colors):
        bar1(count, f'hitnum-xhspnum_{file_label}', 'Excess number of HSPs in hit', 'hits', label, color)

    # 6 NQA PLOTS
    # 6.1 Hit histograms
    hist3([hit['nqa'] for hit in hits], 200, 'hitnum-nqa', 'Number of query aligned', 'hits in OGs', labels, colors)
    for hit, file_label, label, color in zip(hits, file_labels, labels, colors):
        hist1(hit['nqa'], 200, f'hitnum-nqa_{file_label}', 'Number of query aligned', 'hits in OGs', label, color)

    # 6.2 OG histograms
    hist3([OG['nqa'].mean() for OG in OGs], 50, 'OGidnum-nqamean', 'Mean NQA of hits in OG', 'OGs', labels, colors)
    hist3([OG['nqa'].var() for OG in OGs], 50, 'OGidnum-nqavar', 'Variance of NQA of hits in OG', 'OGs', labels, colors)

    # 6.3 OG scatters
    for OG, file_label, label, color in zip(OGs, file_labels, labels, colors):
        hexbin1(OG['nqa'].mean(), OG['nqa'].var(), f'nqavar-nqamean_{file_label}', 'NQA', label)

    # 7 cNQA PLOTS
    # 7.1 Hit histograms
    hist3([hit['cnqa'] for hit in hits], 200, 'hitnum-cnqa', 'Compatible number of query aligned', 'hits in OGs', labels, colors)
    for hit, file_label, label, color in zip(hits, file_labels, labels, colors):
        hist1(hit['cnqa'], 200, f'hitnum-cnqa_{file_label}', 'Compatible number of query aligned', 'hits in OGs', label, color)

    # 7.2 OG histograms
    hist3([OG['cnqa'].mean() for OG in OGs], 50, 'OGidnum-cnqamean', 'Mean cNQA of hits in OG', 'OGs', labels, colors)
    hist3([OG['cnqa'].var() for OG in OGs], 50, 'OGidnum-cnqavar', 'Variance of cNQA of hits in OG', 'OGs', labels, colors)

    # 7.3 OG scatters
    for OG, file_label, label, color in zip(OGs, file_labels, labels, colors):
        hexbin1(OG['cnqa'].mean(), OG['cnqa'].var(), f'cnqavar-cnqamean_{file_label}', 'cNQA', label)

    # 8 FQA PLOTS
    # 8.1 Hit histograms
    hist3([hit['fqa'] for hit in hits], 50, 'hitnum-fqa', 'Fraction of query aligned', 'hits in OGs', labels, colors)
    for hit, file_label, label, color in zip(hits, file_labels, labels, colors):
        hist1(hit['fqa'], 50, f'hitnum-fqa_{file_label}', 'Fraction of query aligned', 'hits in OGs', label, color)

    # 8.2 OG histograms
    hist3([OG['fqa'].mean() for OG in OGs], 50, 'OGidnum-fqamean', 'Mean FQA of hits in OG', 'OGs', labels, colors)
    hist3([OG['fqa'].var() for OG in OGs], 50, 'OGidnum-fqavar', 'Variance of FQA of hits in OG', 'OGs', labels, colors)

    # 8.3 OG scatters
    for OG, file_label, label, color in zip(OGs, file_labels, labels, colors):
        hexbin1(OG['fqa'].mean(), OG['fqa'].var(), f'fqavar-fqamean_{file_label}', 'FQA', label)

    # 9 cFQA PLOTS
    # 9.1 Hit histograms
    hist3([hit['cfqa'] for hit in hits], 50, 'hitnum-cfqa', 'Compatible fraction of query aligned', 'hits in OGs', labels, colors)
    for hit, file_label, label, color in zip(hits, file_labels, labels, colors):
        hist1(hit['cfqa'], 50, f'hitnum-cfqa_{file_label}', 'Compatible fraction of query aligned', 'hits in OGs', label, color)

    # 9.2 OG histograms
    hist3([OG['cfqa'].mean() for OG in OGs], 50, 'OGidnum-cfqamean', 'Mean cFQA of hits in OG', 'OGs', labels, colors)
    hist3([OG['cfqa'].var() for OG in OGs], 50, 'OGidnum-cfqavar', 'Variance of cFQA of hits in OG', 'OGs', labels, colors)

    # 9.3 OG scatters
    for OG, file_label, label, color in zip(OGs, file_labels, labels, colors):
        hexbin1(OG['cfqa'].mean(), OG['cfqa'].var(), f'cfqavar-cfqamean_{file_label}', 'cFQA', label)

    # 10 xFQA PLOTS
    # 10.1 Hit histograms
    hist3([hit['xfqa'] for hit in hits], 50, 'hitnum-xfqa', 'Excess fraction of query aligned', 'hits in OGs', labels, colors)
    for hit, file_label, label, color in zip(hits, file_labels, labels, colors):
        hist1(hit['xfqa'], 50, f'hitnum-xfqa_{file_label}', 'Excess fraction of query aligned', 'hits in OGs', label, color)

    # 10.2 OG histograms
    hist3([OG['xfqa'].mean() for OG in OGs], 50, 'OGidnum-xfqamean', 'Mean xFQA of hits in OG', 'OGs', labels, colors)
    hist3([OG['xfqa'].var() for OG in OGs], 50, 'OGidnum-xfqavar', 'Variance of xFQA of hits in OG', 'OGs', labels, colors)

    # 10.3 OG scatters
    for OG, file_label, label, color in zip(OGs, file_labels, labels, colors):
        hexbin1(OG['xfqa'].mean(), OG['xfqa'].var(), f'xfqavar-xfqamean_{file_label}', 'xFQA', label)

    # 11 FQA-FSA SCATTERS
    for hit, file_label, label in zip(hits, file_labels, labels):
        plt.hist2d(hit['fqa'], hit['fsa'], bins=50, norm=mpl_colors.LogNorm())
        plt.xlabel('Fraction of query aligned')
        plt.ylabel('Fraction of subject aligned')
        plt.title(label)
        plt.colorbar()
        plt.savefig(f'out/blast/hist2d_fsa-fqa_{file_label}.png')
        plt.close()

    # 12 EDGES
    edgenums = [hit[['qppid', 'sppid', 'OGid']].drop_duplicates().groupby('OGid').size() / 2 for hit in hits]
    ppidnums = [OG['qppid'].nunique() for OG in OGs]
    edgefracs = [2 * edgenum / (gnidnum*(gnidnum-1)) for edgenum, gnidnum in zip(edgenums, ppidnums)]

    # 12.1 Edge number histograms
    hist3(edgenums, 50, 'OGidnum-edgenum', 'Number of edges', 'OGs', labels, colors)
    for edgenum, file_label, label, color in zip(edgenums, file_labels, labels, colors):
        hist1(edgenum, 50, f'OGidnum-edgenum_{file_label}', 'Number of edges', 'OGs', label, color)

    # 12.2 Edge fraction histograms
    hist3(edgefracs, 50, 'OGidnum-edgefrac', 'Fraction of possible edges', 'OGs', labels, colors)
    for edgefrac, file_label, label, color in zip(edgefracs, file_labels, labels, colors):
        hist1(edgefrac, 50, f'OGidnum-edgefrac_{file_label}', 'Fraction of possible edges', 'OGs', label, color)

    # 13 CORRELATIONS
    hexbin2(ppidnums[1], OGs[1]['bitscore'].mean(), f'bitscore-ppidnum_{file_labels[1]}', 'Mean bitscore of hits in OG', labels[1])
    hexbin2(ppidnums[1], OGs[1]['fqa'].mean(), f'fqa-ppidnum_{file_labels[1]}', 'Mean FQA of hits in OG', labels[1])
    hexbin2(ppidnums[1], edgenums[1], f'edgenum-ppidnum_{file_labels[1]}', 'Number of edges in OG', labels[1])
    hexbin2(ppidnums[1], edgefracs[1], f'edgefrac-ppidnum_{file_labels[1]}', 'Fraction of possible edges in OG', labels[1])

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