"""Reduce gene clusters to representative polypeptides using pairwise comparisons."""

import os
import re
from itertools import chain, product, combinations

import skbio


def reduce_node(node):
    """Return representative ppids for all gnids associated with terminals in node."""
    is_tip = [child.is_tip() for child in node.children]

    if any(is_tip):
        distances = get_distances(node)
        tips = [child for child in node.children if child.is_tip()]
        nontips = [ppid for child in node.children for ppid in reduce_node(child) if not child.is_tip()]
        rOG = reduce_terminal(tips, nontips, distances)
    else:
        rOG = [ppid for child in node.children for ppid in reduce_node(child)]

    return rOG


def reduce_terminal(tips, nontips, distances):
    """Return list of tuples containing genes at least one terminal node."""
    # Get lists of species ids
    spids_list = [[tip.name for _ in tip.gnids] for tip in tips]

    # Get lists of genes in pOG for each species
    gnids_list = [tip.gnids.keys() for tip in tips]

    # Get sets of proteins corresponding to list of genes for each species
    ppid_sets_list = [tip.gnids.values() for tip in tips]

    # Get all combinations of proteins with each gene represented once for each species
    ppid_prods_list = [product(*ppid_sets) for ppid_sets in ppid_sets_list]

    scores = []
    for ppid_prods in product(*ppid_prods_list):  # Separate species to sum only inter-species hits
        # Add species and gene ids to locate in pOG and calculate products to sum hits bidirectionally
        zs = [list(zip(spids, gnids, ppid_prod)) for spids, gnids, ppid_prod in zip(spids_list, gnids_list, ppid_prods)]

        # Get species combinations then ppid combinations with terminals
        spid_pairs = combinations(zs, 2)
        ppid_pairs1 = [ppid_pair for spid_pair in spid_pairs for ppid_pair in product(*spid_pair)]

        # Get ppid combinations with nonterminals
        ppid_pairs2 = product([ppid for z in zs for ppid in z], nontips)

        score = 0
        for (spid1, gnid1, ppid1), (spid2, gnid2, ppid2) in chain(ppid_pairs1, ppid_pairs2):
            weight = distances[spid1] + distances[spid2]
            try:
                score += pgraph[ppid1][ppid2] / weight
            except KeyError:
                pass
            try:
                score += pgraph[ppid2][ppid1] / weight
            except KeyError:
                pass
        scores.append((score, [ppid for z in zs for ppid in z] + nontips))

    return max(scores)[1]


def get_distances(node, d0=0):
    """Return distances to terminals from current node."""
    distances = {}
    for child in node.children:
        if child.is_tip():
            distances[child.name] = d0 + child.length
        else:
            distances.update(get_distances(child, d0 + child.length))

    return distances


pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+)'}

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, spid, _ = line.split()
        ppid2meta[ppid] = (gnid, spid)

# Load pgraph
pgraph = {}
with open('../../ortho_cluster3/hits2pgraph/out/pgraph2.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        bitscores = {}
        for adj in adjs.split(','):
            adj_node, adj_bitscore = adj.split(':')
            bitscores[adj_node] = float(adj_bitscore)
        pgraph[node] = bitscores

# Load pOGs and tree
pOGs = {}
with open('../../ortho_cluster3/subcluster_pgraph/out/pclusters.txt') as file:
    for line in file:
        _, pOGid, edges = line.rstrip().split(':')
        ppids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        pOGs[pOGid] = ppids

tree_template = skbio.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick', skbio.TreeNode)
tree_template.shear([tip.name for tip in tree_template.tips() if tip.name != 'sleb'])

rOGs = {}
for pOGid, pOG in pOGs.items():
    # Add genes to leaves of tree
    tree = tree_template.deepcopy()
    tips = {tip.name: tip for tip in tree.tips()}
    remain = set()
    for ppid in pOG:
        gnid, spid = ppid2meta[ppid]
        remain.add(spid)
        try:
            tips[spid].gnids[gnid].append(ppid)
        except KeyError:
            tips[spid].gnids[gnid] = [ppid]
        except AttributeError:
            tips[spid].gnids = {gnid: [ppid]}
    tree = tree.shear(remain)

    # Reduce pOG
    rOG = reduce_node(tree)
    rOGs[pOGid] = rOG

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write reduced clusters to file
with open('out/rclusters.tsv', 'w') as outfile:
    outfile.write('pOGid\tspid\tgnid\tppid\n')
    for OGid, rOG in rOGs.items():
        for entry in rOG:
            outfile.write(OGid + '\t' + '\t'.join(entry) + '\n')

"""
DEPENDENCIES
../../ortho_cluster3/hits2pgraph/hits2pgraph.py
    ../../ortho_cluster3/hits2pgraph/out/pgraph2.tsv
../../ortho_cluster3/subcluster_pgraph/subcluster_pgraph.py
    ../../ortho_cluster3/subcluster_pgraph/out/pclusters.txt
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_tree/consensus_tree/consensus_tree.py
    ../../ortho_tree/consensus_tree/out/100red_ni.txt
"""