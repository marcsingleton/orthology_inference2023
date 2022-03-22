"""Reduce gene clusters to representative polypeptides using pairwise comparisons."""

import os
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

    # Get lists of genes in OG for each species
    gnids_list = [tip.gnids.keys() for tip in tips]

    # Get sets of proteins corresponding to list of genes for each species
    ppid_sets_list = [tip.gnids.values() for tip in tips]

    # Get all combinations of proteins with each gene represented once for each species
    ppid_prods_list = [product(*ppid_sets) for ppid_sets in ppid_sets_list]

    scores = []
    for ppid_prods in product(*ppid_prods_list):  # Separate species to sum only inter-species hits
        # Add species and gene ids to locate in OG and calculate products to sum hits bidirectionally
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
                score += graph[ppid1][ppid2] / weight
            except KeyError:
                pass
            try:
                score += graph[ppid2][ppid1] / weight
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


# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, sqid = line.split()
        ppid2meta[ppid] = (gnid, spid, sqid)

# Load graph
graph = {}
with open('../../ortho_cluster3/hits2graph/out/hit_graph.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        bitscores = {}
        for adj in adjs.split(','):
            adj_node, adj_bitscore = adj.split(':')
            bitscores[adj_node] = float(adj_bitscore)
        graph[node] = bitscores

# Load OGs and tree
OGs = {}
with open('../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        sqids = {ppid2meta[node][2] for edge in edges.split('\t') for node in edge.split(',')}
        OGs[OGid] = sqids  # Ensure only representatives are selected for reduced clusters

tree_template = skbio.read('../../ortho_tree/ctree_WAG/out/100red_ni.txt', 'newick', skbio.TreeNode)
tree_template.shear([tip.name for tip in tree_template.tips() if tip.name != 'sleb'])

rOGs = {}
for OGid, OG in OGs.items():
    # Add genes to leaves of tree
    tree = tree_template.deepcopy()
    tips = {tip.name: tip for tip in tree.tips()}
    remain = set()
    for sqid in OG:
        gnid, spid, _ = ppid2meta[sqid]
        remain.add(spid)
        try:
            tips[spid].gnids[gnid].append(sqid)
        except KeyError:
            tips[spid].gnids[gnid] = [sqid]
        except AttributeError:
            tips[spid].gnids = {gnid: [sqid]}
    tree = tree.shear(remain)

    # Reduce OG
    rOG = reduce_node(tree)
    rOGs[OGid] = rOG

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write reduced clusters to file
with open('out/rclusters.tsv', 'w') as outfile:
    outfile.write('OGid\tspid\tgnid\tppid\n')
    for OGid, rOG in rOGs.items():
        for entry in rOG:
            outfile.write(OGid + '\t' + '\t'.join(entry) + '\n')

"""
DEPENDENCIES
../../ortho_cluster3/clique4+_pcommunity/clique4+_pcommunity2.py
    ../../ortho_cluster3/clique4+_pcommunity/out/pgraph2/4clique/pclusters.txt
../../ortho_cluster3/hits2graph/hits2graph.py
    ../../ortho_cluster3/hits2graph/out/hit_graph.tsv
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../../ortho_tree/ctree_WAG/ctree_WAG.py
    ../../ortho_tree/ctree_WAG/out/100red_ni.txt
"""