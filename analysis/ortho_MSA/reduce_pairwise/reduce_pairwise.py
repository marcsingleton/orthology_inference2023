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

    # Get lists of genes in OG for each species
    gnids_list = [tip.gnids for tip in tips]

    # Get sets of proteins corresponding to list of genes for each species
    ppid_sets_list = [[set([ppid for ppid, _, _, in gnid2seqs[gnid]]) for gnid in gnids] for gnids in gnids_list]

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

# Parse genomes
genomes = []
with open('../config/genomes.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, prot_path = line.split()
        genomes.append((spid, source, prot_path))

# Load seq metadata
gnid2spid = {}
ppid2gnid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, spid, _ = line.split()
        gnid2spid[gnid] = spid
        ppid2gnid[ppid] = gnid

# Load seqs
gnid2seqs = {}
for spid0, source, prot_path in genomes:
    with open(prot_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid0 = re.search(pp_regex[source], line).group(1)
                gnid = ppid2gnid[ppid0]
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq0 = ''.join(seqlines)

            try:
                for _, _, seq1 in gnid2seqs[gnid]:
                    if seq0 == seq1:
                        break
                else:
                    gnid2seqs[gnid].append((ppid0, spid0, seq0))
            except KeyError:
                gnid2seqs[gnid] = [(ppid0, spid0, seq0)]

# Load ggraph
ggraph = {}
with open('../../ortho_cluster3/hits2ggraph/out/ggraph2.tsv') as file:
    for line in file:
        node, adjs = line.rstrip('\n').split('\t')
        ggraph[node] = set([adj.split(':')[0] for adj in adjs.split(',')])

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

# Load OGs and tree
OGs = {}
with open('../../ortho_cluster3/clique4+_community/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        OGs[OGid] = set([node for edge in edges.split('\t') for node in edge.split(',')])

tree_template = skbio.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick', skbio.TreeNode)
tree_template.shear([tip.name for tip in tree_template.tips() if tip.name != 'sleb'])

rOGs = {}
for OGid, OG in OGs.items():
    # Add genes to leaves of tree
    tree = tree_template.deepcopy()
    tips = {tip.name: tip for tip in tree.tips()}
    remain = set()
    for gnid in OG:
        spid = gnid2spid[gnid]
        remain.add(spid)
        try:
            tips[spid].gnids.append(gnid)
        except AttributeError:
            tips[spid].gnids = [gnid]
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
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../config/genomes.tsv
../../ortho_cluster3/clique4+_community/clique4+_community2.py
    ../../ortho_cluster3/clique4+_community/out/ggraph2/5clique/gclusters.txt
../../ortho_cluster3/hits2ggraph/hits2ggraph2.py
    ../../ortho_cluster3/hits2ggraph/out/ggraph2.tsv
../../ortho_cluster3/hits2pgraph/hits2pgraph.py
    ../../ortho_cluster3/hits2pgraph/out/pgraph2.tsv
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../../ortho_tree/consensus_tree/consensus_tree.py
    ../../ortho_tree/consensus_tree/out/100red_ni.txt
"""