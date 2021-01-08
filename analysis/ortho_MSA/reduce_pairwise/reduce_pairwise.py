"""Reduce gene clusters to representative polypeptides using pairwise comparisons."""

import os
import re
from copy import deepcopy
from itertools import chain, product, combinations

import Bio.Phylo as Phylo


def reduce_node(node):
    """Return representative ppids for all gnids associated with terminals in node."""
    is_terminal = [leaf.is_terminal() for leaf in node.clades]

    if any(is_terminal):
        distances = get_distances(node)
        terminals = [leaf for leaf in node.clades if leaf.is_terminal()]
        nonterminals = [ppid for leaf in node.clades for ppid in reduce_node(leaf) if not leaf.is_terminal()]
        rOG = reduce_terminal(terminals, nonterminals, distances)
    else:
        rOG = [ppid for leaf in node.clades for ppid in reduce_node(leaf)]

    return rOG


def reduce_terminal(terminals, nonterminals, distances):
    """Return list of tuples containing genes at least one terminal node."""
    # Get lists of species ids
    spids_list = [[terminal.name for _ in terminal.gnids] for terminal in terminals]

    # Get lists of genes in OG for each species
    gnids_list = [terminal.gnids for terminal in terminals]

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
        ppid_pairs2 = product([ppid for z in zs for ppid in z], nonterminals)

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
        scores.append((score, [ppid for z in zs for ppid in z] + nonterminals))

    return max(scores)[1]


def get_distances(node, d0=0):
    """Return distances to terminals from current node."""
    distances = {}
    for leaf in node.clades:
        if leaf.is_terminal():
            distances[leaf.name] = d0 + leaf.branch_length
        else:
            distances.update(get_distances(leaf, d0 + leaf.branch_length))

    return distances


pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+)'}

# Parse parameters
params = []
with open('params.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, prot_path = line.split()
        params.append((spid, source, prot_path))

# Load gn and pp metadata
gnid2spid = {}
ppid2gnid = {}
with open('../../ortho_search/ppid2meta/out/ppid2meta.tsv') as file:
    for line in file:
        ppid, gnid, spid = line.split()
        gnid2spid[gnid] = spid
        ppid2gnid[ppid] = gnid

# Load seqs
gnid2seqs = {}
for spid0, source, prot_path in params:
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
with open('../../ortho_cluster3/clique5+_community/out/ggraph2/5clique/gclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        OGs[OGid] = set([node for edge in edges.split('\t') for node in edge.split(',')])

tree_template = Phylo.read('../../ortho_tree/consensus_tree/out/100red_ni.txt', 'newick')
tree_template.prune('sleb')

rOGs = {}
for OGid, OG in OGs.items():
    # Add genes to leaves of tree
    tree = deepcopy(tree_template)
    terminals = {terminal.name: terminal for terminal in tree.get_terminals()}
    remain = set()
    for gnid in OG:
        spid = gnid2spid[gnid]
        remain.add(spid)
        try:
            terminals[spid].gnids.append(gnid)
        except AttributeError:
            terminals[spid].gnids = [gnid]

    # Remove unused terminals
    remove = set(terminals) - remain
    for terminal in remove:
        tree.prune(terminal)

    # Reduce OG
    rOG = reduce_node(tree.root)  # Pass root since tree object has no clades attribute
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
../../ortho_cluster3/clique5+_community/clique5+_community2.py
    ../../ortho_cluster3/clique5+_community/out/ggraph2/5clique/gclusters.txt
../../ortho_cluster3/hits2ggraph/hits2ggraph2.py
    ../../ortho_cluster3/hits2ggraph/out/ggraph2.tsv
../../ortho_cluster3/hits2pgraph/hits2pgraph.py
    ../../ortho_cluster3/hits2pgraph/out/pgraph2.tsv
../../ortho_search/ppid2meta/ppid2meta.py
    ../../ortho_search/ppid2meta/out/ppid2meta.tsv
../../ortho_tree/consensus_tree/consensus_tree.py
    ../../ortho_tree/consensus_tree/out/100red_ni.txt
./params.tsv
"""