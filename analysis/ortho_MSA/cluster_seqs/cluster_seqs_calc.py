"""Cluster unique sequences in genes."""

import os
import re

from Bio import pairwise2
from Bio.Align import substitution_matrices


def get_pident(alignment):
    aligned1, aligned2 = alignment[:2]
    num_sym = 0
    num_ident = 0
    for sym1, sym2 in zip(aligned1, aligned2):
        num_sym += 1
        if sym1 == sym2:
            num_ident += 1
    return 100 * num_ident / num_sym


pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+)'}
blosum62 = substitution_matrices.load('blosum62')
for key, value in list(blosum62.items()):  # Necessary to prevent RunTimeError due to dictionary size change
    if 'C' in key:
        idx = key.index('C')
        sym = key[(idx+1) % 2]
        blosum62[('U', sym)] = value
        blosum62[(sym, 'U')] = value

# Parse parameters
params = []
with open('params.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, prot_path = line.split()
        params.append((spid, source, prot_path))

# Load pp metadata
ppid2gnid = {}
with open('../../ortho_cluster2/ppid2meta/out/ppid2meta.tsv') as file:
    for line in file:
        ppid, gnid, _ = line.split()
        ppid2gnid[ppid] = gnid

# Load seqs
gnid2seqs = {}
repr2cons = {}  # Mapping of representative PPIDs to constituent sequences
for _, source, prot_path in params:
    with open(prot_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid0 = re.search(pp_regex[source], line).group(1)
                gnid0 = ppid2gnid[ppid0]
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq0 = ''.join(seqlines)

            try:
                for ppid1, seq1 in gnid2seqs[gnid0]:
                    if seq0 == seq1:
                        repr2cons[ppid1].append(ppid0)
                        break
                else:
                    gnid2seqs[gnid0].append((ppid0, seq0))
                    repr2cons[ppid0] = []
            except KeyError:
                gnid2seqs[gnid0] = [(ppid0, seq0)]
                repr2cons[ppid0] = []

# Cluster sequences
clusters = []
for gnid, seqs in list(gnid2seqs.items()):
    seqs = sorted(seqs, key=lambda x: (len(x[1]), x))  # Sort by length then sequence
    gene_clusters = []
    for ppid0, seq0 in seqs:
        pidents = []
        for idx, (ppid1, _, seq1, _) in enumerate(gene_clusters):
            alignments = pairwise2.align.globalds(seq0, seq1, blosum62, -11, -1)
            pident = max([get_pident(alignment) for alignment in alignments])  # Get traceback with best pident
            pidents.append((idx, pident))
        idx, pident = max(pidents, key=lambda x: x[1]) if pidents else (0, 0)
        if pident > 90:
            cluster_ppids = gene_clusters[idx][3]
            cluster_ppids.append((ppid0, pident))
        else:
            gene_clusters.append((ppid0, gnid, seq0, []))
    clusters.extend(gene_clusters)

# Save clusters
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/clusters.tsv', 'w') as file:
    file.write('ppid\tgnid\trepr_id\tcluster_id\tpident\n')  # Identical and cluster PPIDs
    for c_ppid, gnid, _, ppids in clusters:
        file.write('\t'.join([c_ppid, gnid, c_ppid, c_ppid, '100']) + '\n')  # Entry for representative of cluster
        for i_ppid, pident in ppids:
            file.write('\t'.join([i_ppid, gnid, i_ppid, c_ppid, str(pident)]) + '\n')  # Entry for representative of identicals
            for ppid in repr2cons[i_ppid]:
                file.write('\t'.join([ppid, gnid, i_ppid, c_ppid, str(pident)]) + '\n')  # Entry for identical
