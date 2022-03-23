"""Make FASTAs of remaining OGs using longest isoforms."""

import os
import re

import pandas as pd
from src.utils import read_fasta

ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+)'}

# Parse genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, prot_path = line.split()
        genomes.append((spid, source, prot_path))

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, sqid = line.split()
        ppid2meta[ppid] = (gnid, spid, sqid)

# Load seqs
ppid2seq = {}
for spid, source, prot_path in genomes:
    fasta = read_fasta(prot_path)
    for header, seq in fasta:
        ppid = re.search(ppid_regex[header], line).group(1)
        ppid2seq[ppid] = seq

# Load OGs and OG metadata
OGs = {}
with open('../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        _, OGid, _, edges = line.rstrip().split('\t')
        sqids = {ppid2meta[node][2] for edge in edges.split(',') for node in edge.split(':')}
        OGs[OGid] = sqids  # Ensure only representatives are selected for reduced clusters
OG_meta = pd.read_table('../OG_meta/out/OG_meta.tsv')

# Write sequences
if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = OG_meta.loc[~(OG_meta['sqidnum'] == OG_meta['gnidnum']), 'OGid']
for OGid in OGids:
    gnid2sqids = {}
    for sqid in OGs[OGid]:
        gnid, spid, _ = ppid2meta[sqid]
        seq = ppid2seq[sqid]
        try:
            gnid2sqids[gnid].append((sqid, seq))
        except KeyError:
            gnid2sqids[gnid] = [(sqid, seq)]
    records = []
    for gnid, sqids in gnid2sqids.items():
        sqid, seq = max(sqids, key=lambda x: len(x[1]))
        _, spid, _ = ppid2meta[sqid]
        seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]) + '\n'
        records.append((sqid, gnid, spid, seqstring))
    with open(f'out/{OGid}.tfa', 'w') as file:
        for sqid, gnid, spid, seqstring in sorted(records, key=lambda x: x[2]):
            file.write(f'>ppid={sqid}|gnid={gnid}|spid={spid}\n' + seqstring)

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../../ortho_cluster3/cluster4+_graph/cluster4+_graph.py
    ../../ortho_cluster3/cluster4+_graph/out/4clique/clusters.tsv
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../config/genomes.tsv
../OG_meta/OG_meta.py
    ../OG_meta/out/OG_meta.tsv
"""