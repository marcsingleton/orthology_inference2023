"""Make FASTAs of OGs with one unique sequence per gene."""

import os
import re

import pandas as pd
from src.utils import read_fasta

ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+)'}

# Load genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        _, _, source, prot_path, _ = line.rstrip('\n').split('\t')
        genomes.append((source, prot_path))

# Load sequence data
ppid2data = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, _ = line.rstrip('\n').split('\t')
        ppid2data[ppid] = (gnid, spid)

# Load seqs
ppid2seq = {}
for source, prot_path in genomes:
    fasta = read_fasta(prot_path)
    for header, seq in fasta:
        ppid = re.search(ppid_regex[header], line).group(1)
        ppid2seq[ppid] = seq

# Load OGs and OG data
OGs = {}
with open('../cluster4+_graph/out/4clique/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        _, OGid, _, edges = line.rstrip('\n').split('\t')
        ppids = {node for edge in edges.split(',') for node in edge.split(':')}
        OGs[OGid] = ppids
OG_data = pd.read_table('../OG_data/out/OG_data.tsv')

# Write sequences
if not os.path.exists('out/'):
    os.mkdir('out/')

spidnum = len(genomes)
ppid_filter = OG_data['ppidnum'] == spidnum
gnid_filter = OG_data['gnidnum'] == spidnum
spid_filter = OG_data['spidnum'] == spidnum
OGids = OG_data.loc[ppid_filter & gnid_filter & spid_filter, 'OGid']

for OGid in OGids:
    records = []
    for ppid in OGs[OGid]:
        gnid, spid = ppid2data[ppid]
        seq = ppid2seq[ppid]
        seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)])
        records.append((ppid, gnid, spid, seqstring))
    with open(f'out/{OGid}.fa', 'w') as file:
        for ppid, gnid, spid, seqstring in sorted(records, key=lambda x: x[2]):
            file.write(f'>ppid={ppid}|gnid={gnid}|spid={spid}\n{seqstring}\n')

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/fasta/dmel-all-translation-r6.38.fasta
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../config/genomes.tsv
../cluster4+_graph/cluster.py
    ../cluster4+_graph/out/4clique/clusters.tsv
../OG_data/OG_data.py
    ../OG_data/out/OG_data.tsv
"""