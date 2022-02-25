"""Make FASTAs of OGs with one unique sequence per gene."""

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
        _, _, source, prot_path, _ = line.split()
        genomes.append((source, prot_path))

# Load seq metadata
ppid2meta = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, sqid = line.split()
        ppid2meta[ppid] = (gnid, spid, sqid)

# Load seqs
ppid2seq = {}
for source, prot_path in genomes:
    fasta = read_fasta(prot_path)
    for header, seq in fasta:
        ppid = re.search(ppid_regex[header], line).group(1)
        ppid2seq[ppid] = seq

# Load OGs and OG metadata
OGs = {}
with open('../clique4+_pcommunity/out/4clique/pclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        sqids = {ppid2meta[node][2] for edge in edges.split('\t') for node in edge.split(',')}
        OGs[OGid] = sqids
OGs_meta = pd.read_table('../OG_meta/out/OG_meta.tsv')

# Write sequences
num = len(genomes)
gnnum = OGs_meta['gnidnum'] == num
spnum = OGs_meta['spidnum'] == num
sqnum = OGs_meta['sqidnum'] == num
OGids = OGs_meta.loc[gnnum & spnum & sqnum, 'OGid']

if not os.path.exists('out/'):
    os.mkdir('out/')

for OGid in OGids:
    records = []
    for sqid in OGs[OGid]:
        gnid, spid, _ = ppid2meta[sqid]
        seq = ppid2seq[sqid]
        records.append((seq, sqid, gnid, spid))
    with open(f'out/{OGid}.tfa', 'w') as file:
        for seq, sqid, gnid, spid in sorted(records, key=lambda x: x[3]):
            seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]) + '\n'
            file.write(f'>ppid={sqid}|gnid={gnid}|spid={spid}\n' + seqstring)

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
../config/genomes.tsv
../clique4+_pcommunity/clique4+_pcommunity.py
    ../clique4+_pcommunity/out/4clique/pclusters.txt
../OG_meta/OG_meta.py
    ../OG_meta/out/OG_meta.tsv
"""