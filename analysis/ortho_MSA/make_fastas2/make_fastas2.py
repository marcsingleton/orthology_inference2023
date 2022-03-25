"""Make FASTAs of remaining OGs using trees."""

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
        spid, _, source, prot_path = line.split()
        genomes.append((spid, source, prot_path))

# Load sequence data
ppid2data = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        ppid, gnid, spid, sqid = line.split()
        ppid2data[ppid] = (gnid, spid, sqid)

# Load seqs
ppid2seq = {}
for spid, source, prot_path in genomes:
    fasta = read_fasta(prot_path)
    for header, seq in fasta:
        ppid = re.search(ppid_regex[header], line).group(1)
        ppid2seq[ppid] = seq

# Load OGs and OG data
OGs = {}
with open('../reduce_tree/out/clusters.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        OGid, ppids = line.split()
        OGs[OGid] = ppids.split(',')
OG_data = pd.read_table('../OG_data/out/OG_data.tsv')

# Write sequences
if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = OG_data.loc[~(OG_data['sqidnum'] == OG_data['gnidnum']), 'OGid']
for OGid in OGids:
    records = []
    for sqid in OGs[OGid]:
        gnid, spid, _ = ppid2data[sqid]
        seq = ppid2seq[sqid]
        seqstring = '\n'.join([seq[i:i+80] for i in range(0, len(seq), 80)]) + '\n'
        records.append((sqid, gnid, spid, seqstring))
    with open(f'out/{OGid}.fa', 'w') as file:
        for sqid, gnid, spid, seqstring in sorted(records, key=lambda x: x[2]):
            file.write(f'>ppid={sqid}|gnid={gnid}|spid={spid}\n' + seqstring)

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/fasta/dmel-all-translation-r6.38.fasta
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../config/genomes.tsv
../OG_data/OG_data.py
    ../OG_data/out/OG_data.tsv
../reduce_tree/reduce_tree.py
    ../reduce_tree/out/clusters.tsv
"""