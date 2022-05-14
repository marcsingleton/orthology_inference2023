"""Make FASTAs of remaining OGs after removing multiple representatives per gene."""

import os
import re

import pandas as pd
from src.utils import read_fasta

ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+(\.[0-9]+)?)'}

# Load genomes
genomes = []
with open('../../ortho_cluster2/config/genomes.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        genomes.append((fields['spid'], fields['source'], fields['prot_path']))

# Load sequence data
ppid2data = {}
with open('../../ortho_search/sequence_data/out/sequence_data.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppid2data[fields['ppid']] = (fields['gnid'], fields['spid'])

# Load seqs
ppid2seq = {}
for spid, source, prot_path in genomes:
    fasta = read_fasta(prot_path)
    for header, seq in fasta:
        ppid = re.search(ppid_regex[source], header).group(1)
        ppid2seq[ppid] = seq

# Load OGs and OG data
OGs = {}
with open('../reduce_tree/out/clusters.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        OGs[fields['OGid']] = fields['ppids'].split(',')
OG_data = pd.read_table('../OG_data/out/OG_data.tsv')

# Write sequences
if not os.path.exists('out/'):
    os.mkdir('out/')

OGids = OG_data.loc[~(OG_data['ppidnum'] == OG_data['gnidnum']), 'OGid']
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
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.45_FB2022_02/fasta/dmel-all-translation-r6.45.fasta
../../ortho_cluster2/config/genomes.tsv
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
../OG_data/OG_data.py
    ../OG_data/out/OG_data.tsv
../reduce_tree/reduce_tree.py
    ../reduce_tree/out/clusters.tsv
"""