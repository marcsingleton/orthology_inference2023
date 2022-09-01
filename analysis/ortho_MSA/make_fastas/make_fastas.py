"""Make FASTAs of OGs with one unique sequence per gene."""

import os
import re

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
with open('../../ortho_cluster2/add_paralogs/out/clusters.tsv') as file:
    field_names = file.readline().rstrip('\n').split('\t')
    for line in file:
        fields = {key: value for key, value in zip(field_names, line.rstrip('\n').split('\t'))}
        ppids = {node for edge in fields['edges'].split(',') for node in edge.split(':')}
        OGs[fields['OGid']] = ppids

# Write sequences
if not os.path.exists('out/'):
    os.mkdir('out/')

for OGid, OG in OGs.items():
    records = []
    for ppid in OG:
        gnid, spid = ppid2data[ppid]
        seq = ppid2seq[ppid]
        records.append({'ppid': ppid, 'gnid': gnid, 'spid': spid, 'seq': seq})
    with open(f'out/{OGid}.fa', 'w') as file:
        for record in sorted(records, key=lambda x: (x['spid'], x['gnid'], x['ppid'])):
            ppid, gnid, spid, seq = record['ppid'], record['gnid'], record['spid'], record['seq']
            seqstring = '\n'.join([seq[i:i + 80] for i in range(0, len(seq), 80)])
            file.write(f'>ppid={ppid}|gnid={gnid}|spid={spid}\n{seqstring}\n')

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.45_FB2022_02/fasta/dmel-all-translation-r6.45.fasta
../../ortho_cluster2/config/genomes.tsv
../../ortho_cluster2/add_paralogs/add_paralogs.py
    ../../ortho_cluster2/add_paralogs/out/clusters.tsv
../../ortho_search/sequence_data/sequence_data.py
    ../../ortho_search/sequence_data/out/sequence_data.tsv
"""