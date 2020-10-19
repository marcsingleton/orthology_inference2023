"""Make FASTAs of OGs with one unique sequence per gene."""

import os
import re

import pandas as pd

pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+(\.[0-9]+)?)'}
gn_regex = {'FlyBase': r'parent=(FBgn[0-9]+)',
            'NCBI': r'db_xref=GeneID:([0-9]+)'}

# Parse parameters
params = []
with open('params.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, tcds_path = line.split()
        params.append((spid, source, tcds_path))

# Load seqs
gnid2seqs = {}
for spid, source, tcds_path in params:
    with open(tcds_path) as file:
        line = file.readline()
        while line:
            if line.startswith('>'):
                ppid = re.search(pp_regex[source], line).group(1)
                gnid = re.search(gn_regex[source], line).group(1)
                line = file.readline()

            seqlines = []
            while line and not line.startswith('>'):
                seqlines.append(line.rstrip())
                line = file.readline()
            seq0 = ''.join(seqlines)

            try:
                for _, _, seq in gnid2seqs[gnid]:
                    if seq0 == seq:
                        break
                else:
                    gnid2seqs[gnid].append((ppid, spid, seq0))
            except KeyError:
                gnid2seqs[gnid] = [(ppid, spid, seq0)]

# Load OGs and OG metadata
OGs = {}
with open('../../ortho_cluster3/clique5+_community/out/ggraph2/6clique/gclusters.txt') as file:
    for line in file:
        _, OGid, edges = line.rstrip().split(':')
        gnids = set([node for edge in edges.split('\t') for node in edge.split(',')])
        OGs[OGid] = gnids

OGs_meta = pd.read_table('../OGid2meta/out/OGid2meta.tsv')

# Write sequences
gn25 = OGs_meta['gnidnum'] == 25
sp25 = OGs_meta['spidnum'] == 25
sq25 = OGs_meta['sqidnum'] == 25
OGids = OGs_meta.loc[gn25 & sp25 & sq25, 'OGid']

if not os.path.exists('out/'):
    os.mkdir('out/')

for OGid in OGids:
    with open(f'out/{OGid}.tfa', 'w') as file:
        for gnid in OGs[OGid]:
            ppid, spid, seq = gnid2seqs[gnid][0]
            seqstring = '\n'.join(seq[i:i+80] for i in range(0, len(seq), 80)) + '\n'
            file.write(f'>ppid={ppid}|gnid={gnid}|spid={spid}\n')
            file.write(seqstring)

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../../ortho_cluster3/clique5+_community/clique5+_community.py
    ../../ortho_cluster3/clique5+_community/out/ggraph2/6clique/gclusters.txt
../OGid2meta/OGid2meta.py
    ../OGid2meta/out/OGid2meta.tsv
./params.tsv
"""