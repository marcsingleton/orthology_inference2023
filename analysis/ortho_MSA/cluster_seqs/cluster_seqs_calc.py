"""Cluster unique sequences in genes."""

import os
import re
from subprocess import run

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
ppid2gnid = {}
with open('../../ortho_search/seq_meta/out/seq_meta.tsv') as file:
    for line in file:
        ppid, gnid, _ = line.split()
        ppid2gnid[ppid] = gnid

# Load seqs
gnid2seqs = {}
repr2cons = {}  # Mapping of representative PPIDs to constituent sequences
for _, source, prot_path in genomes:
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
                for ppid1, seq1 in gnid2seqs[gnid]:
                    if seq0 == seq1:
                        repr2cons[ppid1].append(ppid0)
                        break
                else:
                    gnid2seqs[gnid].append((ppid0, seq0))
                    repr2cons[ppid0] = []
            except KeyError:
                gnid2seqs[gnid] = [(ppid0, seq0)]
                repr2cons[ppid0] = []

# Save clusters
if not os.path.exists('out/'):
    os.mkdir('out/')

# Cluster sequences
clusters = []
for gnid, seqs in gnid2seqs.items():
    if len(seqs) > 1:
        with open(f'out/{gnid}.fasta', 'w') as file:
            for ppid, seq in seqs:
                file.write(f'>{ppid}\n')
                file.write('\n'.join([seq[i:i+60] for i in range(0, len(seq), 60)]) + '\n')
        cmd = ['../../../bin/cd-hit', '-i', f'out/{gnid}.fasta', '-o', f'out/{gnid}',
               '-c', '0.9',  # sequence identity threshold
               '-n', '5',  # word_length (default as recommended for high identity)
               '-d', '0',  # full header up to first space in .clstr file description
               '-g', '1']  # accurate mode: find best matching cluster
        run(cmd, check=True, capture_output=True)

        with open(f'out/{gnid}.clstr') as file:
            line = file.readline()
            while line:
                if line.startswith('>'):
                    line = file.readline()

                cluster = []
                while line and not line.startswith('>'):
                    fields = line.split()
                    ppid = fields[2][1:-3]
                    if len(fields) == 5:
                        pident = fields[4][:-1]
                        cluster.append((ppid, pident))
                    else:
                        ppid0 = ppid
                    line = file.readline()
                clusters.append((ppid0, gnid, cluster))

        os.remove(f'out/{gnid}.fasta')  # Original sequences
        os.remove(f'out/{gnid}')  # Cluster sequences
        os.remove(f'out/{gnid}.clstr')  # Clusters
    else:
        ppid, _ = seqs[0]
        clusters.append((ppid, gnid, []))

with open('out/clusters.tsv', 'w') as file:
    file.write('ppid\tgnid\trepr_id\tcluster_id\tpident\n')  # Identical and cluster PPIDs
    for cluster_id, gnid, ppids in clusters:
        file.write('\t'.join([cluster_id, gnid, cluster_id, cluster_id, '100.00']) + '\n')  # Entry for representative of cluster
        for repr_id, pident in ppids:
            file.write('\t'.join([repr_id, gnid, repr_id, cluster_id, str(pident)]) + '\n')  # Entry for representative of identicals
            for ppid in repr2cons[repr_id]:
                file.write('\t'.join([ppid, gnid, repr_id, cluster_id, str(pident)]) + '\n')  # Entry for identical

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../config/genomes.tsv
../../ortho_search/seq_meta/seq_meta.py
    ../../ortho_search/seq_meta/out/seq_meta.tsv
"""