"""Get PPIDs of the representative polypeptides for each gene."""

import os
import re

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
ppid2gnid = {}
with open('../../ortho_cluster2/ppid2meta/out/ppid2meta.tsv') as file:
    for line in file:
        ppid, gnid, _ = line.split()
        ppid2gnid[ppid] = gnid

# Load seqs
gnid2seqs = {}
for _, source, prot_path in params:
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
                        break
                else:
                    gnid2seqs[gnid].append((ppid0, seq0))
            except KeyError:
                gnid2seqs[gnid] = [(ppid0, seq0)]

if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/reprs.tsv', 'w') as file:
    for gnid, seqs in gnid2seqs.items():
        ppidstring = ','.join([ppid for ppid, _ in seqs])
        file.write(f'{gnid}\t{ppidstring}\n')

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../../ortho_cluster2/ppid2meta/ppid2meta.py
    ../../ortho_cluster2/ppid2meta/out/ppid2meta.tsv
./params.tsv
"""