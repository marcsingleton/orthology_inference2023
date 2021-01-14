"""Extract and count polypeptide IDs, storing associated metadata as tsv."""

import os
import re

pp_regex = {'FlyBase': r'(FBpp[0-9]+)',
            'NCBI': r'([NXY]P_[0-9]+)'}
gn_regex = {'FlyBase': r'parent=(FBgn[0-9]+)',
            'NCBI': r'db_xref=GeneID:([0-9]+)'}

# Parse genomes
genomes = []
with open('../config/genomes.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        spid, _, source, _, tcds_path = line.split()
        genomes.append((spid, source, tcds_path))

# Extract and count polypeptide IDs
ppid_counts = {}  # Counts for each PPID to find duplicates
ppid2meta = {}  # PPID to gene and species
num_headers = 0
for spid, source, tcds_path in genomes:
    with open(tcds_path) as file:
        for line in file:
            if line.startswith('>'):
                num_headers += 1
                gn_match = re.search(gn_regex[source], line)
                pp_match = re.search(pp_regex[source], line)
                try:
                    # First group is entire line, second is first match
                    gnid = gn_match.group(1)
                    ppid = pp_match.group(1)

                    ppid2meta[ppid] = (gnid, spid)
                    ppid_counts[ppid] = ppid_counts.get(ppid, 0) + 1
                except AttributeError:
                    print(line)

# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Write graph as adjacency list to file
with open('out/ppid2meta.tsv', 'w') as outfile:
    for ppid, meta in ppid2meta.items():
        outfile.write(ppid + '\t' + '\t'.join(meta) + '\n')

print('Total headers:', num_headers)
print('Unique IDs:', len(ppid_counts))

"""
OUTPUT
Total headers: 735719
Unique IDs: 735719

DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../config/genomes.tsv
"""