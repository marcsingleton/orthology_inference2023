"""Extract and count polypeptide IDs, storing associated metadata as tsv."""

import os
import re

pp_regex = {'flybase': r'(FBpp[0-9]+)',
            'ncbi': r'(XP_[0-9]+(\.[0-9]+)?)',
            'YO': r'(YOtr[A-Z0-9]+\|orf[0-9]+)'}
gn_regex = {'flybase': r'parent=(FBgn[0-9]+)',
            'ncbi': r'\[db_xref=GeneID:([0-9]+)\]',
            'YO': r'(YOgn[A-Z]+[0-9]+)'}

# Parse parameters
params = []
with open('params.tsv') as infile:
    fields = infile.readline().split()  # Skip header
    for line in infile:
        params.append(line.split())

# Extract and count polypeptide IDs
ppid_counts = {}  # Counts for each PPID to find duplicates
ppid2meta = {}  # PPID to gene and species
num_headers = 0
for species, _, source, tcds_path in params:
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

                    ppid2meta[ppid] = (gnid, species)
                    ppid_counts[ppid] = ppid_counts.get(ppid, 0) + 1
                except AttributeError:
                    print(line)

# Make output directory
if not os.path.exists(f'out/'):
    os.mkdir(f'out/')

# Write graph as adjacency list to file
with open('out/ppid2meta.tsv', 'w') as outfile:
    for ppid, meta in ppid2meta.items():
        outfile.write(ppid + '\t' + '\t'.join(meta) + '\n')

print('Total headers:', num_headers)
print('Unique IDs:', len(ppid_counts))

"""
OUTPUT
Total headers: 1734870
Unique IDs: 1734870

DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.32_FB2020_01/fasta/dmel-all-translation-r6.32.fasta
../extract_orfs/extract_orfs.py
    ../extract_orfs/out/dpse_translated_orfs.faa
    ../extract_orfs/out/dyak_translated_orfs.faa
./params.tsv
"""