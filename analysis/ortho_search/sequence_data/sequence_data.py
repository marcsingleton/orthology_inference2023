"""Extract protein and gene IDs from FASTAs, storing data as tsv."""

import os
import re

from src.utils import read_fasta


def get_gnid(header, source):
    gnid = None
    if source == 'FlyBase':
        match = re.search(r'parent=(FBgn[0-9]+)', header)
        if match:
            gnid = match.group(1)
    elif source == 'NCBI':
        match = re.search(r'\[db_xref=([^]]+)]', header)
        if match:
            db_xrefs = match.group(1).split(',')
            gnids = [db_xref.split(':')[1] for db_xref in db_xrefs]
            gnid = gnids[0]  # Use first listed GNID if multiple exist
    return gnid


def get_ppid(header, source):
    ppid = None
    match = re.search(tcds_regex[source], header)
    if match:
        ppid = match.group(1)
    return ppid


prot_regex = {'FlyBase': r'>(FBpp[0-9]+)',
              'NCBI': r'>([NXY]P_[0-9]+(\.[0-9]+)?)'}
tcds_regex = {'FlyBase': r'ID=(FBpp[0-9]+)',
              'NCBI': r'protein_id=([NXY]P_[0-9]+(\.[0-9]+)?)'}

# Load genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, prot_path, tcds_path = line.rstrip('\n').split('\t')
        genomes.append((spid, source, prot_path, tcds_path))

# Extract and count IDs
ppid_counts = {}
ppid2data = {}
gnid2seqs = {}
for spid, source, prot_path, tcds_path in genomes:
    # Find parent genes in tcds headers
    tcds_fasta = read_fasta(tcds_path)
    for header, _ in tcds_fasta:
        gnid = get_gnid(header, source)
        ppid = get_ppid(header, source)
        if gnid and ppid:
            ppid2data[ppid] = (gnid, spid)
            ppid_counts[ppid] = ppid_counts.get(ppid, 0) + 1
        else:
            print(f'regex failure in {spid} tcds_path: {header}')

    # Find representative sequences in prot files
    prot_fasta = read_fasta(prot_path)
    for header, seq in prot_fasta:
        ppid_match = re.search(prot_regex[source], header)
        if ppid_match:
            ppid = ppid_match.group(1)
        else:
            print(f'regex failure in {spid} tcds_path: {header}')
            continue

        gnid, spid = ppid2data[ppid]
        if gnid in gnid2seqs:
            seqs = gnid2seqs[gnid]
            if seq in seqs:
                ppid2data[ppid] = (gnid, spid, seqs[seq])
            else:
                seqs[seq] = ppid
                ppid2data[ppid] = (gnid, spid, ppid)
        else:
            gnid2seqs[gnid] = {seq: ppid}
            ppid2data[ppid] = (gnid, spid, ppid)

# Write to file
if not os.path.exists('out/'):
    os.mkdir('out/')

with open('out/sequence_data.tsv', 'w') as file:
    file.write('ppid\tgnid\tspid\tsqid\n')
    for ppid, data in ppid2data.items():
        file.write(ppid + '\t' + '\t'.join(data) + '\n')

output = f"""\
Total headers: {sum(ppid_counts.values())}
Unique PPIDs: {len(ppid_counts)}
"""
with open('out/output.txt', 'w') as file:
    file.write(output)

"""
NOTES
Extracting the IDs from the headers is somewhat complex since a few sequences in the dyak annotation have multiple GNIDs
associated with them. The format of the db_xrefs makes a regex approach insufficiently flexible, so instead a more
manual parsing approach is used where the individual db_xrefs are split on their delimiters.

The PPIDs in the protein files are given as the field immediately following the >, so these are extracted directly with
a regex just to prevent any shenanigans. Even though for dmel the prot and tcds files are actually the same, the PPIDs
are extracted differently in the two cases for consistency with the approach for the NCBI files.

DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_protein.faa
../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.45_FB2022_02/fasta/dmel-all-translation-r6.45.fasta
../config/genomes.tsv
"""