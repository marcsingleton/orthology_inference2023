"""Convert amino acid alignments to nucleotide alignments."""

import os
import re
import sys

from src.utils import read_fasta


def get_codons(nt_seq, aa_seq, file_id):
    delta = len(nt_seq) // 3 - len(aa_seq)
    if abs(delta) > 1:
        print('File:', file_id)
        print('PPID:', ppid)
        print(f'Error: Lengths of translated genomic ({len(nt_seq) // 3}) '
              f'and protein ({len(aa_seq)}) sequences are incompatible.')
        print()
        return
    num_codon = (len(nt_seq) - 3 * max(0, delta)) // 3  # Exclude last codon only if delta == 1

    codons = []  # Translated symbols
    for i in range(num_codon):
        codon = nt_seq[3*i:3*i+3]
        tr = 'X' if 'N' in codon else codon2aa[codon]
        aa = aa_seq[i]
        if tr != aa:
            if i == 0 and aa == 'M':  # Correct non-standard start codon
                codons.append(codon)
                print('File:', file_id)
                print('PPID:', ppid)
                print('Warning: Non-standard start codon detected. '
                      f'Translated residue M{i+1} converted to observed residue {aa}.')
                print()
            elif tr == 'X':  # Correct unknown amino acid
                codons.append(codon)
                print('File:', file_id)
                print('PPID:', ppid)
                print(f'Warning: Converted unknown translated residue X{i+1} to observed residue {aa}.')
                print()
            elif tr == '*':
                codons.append(codon)
                print('File:', file_id)
                print('PPID:', ppid)
                print(f'Warning: Converted stop codon *{i+1} to observed residue {aa}.')
                print()
            else:  # Mismatch
                print('File:', file_id)
                print('PPID:', ppid)
                print(f'Error: Mismatched translated residue {tr}{i+1} and observed residue {aa}.')
                print()
                return
        else:
            codons.append(codon)
    if delta == -1:
        codons.append(nt_seq[3*num_codon:] + (3-len(nt_seq) % 3) * '-')
        print('File:', file_id)
        print('PPID:', ppid)
        print(f'Warning: Missing final residue {aa_seq[-1]} appended to translated sequence.')
        print()
    return codons


ppid_regex = {'FlyBase': r'(FBpp[0-9]+)',
              'NCBI': r'([NXY]P_[0-9]+)'}

# Load genomes
genomes = []
with open('../config/genomes.tsv') as file:
    file.readline()  # Skip header
    for line in file:
        spid, _, source, _, cds_path = line.rstrip('\n').split('\t')
        genomes.append((spid, source, cds_path))

# Load codon table
codon2aa = {}
with open('codons.txt') as file:
    lines = [line.rstrip('\n').split(' = ')[1] for line in file]
    for i in range(len(lines[0])):
        aa = lines[0][i]
        codon = ''.join([lines[j][i] for j in range(2, 5)])
        codon2aa[codon] = aa

# Load CDSs
ppid2cds = {}
for spid, source, cds_path in genomes:
    fasta = read_fasta(cds_path)
    for header, seq in fasta:
        ppid = re.search(ppid_regex[header], line).group(1)
        ppid2cds[ppid] = seq

if not os.path.exists('out/'):
    os.mkdir('out/')

sys.stdout = open('out/out.txt', 'w')  # Redirect stdout to file
for path in [path for path in os.listdir('../align_fastas/out/') if path.endswith('.afa')]:
    file_id = path.removesuffix('.afa')

    # Translate and write CDS
    nt_aligns = []
    for header, aa_align in read_fasta('../align_fastas/out/' + path):
        ppid = re.search(r'ppid=([A-Za-z0-9_]+)', header).group(1)
        aa_seq = aa_align.replace('-', '')
        nt_seq = ppid2cds[ppid]

        codons = get_codons(nt_seq, aa_seq, file_id)
        if not codons:
            print(f'Conversion of {path} failed due to {ppid}.')
            print()
            break

        nt_align = []
        for aa in aa_align:
            if aa == '-':
                nt_align.append('---')
            else:
                nt_align.append(codons.pop(0))
        nt_aligns.append((header, ''.join(nt_align)))
    else:
        with open(f'out/{file_id}.afa', 'w') as file:
            for header, nt_align in nt_aligns:
                seqstring = '\n'.join(nt_align[i:i+80] for i in range(0, len(nt_align), 80))
                file.write(f'{header}\n{seqstring}\n')
sys.stdout.close()

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_cds_from_genomic.fna
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.38_FB2021_01/fasta/dmel-all-CDS-r6.38.fasta
../config/genomes.tsv
../align_fastas/align_fastas.py
    ../align_fastas/out/*.afa
./codons.txt
"""