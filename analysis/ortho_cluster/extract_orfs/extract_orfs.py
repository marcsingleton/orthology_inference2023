"""Extract and translate ORFs from transcripts."""

import Bio.SeqIO as SeqIO
import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def find_orfs(seq, trans_table=1, min_protein_length=1):
    orfs = []
    for frame in range(3):
        trans = str(seq[frame:len(seq) - (len(seq) - frame) % 3].translate(trans_table).seq)
        aa_start = trans.find('M')
        stop = True
        while aa_start != -1:
            aa_end = trans.find('*', aa_start)
            if aa_end == -1:
                aa_end = len(trans)
                stop = False

            if aa_end - aa_start >= min_protein_length:
                nt_start = 3 * aa_start + frame
                nt_end = 3 * (aa_end + int(stop)) + frame  # Add extra codon to end if stop codon exists
                orfs.append((nt_start, nt_end, stop, trans[aa_start:aa_end]))
            aa_start = trans.find('M', aa_end + 1)
    return orfs


# Make output directory
if not os.path.exists('out/'):
    os.mkdir('out/')

# Extract and translate ORFs
for path in filter(lambda x: x.endswith('.fna'), os.listdir('../extract_transcripts/out/')):
    trans = []
    for record in SeqIO.parse('../extract_transcripts/out/' + path, 'fasta'):
        orfs = find_orfs(record)
        ids = record.description.split()
        for i, orf in enumerate(orfs):
            trans.append(SeqRecord(Seq(orf[3]), id=f'{ids[1][7:]}|{ids[0]}|orf{i}', description=f'nt_start={orf[0]} nt_end={orf[1]} stop={orf[2]}'))
    SeqIO.write(trans, f'out/{path[0:4]}_translated_orfs.faa', 'fasta')

"""
DEPENDENCIES
../extract_transcripts/extract_transcripts.sh
    ../extract_transcripts/out/*
"""