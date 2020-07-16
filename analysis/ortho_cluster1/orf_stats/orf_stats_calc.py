"""Calculate statistics for ORFs in transcripts from YO annotations."""

from Bio import SeqIO
import os
import pandas as pd


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

for path in filter(lambda x: x.endswith('.tsv'), os.listdir('../extract_transcripts/out/')):
    max_orfs = []
    for record in SeqIO.parse('../extract_transcripts/' + path, 'fasta'):
        orfs = sorted(find_orfs(record), key=lambda x: x[1] - x[0], reverse=True)

        # Attempt to extract longest two transcripts
        try:
            orf0 = orfs[0]
            len0 = orf0[1] - orf0[0]
            stop0 = orf0[2]
        except IndexError:
            len0 = stop0 = 'NA'
        try:
            orf1 = orfs[1]
            len1 = orf1[1] - orf1[0]
            stop1 = orf1[2]
        except IndexError:
            len1 = stop1 = 'NA'

        # Compile stats
        max_orfs.append({'name': record.name, 'len': len(record),
                         'len0': len0, 'stop0': stop0,
                         'len1': len1, 'stop1': stop1})

    # Convert to dataframe and save
    pd.DataFrame(max_orfs).to_csv(f'out/{path[:4]}_maxorfs.tsv', sep='\t', index=False)

"""
DEPENDENCIES
../extract_transcripts/extract_transcripts.sh
    ../extract_transcripts/out/*.fna
"""