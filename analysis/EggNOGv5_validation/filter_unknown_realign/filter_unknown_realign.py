"""Remove sequences containing unknown amino acids and re-align."""

import gzip
import os
import shutil
import subprocess
from Bio import AlignIO

# Parameters
path_inmembers = '../../../data/EggNOGv5/drosophilidae/7214_members.tsv'
path_outmembers = '7214_noX_members.tsv'
resume = False  # Run ete3 if directory for exists for re-alignment

# Make output directories
if not os.path.exists('align'):
    os.mkdir('align')
if not os.path.exists('raw'):
    os.mkdir('raw')

with open(path_outmembers, 'w') as file_outmembers:
    with open(path_inmembers) as file_inmembers:
        for line in file_inmembers:
            fields = line.rstrip().split('\t')
            alignment_id = fields[1]
            with gzip.open('../../../data/EggNOGv5/drosophilidae/7214/' + alignment_id + '.raw_alg.faa.gz', 'rt') as file_MSA:  # read, text mode
                MSA = AlignIO.read(file_MSA, 'fasta')
                unknown = False  # If alignment contains sequence with unknown amino acids
                proteins = []  # Proteins in filtered alignment
                species = set()  # Species in filtered alignment
                records = []  # Records to generate filtered alignment (for later re-alignment)

                # Iterate through records
                for i, record in enumerate(MSA):
                    if 'X' in record.seq.upper():
                        unknown = True
                    else:
                        proteins.append(record.name)
                        species.add(record.name.split('.')[0])
                        records.append(record)
                MSA = AlignIO.MultipleSeqAlignment(records)  # Re-assign MSA to remove sequences with unknown amino acids

                if unknown and len(MSA) > 1:  # Re-align if unknown flag is True and more than sequence is present
                    args = ['xvfb-run', '/netdata/singlemd/miniconda3/envs/ete3/bin/ete3', 'build',
                            '-w', 'eggnog41', '-a', f'./raw/{alignment_id}/unaligned.fa', '-o', f'./raw/{alignment_id}', '--dealign', '--cpu', '4']
                    path_aligned = f'./raw/{alignment_id}/metaligner_trimmed-trimal01-prottest_default-phyml_default/unaligned.fa.final_tree.fa'

                    # Create directory to store ete3 input and output
                    if not os.path.exists(f'./raw/{alignment_id}'):  # If output directory does not exist, make and execute from scratch
                        os.mkdir(f'./raw/{alignment_id}')
                    elif resume or not os.path.exists(path_aligned):  # If out directory already exists or if final output file does not exist, resume execution
                        args.append('--resume')
                    else:  # Otherwise output exists and resume is false, so do nothing
                        args = ['true']

                    # Save unaligned sequences to file
                    with open(f'./raw/{alignment_id}/unaligned.fa', 'w') as file_unaligned:
                        AlignIO.write(MSA, file_unaligned, 'fasta')

                    # Execute re-alignment in ete3
                    subprocess.run(args, check=True)  # check raises exception if process exits with non-zero code

                    # Save zipped copy of output to alignment directory
                    with open(path_aligned, 'rb') as file_aligned:
                        with gzip.open(f'./align/{alignment_id}.raw_alg.faa.gz', 'wb') as file_zipped:
                            shutil.copyfileobj(file_aligned, file_zipped)
                elif len(MSA) == 1:  # Save sequence without gaps if only one sequence is present
                    # Create directory to store ete3 input and output
                    if not os.path.exists(f'./raw/{alignment_id}'):
                        os.mkdir(f'./raw/{alignment_id}')

                    # Save unaligned sequences to file with gaps removed
                    with open(f'./raw/{alignment_id}/unaligned.fa', 'w') as file_unaligned:
                        file_unaligned.write('>' + MSA[0].name + '\n')  # Header
                        file_unaligned.write(str(record.seq.upper()).translate({ord('-'): None}) + '\n') # Sequence with gaps removed and uppercase

                    # Save zipped copy of output to alignment directory
                    with open(f'./raw/{alignment_id}/unaligned.fa', 'rb') as file_aligned:
                        with gzip.open(f'./align/{alignment_id}.raw_alg.faa.gz', 'wb') as file_zipped:
                            shutil.copyfileobj(file_aligned, file_zipped)
                else:  # Otherwise create hard link to original data in output directory
                    if not os.path.exists(f'align/{alignment_id}.raw_alg.faa.gz'):
                        os.link('../../../data/EggNOGv5/drosophilidae/7214/' + alignment_id + '.raw_alg.faa.gz', f'./align/{alignment_id}.raw_alg.faa.gz')

                # Write updated members to file
                str_genes = ','.join(sorted(proteins))
                str_species = ','.join(sorted(species))
                outmember = '\t'.join([fields[0], fields[1], str(len(proteins)), str(len(species)), str_genes, str_species, str(unknown)])
                file_outmembers.write(outmember + '\n')

"""
NOTES
Must run in bash shell in ete3 conda environment to correctly execute shell commands

DEPENDENCIES
../../../data/EggNOGv5/drosophilidae/
    ../../../data/EggNOGv5/drosophilidae/7214/
"""