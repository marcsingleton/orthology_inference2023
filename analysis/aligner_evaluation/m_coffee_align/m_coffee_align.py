"""Align BAliBase sequences using m-coffee."""

import multiprocessing as mp
import os
import re
from subprocess import run
from time import time_ns


def get_cmd(aln_paths, label):
    cmd = (f'{bin_path}t_coffee {data_path}RV{ref}/{file_id}.tfa '
           f'-aln {aln_paths} '
           f'-thread=1 -newtree=out/m_coffee{label}/RV{ref}/{file_id}.dnd '
           f'-output=fasta_aln -outfile=out/m_coffee{label}/RV{ref}/{file_id}.mfa '
           f'&> out/m_coffee{label}/RV{ref}/{file_id}.err')
    return cmd


def run_cmd(aligner, ref, file_id, cmd):
    t0 = time_ns()
    run(cmd, shell=True, check=True)
    t1 = time_ns()

    return aligner, ref, file_id, str(t1-t0)


num_processes = int(os.environ['SLURM_CPUS_ON_NODE'])
bin_path = '../../../bin/'
data_path = '../../../data/bb3_release/'

if __name__ == '__main__':
    # Make list of commands
    cmds = []
    for ref in ['11', '12', '20', '30', '40', '50']:
        regex = f'(BB{ref}' + r'[0-9]{3}).tfa'
        for file in os.listdir(f'../../../data/bb3_release/RV{ref}/'):
            match = re.match(regex, file)
            if match:
                file_id = match.group(1)
                aligners = ['clustalo1', 'mafft', 'probalign', 'probcons', 't_coffee']
                aln_paths1 = ' '.join([f'../bm_align/out/{aligner}/RV{ref}/{file_id}.mfa' for aligner in aligners])
                aln_paths2 = ' '.join([f'../bm_align/out/{aligner}/RV{ref}/{file_id}.mfa' for aligner in aligners[:-1]])

                cmd1 = get_cmd(aln_paths1, 1)
                cmd2 = get_cmd(aln_paths2, 2)

                if not os.path.exists(f'out/m_coffee1/RV{ref}/'):
                    os.makedirs(f'out/m_coffee1/RV{ref}/')
                if not os.path.exists(f'out/m_coffee2/RV{ref}/'):
                    os.makedirs(f'out/m_coffee2/RV{ref}/')

                cmds.append(('m_coffee1', ref, file_id, cmd1))
                cmds.append(('m_coffee2', ref, file_id, cmd2))

    # Run commands
    with mp.Pool(processes=num_processes) as pool:
        rows = pool.starmap(run_cmd, cmds)

    # Save results
    with open('out/times.tsv', 'w') as file:
        file.write('aligner\tref\tfile_id\ttime_ns\n')
        for row in rows:
            file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../../../data/bb3_release/
../bm_align/bm_align.py
    ../bm_align/out/*
"""