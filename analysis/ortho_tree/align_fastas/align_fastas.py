"""Align FASTAs of OGs with one unique sequence per gene."""

import multiprocessing as mp
import os
from subprocess import run, CalledProcessError
from time import time_ns


def run_cmd(file_id):
    aligners = [('probalign', f'../../../bin/probalign ../make_fastas/out/{file_id}.tfa '
                              f'1> out/{file_id}.mfa 2> out/{file_id}.err'),
                ('probcons', f'../../../bin/probcons ../make_fastas/out/{file_id}.tfa '
                             f'1> out/{file_id}.mfa 2> out/{file_id}.err'),
                ('mafft', f'../../../bin/mafft --genafpair --maxiterate 1000 --thread 1 '
                          f'../make_fastas/out/{file_id}.tfa '
                          f'1> out/{file_id}.mfa 2> out/{file_id}.err'),
                ('clustalo', f'../../../bin/clustalo -i ../make_fastas/out/{file_id}.tfa '
                             f'--threads=1 -o out/{file_id}.mfa')]
    for (aligner, cmd) in aligners:
        try:
            t0 = time_ns()
            run(cmd, shell=True, check=True)
            t1 = time_ns()
            return file_id, aligner, str(t1-t0)
        except CalledProcessError:
            pass
    else:
        return file_id, 'ALIGN_ERROR', 'NaN'


num_processes = int(os.environ['SLURM_NTASKS'])

if __name__ == '__main__':
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with mp.Pool(processes=num_processes) as pool:
        file_ids = [file[:-4] for file in os.listdir('../make_fastas/out/')]
        rows = pool.map(run_cmd, file_ids)

    with open('out/times.tsv', 'w') as file:
        file.write('file_id\taligner\ttime_ns\n')
        for row in rows:
            file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../make_fastas/make_fastas.py
    ../make_fastas/out/*.tfa
"""