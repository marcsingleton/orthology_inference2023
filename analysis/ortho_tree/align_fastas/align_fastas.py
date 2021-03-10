"""Align FASTAs of OGs with one unique sequence per gene."""

import multiprocessing as mp
import os
from subprocess import run, CalledProcessError
from time import time_ns


def run_cmd(file_id):
    cmd = (f'../../../bin/mafft --globalpair --maxiterate 1000 '
           f'--thread 1 --anysymbol --allowshift --unalignlevel 0.4 --leavegappyregion '
           f'../make_fastas/out/{file_id}.tfa '
           f'1> out/{file_id}.mfa 2> out/{file_id}.err')
    try:
        t0 = time_ns()
        run(cmd, shell=True, check=True)
        t1 = time_ns()
        return file_id, str(t1-t0)
    except CalledProcessError:
        return file_id, 'NaN'


num_processes = int(os.environ['SLURM_NTASKS'])

if __name__ == '__main__':
    if not os.path.exists('out/'):
        os.mkdir('out/')

    with mp.Pool(processes=num_processes) as pool:
        file_ids = [file[:-4] for file in os.listdir('../make_fastas/out/')]
        rows = pool.map(run_cmd, file_ids)

    with open('out/times.tsv', 'w') as file:
        file.write('file_id\ttime_ns\n')
        for row in rows:
            file.write('\t'.join(row) + '\n')

"""
DEPENDENCIES
../make_fastas/make_fastas.py
    ../make_fastas/out/*.tfa
"""