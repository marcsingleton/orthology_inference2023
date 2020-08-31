"""BLAST all-against-all using annotated Drosophila protein sequences and databases."""

import os
from sys import argv
from subprocess import run
from time import asctime

num_threads = str(os.environ['SLURM_CPUS_ON_NODE'])
query_species = argv[1]
tcds_path = argv[2]
blast_path = argv[3]

# Parse parameters
params = []
with open('../blast_dbs/params.tsv') as file:
    fields = file.readline().split()  # Skip header
    for line in file:
        params.append(line.split())

# Make output directory
if not os.path.exists(f'out/{query_species}/'):
    os.makedirs(f'out/{query_species}/')  # Recursive folder creation

# Execute BLASTs
for db_species, _, _ in params:
    if not os.path.exists(f'out/{query_species}/{db_species}.blast'):
        # Generate args
        input_args = [blast_path, '-query', tcds_path]
        output_args = ['-out', f'out/{query_species}/{db_species}.blast']
        search_args = ['-db', f'../blast_dbs/out/{db_species}_blastdb', '-evalue', '1', '-num_threads', num_threads]
        format_args = ['-outfmt', '7 qacc sacc length nident gaps qlen qstart qend slen sstart send evalue bitscore']

        # Execute command
        t0 = asctime()
        run(input_args + output_args + search_args + format_args, check=True)
        t1 = asctime()

        # Manually write output to file since direction while in background does not immediately write to file
        with open('out/blast_AAA.out', 'a') as outfile:
            outfile.write('\t'.join([query_species, db_species, t0, t1 + '\n']))

"""
DEPENDENCIES
../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
../blast_dbs/blast_dbs.py
    ../blast_dbs/out/*
./params.tsv
"""