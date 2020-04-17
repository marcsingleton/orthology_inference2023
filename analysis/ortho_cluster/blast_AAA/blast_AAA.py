"""BLAST all-against-all using annotated Drosophila protein sequences and databases."""

import os
from subprocess import run
from time import asctime

num_threads = '20'

# Parse parameters
params = []
with open('params.tsv') as infile:
    fields = infile.readline().split()  # Skip header
    for line in infile:
        params.append(line.split())

# Execute BLASTs
for query_species, _, tcds_path in params:
    # Make output directory
    if not os.path.exists(f'out/{query_species}'):
        os.makedirs(f'out/{query_species}')

    for db_species, _, _ in params:
        if not os.path.exists(f'out/{query_species}/{db_species}.blast'):
            # Generate args
            input_args = ['../../../bin/ncbi-blast-2.10.0+/bin/blastp', '-query', tcds_path]
            output_args = ['-out', f'out/{query_species}/{db_species}.blast']
            search_args = ['-db', f'../blast_dbs/out/{db_species}_blastdb', '-evalue', '1', '-num_threads', num_threads]
            format_args = ['-outfmt', '7 qacc sacc length pident nident gaps gapsopen qlen qstart qend slen sstart send evalue bitscore']

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
../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.32_FB2020_01/fasta/dmel-all-translation-r6.32.fasta
../extract_orfs/extract_orfs.py
    ../extract_orfs/out/dpse_translated_orfs.faa
    ../extract_orfs/out/dyak_translated_orfs.faa
../blast_dbs/blast_dbs.py
    ../blast_dbs/out/*
./params.tsv
"""