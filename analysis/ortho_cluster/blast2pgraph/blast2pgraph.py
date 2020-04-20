"""Parse BLAST results as transcript graph."""

import os
import re
from itertools import permutations

pp_regex = {'flybase': r'(FBpp[0-9]+)',
            'ncbi': r'(XP_[0-9]+(\.[0-9]+)?)',
            'YO': r'(YOtr[A-Z0-9]+\|orf[0-9]+)'}

# Parse parameters
params = {}
with open('params.tsv') as infile:
    fields = infile.readline().split()  # Skip header
    for line in infile:
        species, _, source = line.split()
        params[species] = source

# Parse BLAST results
pgraph = {}
for query_species, db_species in permutations(params.keys(), 2):
    with open(f'../blast_AAA/out/{query_species}/{db_species}.blast') as file:
        query_ppid, subjects = None, []
        line = file.readline()
        while line:
            while line.startswith('#'):
                if line == '# BLASTP 2.10.0+\n' and query_ppid is not None:  # Skip for first line
                    pgraph[query_ppid] = pgraph.get(query_ppid, [])
                elif line.startswith('# Query:'):
                    query_ppid = re.search(pp_regex[params[query_species]], line).group(1)
                line = file.readline()

            subjects = []  # Or "hits," but using the BLAST jargon here
            while line and not line.startswith('#'):
                subjects.append(line.split())
                line = file.readline()
            BH_ppid = [re.search(pp_regex[params[db_species]], subjects[0][1]).group(1)] if subjects else []
            try:
                pgraph[query_ppid].extend(BH_ppid)
            except KeyError:
                pgraph[query_ppid] = BH_ppid

# Make output directory
if not os.path.exists(f'out/'):
    os.mkdir(f'out/')

# Write graph as adjacency list to file
with open('out/pgraph.tsv', 'w') as outfile:
    for query_ppid, BH_ppids in pgraph.items():
        outfile.write(query_ppid + '\t' + ','.join(BH_ppids) + '\n')

"""
DEPENDENCIES
../blast_AAA/blast_AAA.py
    ../blast_AAA/out/*
./params.tsv
"""