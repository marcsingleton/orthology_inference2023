#!/bin/zsh
# Extract sequence IDs from dsim NCBI and FlyBase genomes and compare

cd "$(dirname "$0")"
if [ ! -d out ]; then
  mkdir out
fi

zgrep ">" ../../../data/flybase_genomes/Drosophila_simulans/dsim_r2.02_FB2017_04/fasta/dsim-all-chromosome-r2.02.fasta.gz | cut -f 4 -d ';' | cut -c 12- | sort > out/dsim_r2.02_FB2017_04.out
zgrep ">" ../../../data/ncbi_genomes/fna/GCA_000754195.3_ASM75419v3_genomic.fna.gz | cut -f 1 -d '.' | cut -c 2- | sort > out/GCA_000754195.3_ASM75419v3_genomic.out

diff dsim_r2.02_FB2017_04.out GCA_000754195.3_ASM75419v3_genomic.out

# NOTES
# The only difference between the files is an entry corresponding to the mitochondrial genome

# DEPENDENCIES
# ../../../data/flybase_genomes/Drosophila_simulans/dsim_r2.02_FB2017_04/fasta/dsim-all-chromosome-r2.02.fasta.gz
# ../../../data/ncbi_genomes/fna/GCA_000754195.3_ASM75419v3_genomic.fna.gz