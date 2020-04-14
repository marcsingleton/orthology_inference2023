#!/bin/zsh
# Count sequences in NCBI and FlyBase genomes

cd "$(dirname "$0")"
if [ ! -d out ]; then
  mkdir out
fi

cat fb2ncbi.tsv | while read species ncbi flybase
do
	if [[ $species != \#* ]]; then  # Double brackets is expanded syntax for tests
		echo "$species"
		zgrep ">" ../../../data/ncbi_genomes/fna/"$ncbi" | wc
		zgrep ">" ../../../data/flybase_genomes/"$flybase" | wc
	fi
done

# NOTES
# The identical counts between files indicates that despite differences in names, the genomes are identical.
# The additional sequence in the flybase dsim genome is a mitochondrial sequence as shown in dsim_seqs.sh

# DEPENDENCIES
# ../../../data/flybase_genomes/*
# ../../../data/ncbi_genomes/fna/*
# ./fb2ncbi.tsv