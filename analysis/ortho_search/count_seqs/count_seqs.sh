#!/bin/zsh
# Count sequences in protein and translated CDS files of NCBI genomes

cd "$(dirname "$0")"
if [ ! -d out ]; then
  mkdir out
fi

if [ -f out/count_seqs.out ]; then
  rm out/count_seqs.out
fi

cat params.tsv | while read species taxid
do
	if [[ $species != \#* && $species != dmel ]]; then;  # Double brackets is expanded syntax for tests
		echo "$species" >> out/count_seqs.out
		grep ">" $(find ../../../data/ncbi_annotations/"${taxid}" -type f -name "*protein.faa") | wc >> out/count_seqs.out
		grep ">" $(find ../../../data/ncbi_annotations/${taxid} -type f -name "*translated_cds.faa") | wc >> out/count_seqs.out
	fi
done

# DEPENDENCIES
# ../../../data/flybase_genomes/*
# ../../../data/ncbi_genomes/fna/*
# ./params.tsv