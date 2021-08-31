#!/bin/zsh
# Count sequences in protein and translated CDS files of NCBI genomes

cd "$(dirname "$0")"
if [ ! -d out/ ]; then
  mkdir out/
fi

if [ -f out/count_seqs.out ]; then
  rm out/count_seqs.out
fi

cat ../config/genomes.tsv | while read spid txid source prot_path tcds_path
do
	if [[ ${spid} != \#* && ${spid} != dmel ]]; then  # Double brackets is expanded syntax for tests
		echo "${spid}" >> out/count_seqs.out
		grep ">" ${prot_path} | wc >> out/count_seqs.out
		grep ">" ${tcds_path} | wc >> out/count_seqs.out
	fi
done

# DEPENDENCIES
# ../../../data/ncbi_annotations/*/*/*/*_protein.faa
# ../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
# ../config/genomes.tsv