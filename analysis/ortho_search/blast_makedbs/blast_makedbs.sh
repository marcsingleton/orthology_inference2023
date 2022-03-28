#!/bin/zsh
# Construct BLAST databases for the translated Drosophila coding sequences

cd "$(dirname "$0")"
if [ ! -d out/ ]; then
  mkdir out/
fi

cat ../config/genomes.tsv | while read spid txid source prot_path tcds_path
do
  if [[ ${spid} != \#* ]]; then  # Double brackets is expanded syntax for tests
    ../../../bin/ncbi-blast-2.10.1+/bin/makeblastdb -in "../remove_duplicates/out/${spid}.fa" -dbtype prot -title "${txid}_${spid}_blastdb" -parse_seqids -out "out/${spid}_blastdb" -taxid "${txid}" -logfile "out/${spid}_blastdb.out"
  fi
done

# DEPENDENCIES
# ../config/genomes.tsv
# ../remove_duplicates/remove_duplicates.py
# ../remove_duplicates/remove_duplicates/out/*.fa