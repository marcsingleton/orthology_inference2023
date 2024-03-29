#!/bin/zsh
# Construct BLAST databases for the translated Drosophila coding sequences

cd "$(dirname "$0")"
if [ ! -d out/ ]; then
  mkdir out/
fi

tail -n +2 ../config/genomes.tsv | while read spid txid source prot_path tcds_path
do
  ../../../bin/ncbi-blast-2.13.0+/bin/makeblastdb -in "../remove_duplicates/out/${spid}.fa" -dbtype prot -title "${txid}_${spid}_blastdb" -parse_seqids -out "out/${spid}_blastdb" -taxid "${txid}" -logfile "out/${spid}_blastdb.out"
done

# DEPENDENCIES
# ../config/genomes.tsv
# ../remove_duplicates/remove_duplicates.py
#   ../remove_duplicates/remove_duplicates/out/*.fa