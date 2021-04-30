#!/bin/zsh
# Construct BLAST databases for the translated Drosophila coding sequences

cd "$(dirname "$0")"
if [ ! -d out/ ]; then
  mkdir out/
fi

cat ../config/genomes.tsv | while read spid txid source prot_path tcds_path
do
  if [[ ${spid} != \#* ]]; then  # Double brackets is expanded syntax for tests
    ../../../bin/ncbi-blast-2.10.1+/bin/makeblastdb -in "${prot_path}" -dbtype prot -title "${txid}_${spid}_blastdb" -parse_seqids -out "out/${spid}_blastdb" -taxid "${txid}" -logfile "out/${spid}_blastdb.out"
  fi
done

# DEPENDENCIES
# ../../../data/ncbi_annotations/*/*/*/*_protein.faa
# ../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.34_FB2020_03/fasta/dmel-all-translation-r6.34.fasta
# ../config/genomes.tsv