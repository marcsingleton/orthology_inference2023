#!/bin/zsh
# Construct BLAST databases for the translated Drosophila coding sequences

cd "$(dirname "$0")"
if [ ! -d out ]; then
  mkdir out
fi

cat params.tsv | while read species taxid tcds_path;
do
  if [[ $species != \#* ]]; then  # Double brackets is expanded syntax for tests
    ../../../bin/ncbi-blast-2.10.0+/bin/makeblastdb -in "$tcds_path" -dbtype prot -title "${taxid}_${species}_blastdb" -parse_seqids -out "out/${species}_blastdb" -taxid "$taxid" -logfile "out/${species}_blastdb.out"
  fi
done

# DEPENDENCIES
# ../../../data/ncbi_annotations/*/*/*/*_translated_cds.faa
# ../../../data/flybase_genomes/Drosophila_melanogaster/dmel_r6.32_FB2020_01/fasta/dmel-all-translation-r6.32.fasta
# ../extract_orfs/extract_orfs.py
#   ../extract_orfs/out/dpse_translated_orfs.faa
#   ../extract_orfs/out/dyak_translated_orfs.faa
# ./params.tsv