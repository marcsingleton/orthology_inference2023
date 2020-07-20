#!/bin/zsh
# Extract transcripts from YO annotations

cd "$(dirname "$0")"
if [ ! -d out ]; then
  mkdir out
fi

cat params.tsv | while read species genome_path
do
	if [[ $species != \#* ]]  # Double brackets is expanded syntax for tests
	then
		../../../bin/gffread -w out/"$species"_transcript.fna -g ../../../data/flybase_genomes/"$genome_path" -F ../../../data/YO_annotations/"$species".YO.gff3
	fi
done

# DEPENDENCIES
# ../../../data/flybase_genomes/*
# ../../../data/YO_annotations/*.YO.gff3
# ./params.tsv