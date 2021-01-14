#!/bin/zsh
# Extract transcripts from YO annotations

cd "$(dirname "$0")"
if [ ! -d out/ ]; then
  mkdir out/
fi

cat genomes.tsv | while read spid genome_path
do
	if [[ ${spid} != \#* ]]  # Double brackets is expanded syntax for tests
	then
		../../../bin/gffread -w out/"${spid}"_transcript.fna -g ../../../data/flybase_genomes/"${genome_path}" -F ../../../data/YO_annotations/"${spid}".YO.gff3
	fi
done

# DEPENDENCIES
# ../../../data/flybase_genomes/*
# ../../../data/YO_annotations/*.YO.gff3
# ./genomes.tsv