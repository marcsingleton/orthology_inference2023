#!/bin/bash
# Key parameters
#SBATCH --account=fc_eisenlab
#SBATCH --partition=savio2
#SBATCH --time=48:00:00
#SBATCH --qos=savio_normal
#
# Process parameters
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#
# Reporting parameters
#SBATCH --job-name=iqtree_GTR
#SBATCH --output=iqtree_GTR.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
module load gnu-parallel

for folder in $(ls ../make_metaNT/out/); do
  if [[ $folder == *_ni ]]; then
    invariant=""
  else
    invariant="+I"
  fi
  for file in $(ls ../make_metaNT/out/${folder}/); do
    if [[ $file == *.fasta ]]; then
      echo $folder $(basename $file .fasta) $invariant
    fi
  done
done | parallel --jobs $SLURM_CPUS_ON_NODE --colsep ' ' bash iqtree_GTR.sh {1} {2} {3}