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
#SBATCH --job-name=iqtree_LG
#SBATCH --output=iqtree_LG.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
module load gnu-parallel

for folder in $(ls ../make_meta_AA/out/); do
  if [[ ! -d ../make_meta_AA/out/${folder}/ ]]; then
    continue
  fi
  if [[ $folder == *_NI ]]; then
    invariant=""
  else
    invariant="+I"
  fi
  for file in $(ls ../make_meta_AA/out/${folder}/); do
    if [[ $file == *.afa ]]; then
      echo $folder $(basename $file .afa) $invariant
    fi
  done
done | parallel --jobs $SLURM_CPUS_ON_NODE --colsep ' ' bash iqtree_LG.sh {1} {2} {3}
