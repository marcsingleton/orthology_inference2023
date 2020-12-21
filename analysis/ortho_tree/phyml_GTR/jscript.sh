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
#SBATCH --job-name=phyml_GTR
#SBATCH --output=phyml_GTR.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
module load gcc
module load gnu-parallel

for folder in $(ls ../make_meta/out/); do
  if [[ $folder == *_ni ]]; then
    pinv=0
  else
    pinv=e
  fi
  for file in $(ls ../make_meta/out/${folder}/); do
    if [[ $file == *.phy ]]; then
      echo $folder $(basename $file .phy) $pinv
    fi
  done
done | parallel --jobs $SLURM_CPUS_ON_NODE --colsep ' ' bash phyml_GTR.sh {1} {2} {3}
