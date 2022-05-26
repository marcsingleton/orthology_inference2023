#!/bin/bash
# Key parameters
#SBATCH --account=fc_eisenlab
#SBATCH --partition=savio
#SBATCH --time=24:00:00
#SBATCH --qos=savio_normal
#
# Process parameters
#SBATCH --nodes=1
#
# Reporting parameters
#SBATCH --job-name=decode
#SBATCH --output=decode.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
# Link to output in scratch
if [ ! -d out ]; then
  out_dir=/global/scratch/users/singleton/IDREvoDevo/analysis/ortho_MSA/deletion_decode/out/
  if [ ! -d ${out_dir} ]; then
    mkdir -p ${out_dir}  # -p makes intermediate directory if they do not exist
  fi
  ln -s ${out_dir} out
fi

source /global/home/users/singleton/.bashrc
conda activate IDREvoDevo
python decode.py
