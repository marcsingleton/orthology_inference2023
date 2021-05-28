#!/bin/bash
# Key parameters
#SBATCH --account=fc_eisenlab
#SBATCH --partition=savio
#SBATCH --time=48:00:00
#SBATCH --qos=savio_normal
#
# Process parameters
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#
# Reporting parameters
#SBATCH --job-name=segment_hmm3_calc
#SBATCH --output=segment_hmm3_calc.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
source /global/home/users/singleton/.bashrc
conda activate IDREvoDevo
python segment_hmm3_calc.py
