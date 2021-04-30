#!/bin/bash
# Key parameters
#SBATCH --account=fc_eisenlab
#SBATCH --partition=savio2
#SBATCH --time=24:00:00
#SBATCH --qos=savio_normal
#
# Process parameters
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#
# Reporting parameters
#SBATCH --job-name=mixture_rate_calc
#SBATCH --output=mixture_rate.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=marcsingleton@berkeley.edu
#
# Command(s) to run:
module load python
python mixture_rate_calc.py
