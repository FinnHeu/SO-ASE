#!/bin/bash
#SBATCH --account=nwg_so-ase.nwg_so-ase
#SBATCH --job-name=distribute_runoff
#SBATCH --partition=mpp
#SBATCH --mem=64GB
#SBATCH --time=12:00:00
#SBATCH --qos=12h
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

module load matlab

matlab -nodisplay -nosplash -r "run('./../matlab/distribute_runoff.m'); exit;" > river2coast.log 2>&1