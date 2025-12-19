#!/bin/bash
#SBATCH --account=nwg_so-ase.nwg_so-ase
#SBATCH --job-name=smooth_mesh_2
#SBATCH --partition=mpp
#SBATCH --mem=40GB
#SBATCH --time=12:00:00
#SBATCH --qos=12h
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

start_time=$SECONDS

module load netcdf-fortran/4.5.4-oneapi2022.1.0 
module load intel-oneapi-compilers/2022.1.0 
module load intel-oneapi-mpi/2021.6.0 

./smooth_depth_z_new.x > smooth_depth_cav_v4.log 2>&1

