#!/bin/bash
#SBATCH --account=oze.oze
#SBATCH --job-name=unrotate_parallel
#SBATCH --partition=mpp
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36

# Load your environment
#micromamba /albedo/home/fheukamp/miniconda3/envs/BarentsSea

# Output directory for logs
LOG_DIR=logs_unrotate
mkdir -p "$LOG_DIR"

############ Define Parameters ############
base_path_model="/albedo/work/user/ufoezk001/result/1stcycle/" # Directory containing unrotated model velocity output (u.fesom.*.nc, v.fesom.*.nc), fesom.mesh.diag.nc must be in this path
year_start=1979 # Starting year for unrotation
year_end=2014 # Ending year for unrotation

# Loop over years and run in parallel with logging
for year in $(seq $year_start $year_end); do
    echo "Launching unrotation for year $year"
    ~/micromamba/envs/pyfesom2_dev/bin/python unrotate_velocities_on_nodes.py "$year" "$base_path_model" > "$LOG_DIR/unrotate_$year.log" 2>&1 &
done

# Wait for all background jobs to complete
wait

echo "All unrotation jobs completed."
