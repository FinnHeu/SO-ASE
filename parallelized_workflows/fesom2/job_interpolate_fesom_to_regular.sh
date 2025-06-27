#!/bin/bash
#SBATCH --account=oze.oze
#SBATCH --job-name=interpolate_parallel
#SBATCH --partition=mpp
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36

# Output directory for logs
LOG_DIR=logs_interpolate
mkdir -p "$LOG_DIR"

############ Define Parameters ############
base_path_model="/albedo/work/user/ufoezk001/result/1stcycle/" # Directory containing model output
base_path_mesh="~/mesh/Arc01/" # Directory containing FESOM2 mesh files
dest_path="/albedo/work/user/ufoezk001/result/1stcycle/interpolated/" # Destination directory for interpolated data
year_start=1979 # Starting year for unrotation
year_end=2014 # Ending year for unrotation
variables_to_interpolate="unod" "vnod" "temp" "salt" # Variables on nodes to interpolate

# New parameters for latitude and longitude bounds and increments
min_lat=-90
max_lat=90
min_lon=-180
max_lon=180
lat_increment=1.0
lon_increment=1.0

# Loop over years and run in parallel with logging
for year in $(seq $year_start $year_end); do
    echo "Launching unrotation for year $year"
    ~/micromamba/envs/pyfesom2_dev/bin/python interpolate_fesom_to_regular.py "$year" "$base_path_model" "$base_path_mesh" "$min_lat" "$max_lat" "$min_lon" "$max_lon" "$lat_increment" "$lon_increment" "$variables_to_interpolate" > "$LOG_DIR/interpolate_$year.log" 2>&1 &
done

# Wait for all background jobs to complete
wait

echo "All interpolation jobs completed."
