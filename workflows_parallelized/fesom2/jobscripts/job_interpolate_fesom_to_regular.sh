#!/bin/bash
#SBATCH --account=nwg_so-ase.nwg_so-ase
#SBATCH --job-name=interpolate_parallel
#SBATCH --partition=mpp
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5

# Set library path so pyproj finds compatible libstdc++
export LD_LIBRARY_PATH=$HOME/.conda/envs/pyfesom2/lib:$LD_LIBRARY_PATH

# Output directory for logs
LOG_DIR=logs_interpolate
mkdir -p "$LOG_DIR"

############ Define Parameters ############
base_path_model="/albedo/work/user/fheukamp/PostDoc1/FramStraitWSC/results_CTRL_cvmix/"
base_path_mesh="$HOME/mesh/Arc01/" 
year_start=2000
year_end=2005
variables_to_interpolate=("unod" "vnod" "temp" "salt")  
min_lat=-90
max_lat=90
min_lon=-180
max_lon=180
lat_increment=1.0
lon_increment=1.0

LOG_DIR="./logs/interpolate/"
mkdir -p "$LOG_DIR"  # Ensure log directory exists

# Loop over years
for year in $(seq $year_start $year_end); do
    echo "Launching interpolation for year $year"
    ~/.conda/envs/pyfesom2/bin/python ./../python/interpolate_fesom_to_regular.py \
        "$year" "$base_path_model" "$base_path_mesh" \
        "$min_lat" "$max_lat" "$min_lon" "$max_lon" \
        "$lat_increment" "$lon_increment" \
        "${variables_to_interpolate[@]}" \
        > "$LOG_DIR/interpolate_$year.log" 2>&1 &
done

wait  # Wait for all background jobs to finish

# Wait for all background jobs to complete
wait

echo "All interpolation jobs completed."
