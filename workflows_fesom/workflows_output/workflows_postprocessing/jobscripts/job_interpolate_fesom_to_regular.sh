#!/bin/bash
#SBATCH --job-name=interpolate_fesom
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00
#SBATCH -o slurm-out.out
#SBATCH -e slurm-err.out
#SBATCH -A ab0995

# Set library path so pyproj finds compatible libstdc++
#export LD_LIBRARY_PATH=$HOME/.conda/envs/pyfesom2/lib:$LD_LIBRARY_PATH

# Output directory for logs
LOG_DIR=logs_interpolate
mkdir -p "$LOG_DIR"

############ Define Parameters ############
base_path_model="/work/ab0995/a270186/esm_tools/runtime/awicm3-develop/CORE2ice_TEST/outdata/fesom/"
base_path_mesh="/work/ab0995/a270186/model_inputs/fesom2/mesh/CORE2ice/" 
year_start=1859
year_end=1859
variables_to_interpolate=("unod" "vnod" "temp" "salt" "fh" "fw" "sst")  
min_lat=-90
max_lat=-50
min_lon=-180
max_lon=180
lat_increment=0.25
lon_increment=0.25

LOG_DIR="./logs/interpolate/"
mkdir -p "$LOG_DIR"  # Ensure log directory exists

# Loop over years
for year in $(seq $year_start $year_end); do
    echo "Launching interpolation for year $year"
    ~/.conda/envs/pyfesom2_dev/bin/python ./../python/interpolate_fesom_to_regular.py \
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
