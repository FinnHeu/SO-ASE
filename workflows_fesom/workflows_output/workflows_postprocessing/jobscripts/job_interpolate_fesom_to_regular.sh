#!/bin/bash
#SBATCH --job-name=interpolate_fesom
#SBATCH -p compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH -o slurm-out.out
#SBATCH -e slurm-err.out
#SBATCH -A ab0995
#SBATCH --array=1650-1849%128

# Set library path so pyproj finds compatible libstdc++
#export LD_LIBRARY_PATH=$HOME/.conda/envs/pyfesom2/lib:$LD_LIBRARY_PATH

############ Define Parameters ############
base_path_model="/work/ab0995/a270186/esm_tools/runtime/awicm3-develop/production_OSM26/PICTRL_CORE2_1650_1850_2/outdata/fesom/"
base_path_mesh="/work/ab0995/a270186/model_inputs/fesom2/mesh/CORE2/" 
variables_to_interpolate=("a_ice" "sst" "sss" "MLD1")  
min_lat=-90
max_lat=90
min_lon=-180
max_lon=180
lat_increment=0.25
lon_increment=0.25
radius_of_influence=100000 # meters

LOG_DIR="./logs/interpolate/"
mkdir -p "$LOG_DIR"  # Ensure log directory exists

# Get the current year from SLURM array task ID
year=$SLURM_ARRAY_TASK_ID

echo "Launching interpolation for year $year"
~/.conda/envs/so_ase/bin/python ./../python/interpolate_fesom_to_regular.py \
    "$year" "$base_path_model" "$base_path_mesh" \
    "$min_lat" "$max_lat" "$min_lon" "$max_lon" \
    "$lat_increment" "$lon_increment" "$radius_of_influence" \
    "${variables_to_interpolate[@]}" \
    > "$LOG_DIR/interpolate_$year.log" 2>&1

echo "Interpolation job for year $year completed."
