#!/bin/bash
#SBATCH --job-name=interpolate_fesom
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --mem=10G
#SBATCH --exclusive
#SBATCH --time=00:15:00
#SBATCH --account=ab0995

# Set library path so pyproj finds compatible libstdc++
#export LD_LIBRARY_PATH=$HOME/.conda/envs/pyfesom2/lib:$LD_LIBRARY_PATH

############ Define Parameters ############
conda_env=/work/ab0995/a270301/myconda/envs/so_ase/
base_path_model="/work/ba1550/a270301/runtime/awiesm3-v3.4.1/branchoff_DARS2cav/outdata/fesom/"
base_path_mesh="/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2cav/" 
dest_path="/work/ab0995/a270301/SO-ASE/workflows_fesom/workflows_output/workflows_postprocessing/"
variables_to_interpolate=("a_ice" "temp")  #"sst" "sss" "MLD1" "MLD2" "MLD3" "thdgrice" "dyngrice" "fw_ice"  
min_lat=-90
max_lat=-40
min_lon=-180
max_lon=180
lat_increment=0.5
lon_increment=0.5
radius_of_influence=100000 # meters
year_start=1600
year_end=1601

LOG_DIR="./logs/interpolate/"
mkdir -p "$LOG_DIR"  # Ensure log directory exists

# Process all years in parallel using srun  
MAX_PARALLEL=128 # Make sure this is not higher than the number of tasks since the job will get stuck otherwise
i=0

for ((year=year_start; year <= year_end; year++)); do
    echo "$year"
    srun --exclusive --ntasks=1 --nodes=1 \
        $conda_env/bin/python /work/ab0995/a270301/SO-ASE/workflows_fesom/workflows_output/workflows_postprocessing/python/interpolate_fesom_to_regular.py \
        "$year" "$base_path_model" "$base_path_mesh" \
        "$min_lat" "$max_lat" "$min_lon" "$max_lon" \
        "$lat_increment" "$lon_increment" "$radius_of_influence" \
        "$dest_path" "${variables_to_interpolate[@]}" \
        > "${LOG_DIR}/log_${year}.log" 2>&1 &

    ((i++))
    if (( i % MAX_PARALLEL == 0 )); then
        wait
    fi
done
wait

echo "End of script"
