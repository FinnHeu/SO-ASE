
# =============================================================================
# ============================== READ MODULES =================================
# =============================================================================
import warnings
warnings.filterwarnings("ignore")

import os
import so_ase as so


# =============================================================================
# ============================== USER SETTINGS ================================
# =============================================================================

# Mesh directories
path_mesh = "/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2cav/"

# Restart Source (already generated restarts on DARS2cav grid for the branchoff initialization run)
restart_year = 1599
path_restart_src_oce = f"/work/ba1550/a270186/simulations/awiesm3-v3.4.1/restarts/DARS2_to_DARS2cav/CAVini/fesom.{restart_year}.oce.restart/"
path_restart_src_ice = f"/work/ba1550/a270186/simulations/awiesm3-v3.4.1/restarts/DARS2_to_DARS2cav/CAVini/fesom.{restart_year}.ice.restart/"

# Cavity Restart (restarts from which to take the cavity values)
path_restart_cavity_fill_oce = "/work/ba1550/a270301/runtime/awiesm3-v3.4.1/branchoff_DARS2cav/restart/fesom/fesom.1611.oce.restart/"
path_restart_cavity_fill_ice = "/work/ba1550/a270301/runtime/awiesm3-v3.4.1/branchoff_DARS2cav/restart/fesom/fesom.1611.ice.restart/"

# Restart Destination (location where restart files are going to be stored)
path_restart_dst_oce = f"/work/ba1550/a270186/simulations/awiesm3-v3.4.1/restarts/DARS2_to_DARS2cav/CAV/fesom.{restart_year}.oce.restart/"
path_restart_dst_ice = f"/work/ba1550/a270186/simulations/awiesm3-v3.4.1/restarts/DARS2_to_DARS2cav/CAV/fesom.{restart_year}.ice.restart/"

# Plots
plot = True
path_dst_plots = "./plots/refill/"

if __name__ == "__main__":
    # Create destination directories if they do not exist
    if not os.path.isdir(path_restart_dst_oce):
        os.makedirs(path_restart_dst_oce)
    if not os.path.isdir(path_restart_dst_ice):
        os.makedirs(path_restart_dst_ice)
    if plot and not os.path.isdir(path_dst_plots):
        os.makedirs(path_dst_plots)

    # Read mesh coordinates for plotting
    node_lon, node_lat, node_idx, _ = so.read_nodes(path_mesh)
    elem, elem_lon, elem_lat = so.read_elements(path_mesh, return_coordinates=True)

    # Variables to process/copy
    vars_oce = ['salt', 'temp', 'temp_AB', 'salt_AB', 'temp_M1', 'salt_M1']
    vars_oce_copy = ['ssh', 'ssh_rhs_old', 'hbar', 'hnode', 'u', 'v', 'vrhs_AB', 'urhs_AB', 'urhs_AB3', 'vrhs_AB3', 'w', 'w_impl', 'w_expl'] 
    vars_ice_copy = ['area', 'hsnow', 'hice', 'uice', 'vice', 'ice_temp', 'ice_albedo']
    

    # Process required ocean variables which need cavity refill
    for var in vars_oce:
        so.fill_cavities_from_existing_restart(var, path_restart_src_oce, path_restart_dst_oce, path_restart_cavity_fill_oce, path_mesh)
        if plot:
            so.plot_refill_comparison(var, path_restart_src_oce, path_restart_dst_oce, path_restart_cavity_fill_oce,
                                   node_lon, node_lat, path_dst_plots)

    # Copy variables that don't need cavity filling
    for var in vars_oce_copy:
        print(f'Copying {var}...')
        os.system(f'cp {path_restart_src_oce}/{var}.nc {path_restart_dst_oce}/{var}.nc')
    
    for var in vars_ice_copy:
        print(f'Copying {var}...')
        os.system(f'cp {path_restart_src_ice}/{var}.nc {path_restart_dst_ice}/{var}.nc')

    print('=============================================================================')
    print('============================ REFILL COMPLETE ================================')
    print('=============================================================================')  