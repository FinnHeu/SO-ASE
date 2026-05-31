""" 
# =============================================================================
# DARS2/CORE2 to DARS2cav/CORE2ice: Generate Restarts
# =============================================================================
# Author: Finn Heukamp, 2025
#
# Here, FESOM2 restarts are taken from a DARS2/CORE2 simulation and extrapolated
# onto an existing, however empty, DARS2cav/CORE2ice restart.
#
# The aim is to avoid a full cold start with the DARS2cav mesh and instead
# branch off a DARS2 simulation at the end of the restart.
#
# Extrapolation strategy:
#
# Cavity:
#   - Temperature and Salinity: nearest neighbor
#   - Velocity: set to 0 (cold start) in all cavity elements
#   - If desired, cavities are filled from existing restart files to reduce 
#     spin-up time and avoid a freshwater shock
#
# Open Ocean:
#   - Temperature, Salinity, Velocity: nearest neighbor
#
# List of all restart variables
# ----------------------------------------------------------------------------- 
#
# Ocean:
#
# 2D node: (time, node)
#   - ssh.nc          (time, node) 
#   - ssh_rhs_old.nc  (time, node)
#   - hbar.nc         (time, node)
#
# 3D node nz1: (time, nz_1, node)
#   - hnode.nc    (time, nz_1, node)
#   - salt.nc     (time, nz_1, node)
#   - temp.nc     (time, nz_1, node)
#   - temp_AB.nc  (time, nz_1, node)
#   - salt_AB.nc  (time, nz_1, node)
#   - temp_AB3.nc (time, nz_1, node)
#   - salt_AB3.nc (time, nz_1, node)
#   - temp_M1.nc  (time, nz_1, node)
#   - salt_M1.nc  (time, nz_1, node)
#
# 3D node nz: (time, nz, node)
#   - w_impl.nc  (time, nz, node)
#   - w_expl.nc  (time, nz, node)
#   - w.nc       (time, nz_1, node)
#
# 3D element: (time, nz_1, elem)
#   - u.nc        (time, nz_1, elem)
#   - v.nc        (time, nz_1, elem)
#   - vrhs_AB.nc  (time, nz_1, elem)
#   - urhs_AB.nc  (time, nz_1, elem)
#
# Ice:
#   - area.nc
#   - hice.nc
#   - hsnow.nc
#   - uice.nc
#   - vice.nc
#
# =============================================================================
"""

# =============================================================================
# ============================== READ MODULES =================================
# =============================================================================
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import so_ase as so
import os
import sys



# =============================================================================
# ============================== USER SETTINGS ================================
# =============================================================================

# Mesh directories
path_mesh_src = "/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2/"
path_mesh_tgt = "/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2cav/"

# Restart files on DARS2 mesh which are supposed to be used on DARS2cav mesh
path_restart_src_oce = f"/work/bb1469/a270089/runtime/awiesm3-v3.4.1/AWI-ESM3-VEG-HR-CMIP7-Spinup_cont2/restart/fesom/fesom.1599.oce.restart/"
path_restart_src_ice = f"/work/bb1469/a270089/runtime/awiesm3-v3.4.1/AWI-ESM3-VEG-HR-CMIP7-Spinup_cont2/restart/fesom/fesom.1599.ice.restart/"

# Restart files on DARS2cav mesh (template files of which the file structure is taken)
path_restart_tgt_oce = f"/work/ab0995/a270186/model_inputs/awicm3/pool/restarts/templates/DARS2cav/v2.7.1/fesom.oce.restart/"
path_restart_tgt_ice = f"/work/ab0995/a270186/model_inputs/awicm3/pool/restarts/templates/DARS2cav/v2.7.1/fesom.ice.restart/"

# Restart Destination (destination of the generated restart files to be used for the run)
restart_year = 1599
path_restart_dst_oce = f"/work/ba1550/a270186/simulations/awiesm3-v3.4.1/restarts/DARS2_to_DARS2cav/CAVini/fesom.{restart_year}.oce.restart/"
path_restart_dst_ice = f"/work/ba1550/a270186/simulations/awiesm3-v3.4.1/restarts/DARS2_to_DARS2cav/CAVini/fesom.{restart_year}.ice.restart/"

# Plots Destination
plot = True
path_dst_plots = "./plots/"

# Coupled Model (AWI-CM3 with FESOM2.7) also requires ice_temp.nc and ice_albedo.nc
is_coupled = True

# Path to restart files on target mesh ( ---> fesom v2.7 <---  ) for masking cavities
path_restart_tgt_oce_v27 = path_restart_tgt_oce

# =============================================================================
# ============================ SET LOG FILES ==================================
# =============================================================================

# Create destination directories if they do not exist
if not os.path.isdir(path_restart_dst_oce):
    os.makedirs(path_restart_dst_oce)
if not os.path.isdir(path_restart_dst_ice):
    os.makedirs(path_restart_dst_ice)

# Determine the parent folder of path_restart_dst_oce
log_file = os.path.join(os.path.dirname(path_restart_dst_oce.rstrip('/')), 'extrapolate_cavity_restart.log')

# Tee class to write to both terminal and log file
class Tee:
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()

# Redirect stdout and stderr to both terminal and log file
log_file_handle = open(log_file, 'w')
sys.stdout = Tee(sys.__stdout__, log_file_handle)
sys.stderr = Tee(sys.__stderr__, log_file_handle)

if __name__ == "__main__":
    print('=============================================================================')
    print('======================== EXTRAPOLATE CAVITY RESTARTS ========================')
    print('=============================================================================')
    print(' ')

    # =============================================================================
    # ============================ READ MESHES =================================
    # =============================================================================

    # DARS2/CORE2
    node_lon_src, node_lat_src, node_id_src, node_coastmask_src = so.read_nodes(path_mesh_src)
    elem_src = so.read_elements(path_mesh_src)
    elem_lon_src = node_lon_src[elem_src].mean(axis=1)
    elem_lat_src = node_lat_src[elem_src].mean(axis=1)

    # DARS2cav/CORE2ice
    node_lon_tgt, node_lat_tgt, node_id_tgt, node_coastmask_tgt = so.read_nodes(path_mesh_tgt)
    elem_tgt = so.read_elements(path_mesh_tgt)
    elem_lon_tgt = node_lon_tgt[elem_tgt].mean(axis=1)
    elem_lat_tgt = node_lat_tgt[elem_tgt].mean(axis=1)


    # =============================================================================
    # ============================ BUILD MAPPERS ==================================
    # =============================================================================
    print("###---> Building Nearest Neighbor Mappers...")
    print(" ")
    mapper_elements, distance_elements = so.build_spherical_nn_mapper(elem_lon_src, elem_lat_src, elem_lon_tgt, elem_lat_tgt)
    mapper_nodes, distance_nodes = so.build_spherical_nn_mapper(node_lon_src, node_lat_src, node_lon_tgt, node_lat_tgt)

    if plot:
        print("###---> Plotting Nearest Neighbor Mappers...")
        print(" ")
        so.plot_mapper(mapper_elements, elem_lon_src, elem_lat_src, elem_lon_tgt, elem_lat_tgt, 'elem', path_dst_plots)
        so.plot_mapper(mapper_nodes, node_lon_src, node_lat_src, node_lon_tgt, node_lat_tgt, 'node', path_dst_plots)

    # =============================================================================
    # ======================== 2D nodal fiels (time, node) ========================
    # =============================================================================
    varnames_2D_node = ['ssh', 'ssh_rhs_old', 'hbar']
    mask_file = f"{path_restart_tgt_oce_v27}ssh.nc"

    for varname in varnames_2D_node:
        so.interpolate_extrapolate_2D(varname, path_restart_src_oce, path_restart_tgt_oce, mapper_nodes, path_restart_dst_oce, mask_file, mask_file_varname='ssh', t_step=-1, verbose=True)
        if plot:
            so.plot_interpolated_extrapolated_field(path_restart_src_oce, path_restart_dst_oce, path_restart_tgt_oce, varname, node_lon_src, node_lat_src, node_lon_tgt, node_lat_tgt, path_dst_plots)

    varnames_2D_node = ['area', 'hsnow', 'hice', 'uice', 'vice']
    if is_coupled:
        varnames_2D_node.extend(['ice_temp', 'ice_albedo'])

    for varname in varnames_2D_node:
        so.interpolate_extrapolate_2D(varname, path_restart_src_ice, path_restart_tgt_ice, mapper_nodes, path_restart_dst_ice, mask_file, mask_file_varname='ssh', t_step=-1, verbose=True)
        if plot:
            so.plot_interpolated_extrapolated_field(path_restart_src_ice, path_restart_dst_ice, path_restart_tgt_ice, varname, node_lon_src, node_lat_src, node_lon_tgt, node_lat_tgt, path_dst_plots)

    # =============================================================================
    # ======================== 3D nodal/elem fiels (time, nz_1/nz, node/elem) =====
    # =============================================================================
    varnames_3D_node = ['hnode', 'salt', 'temp', 'temp_AB', 'salt_AB', 'temp_M1', 'salt_M1', 'w', 'w_impl', 'w_expl']

    for varname in varnames_3D_node:
        so.interpolate_extrapolate_3D(varname, path_restart_src_oce, path_restart_tgt_oce, mapper_nodes, path_restart_dst_oce, t_step=-1, verbose=True)
        if plot:
            so.plot_interpolated_extrapolated_field(path_restart_src_oce, path_restart_dst_oce, path_restart_tgt_oce, varname, node_lon_src, node_lat_src, node_lon_tgt, node_lat_tgt, path_dst_plots)


    varnames_3D_element = ['u', 'v', 'vrhs_AB', 'urhs_AB', 'urhs_AB3', 'vrhs_AB3']
    for varname in varnames_3D_element:
        so.interpolate_extrapolate_3D(varname, path_restart_src_oce, path_restart_tgt_oce, mapper_elements, path_restart_dst_oce, t_step=-1, verbose=True)
        if plot:
            so.plot_interpolated_extrapolated_field(path_restart_src_oce, path_restart_dst_oce, path_restart_tgt_oce, varname, elem_lon_src, elem_lat_src, elem_lon_tgt, elem_lat_tgt, path_dst_plots)

    
    print('=============================================================================')
    print('======================== EXTRAPOLATION COMPLETE ============================')
    print('=============================================================================')