"""
FESOM2 Postprocessing Pipeline

This script combines all preprocessing/postprocessing functions from so_ase.fesom2
into a single configurable pipeline. Each processing step can be enabled/disabled
via switches at the top of the configuration section.

Author: Auto-generated
Updated: July 2026

Usage:
    python fesom_postprocessing.py <YEAR>

Arguments:
    <YEAR>  Year to process (passed from SLURM array job)

Dependencies:
    - so_ase
    - numpy
    - xarray
"""

import sys
import os
import numpy as np

# Add parent directory to path if needed
sys.path.insert(0, os.path.expanduser("~/python_modules/SO-ASE"))

import so_ase as so

# ============================================================================
#                         PROCESSING SWITCHES
# ============================================================================
# Set to True to enable, False to disable each processing step

# Ocean processing
DO_EXTRACT_BOTTOM_TEMP = False
DO_EXTRACT_BOTTOM_SALT = False
DO_EXTRACT_BOTTOM_SIGMA0 = False
DO_DEPTH_AVG_TEMP_0_50 = False
DO_DEPTH_AVG_TEMP_50_500 = False
DO_DEPTH_AVG_SALT_0_50 = False
DO_DEPTH_AVG_SALT_50_500 = False
DO_TOTAL_KINETIC_ENERGY = False
DO_DEPTH_EXTRACT = [
    {'varname': 'temp', 'depth': 5},
    {'varname': 'temp', 'depth': 50},
    {'varname': 'temp', 'depth': 250},
    {'varname': 'temp', 'depth': 500},
    {'varname': 'temp', 'depth': 1000},
    {'varname': 'temp', 'depth': 4000},
    {'varname': 'salt', 'depth': 5},
    {'varname': 'salt', 'depth': 50},
    {'varname': 'salt', 'depth': 250},
    {'varname': 'salt', 'depth': 500},
    {'varname': 'salt', 'depth': 1000},
    {'varname': 'salt', 'depth': 4000},
]
# Cavity/ice shelf processing
DO_SUBSHELF_FRESHWATERFLUX = False
DO_SUBSHELF_HEATFLUX = False
DO_FRESHWATERFLUX_TO_GTY = False
DO_FRESHWATERFLUX_TO_GTM = False
DO_SUBSHELF_HYDROGRAPHY_TEMP = False
DO_SUBSHELF_HYDROGRAPHY_SALT = False

# Sea ice processing
DO_SEA_ICE_AREA = True
DO_SEA_ICE_VOLUME = False
DO_SEA_ICE_EXTENT = False

# Region mean vertical profiles
DO_REGION_MEAN_PROFILE = False
# DO_REGION_MEAN_PROFILE = [
#     {'varname': 'temp'}, 
#     {'varname': 'salt'}
#     ]  # Set to False to disable

# ============================================================================
#                         PATH CONFIGURATION
# ============================================================================

# Input paths
SRC_PATH_FESOM = "/work/bb1469/a270089/runtime/awiesm3-v3.4.2/AWI-ESM3-VEG-HR-CMIP7-abrupt4xCO2_1949/outdata/fesom/"  # FESOM output files
SRC_PATH_FESOM_PREPROCESSED = "/work/ba1550/a270186/simulations/awiesm3-v3.4.2-CAV-ICB/production/CMIP7-4CO2/outdata/fesom/"  # FESOM output files
MESH_PATH = "/work/ab0995/a270186/model_inputs/fesom2/mesh/DARS2/"  # Mesh files
MESH_DIAG_PATH = MESH_PATH  # fesom.mesh.diag.nc

# Output paths
DEST_PATH_BASE = "/work/ba1550/a270186/simulations/awiesm3-v3.4.2-CAV-ICB/production/CMIP7-4CO2/analysis/fesom/"
DEST_PATH_BOTTOM = DEST_PATH_BASE + "bottom_fields/"
DEST_PATH_KE = DEST_PATH_BASE + "kinetic_energy/"
DEST_PATH_SUBSHELF = DEST_PATH_BASE + "subshelf/"
DEST_PATH_SUBSHELF_GTY = DEST_PATH_BASE + "subshelf_gty/"
DEST_PATH_SUBSHELF_GTM = DEST_PATH_BASE + "subshelf_gtm/"
DEST_PATH_SUBSHELF_HYDRO = DEST_PATH_BASE + "subshelf_hydrography/"
DEST_PATH_SEA_ICE = DEST_PATH_BASE + "sea_ice/"
DEST_PATH_DEPTH_AVG = DEST_PATH_BASE + "depth_averages/"
DEST_PATH_DEPTH_EXTRACT = DEST_PATH_BASE + "depth_extracts/"
DEST_PATH_REGION_AVG = DEST_PATH_BASE + "region_averages/"

# Mask files for region averaging (list of mask file paths)
MASK_FILES_REGION = [
    MESH_PATH + "mask.AAshelf.max1000m.nc",
]

# KML path for regional masks (cavity processing)
KML_PATH = os.path.expanduser("~/python_modules/SO-ASE/files/kml/")

# ============================================================================
#                         PROCESSING PARAMETERS
# ============================================================================

# Sea ice parameters
SEA_ICE_BOX_SO = [-180, 180, -90, -45]  # [lon_min, lon_max, lat_min, lat_max]
SEA_ICE_BOX_ARCTIC = [-180, 180, 45, 90]  # [lon_min, lon_max, lat_min, lat_max]
SEA_ICE_BOX_ARCMED = [-40, 100, 40, 80]  # [lon_min, lon_max, lat_min, lat_max]
SICONC_THRESHOLD = 0.15

# Cavity mask configuration
# List of all cavity masks to process (set to empty list [] to skip regional processing)
# Note: kml_path should be the directory, not the full file path (filename is built from name)
CAVITY_MASKS = [
    {'name': 'all'},  # All cavities combined
    {'name': 'Amery', 'kml_path': KML_PATH},
    {'name': 'Amundsen', 'kml_path': KML_PATH},
    {'name': 'AustralianSector', 'kml_path': KML_PATH},
    {'name': 'Bellingshausen', 'kml_path': KML_PATH},
    {'name': 'EasternWeddell', 'kml_path': KML_PATH},
    {'name': 'FilchnerRonne', 'kml_path': KML_PATH},
    {'name': 'Larsen', 'kml_path': KML_PATH},
    {'name': 'Ross', 'kml_path': KML_PATH},
]

# Kinetic energy mask ('all', 'cavity', 'open_ocean')
KE_MASK = 'cavity'

# Subshelf hydrography variables
SUBSHELF_HYDRO_VARIABLES = ['temp', 'salt']

# Logging
LOG = True

# ============================================================================
#                         MAIN PROCESSING FUNCTION
# ============================================================================

def run_postprocessing(year):
    """Run all enabled postprocessing steps for a single year."""
    
    print("=" * 78)
    print(f"{'FESOM2 Postprocessing Pipeline':^78}")
    print("=" * 78)
    print(f"Processing year: {year}")
    print(f"Source path: {SRC_PATH_FESOM}")
    print(f"Mesh path: {MESH_PATH}")
    print("=" * 78)
    
    # Create output directories
    for path in [DEST_PATH_BOTTOM, DEST_PATH_KE, DEST_PATH_SUBSHELF,
                 DEST_PATH_SUBSHELF_GTY, DEST_PATH_SUBSHELF_GTM,
                 DEST_PATH_SUBSHELF_HYDRO, DEST_PATH_SEA_ICE, DEST_PATH_DEPTH_AVG,
                 DEST_PATH_DEPTH_EXTRACT, DEST_PATH_REGION_AVG]:
        os.makedirs(path, exist_ok=True)
    
    years = (year, year)  # Process single year
    
    # -------------------------------------------------------------------------
    # Ocean processing
    # -------------------------------------------------------------------------
    
    if DO_EXTRACT_BOTTOM_TEMP:
        print("\n>>> Extracting bottom temperature...")
        so.fesom_extract_bottom_values(
            src_path=SRC_PATH_FESOM,
            dest_path=DEST_PATH_BOTTOM,
            varname='temp',
            years=years,
            log=LOG
        )
    
    if DO_EXTRACT_BOTTOM_SALT:
        print("\n>>> Extracting bottom salinity...")
        so.fesom_extract_bottom_values(
            src_path=SRC_PATH_FESOM,
            dest_path=DEST_PATH_BOTTOM,
            varname='salt',
            years=years,
            log=LOG
        )
    
    if DO_EXTRACT_BOTTOM_SIGMA0:
        print("\n>>> Extracting bottom density (sigma0)...")
        so.fesom_extract_bottom_values(
            src_path=SRC_PATH_FESOM,
            dest_path=DEST_PATH_BOTTOM,
            varname='sigma0',
            years=years,
            log=LOG
        )
    
    if DO_DEPTH_AVG_TEMP_0_50:
        print("\n>>> Computing depth-averaged temperature (0-50m)...")
        so.fesom_depth_average(
            src_path=SRC_PATH_FESOM,
            dest_path=DEST_PATH_DEPTH_AVG,
            mesh_diag_path=MESH_DIAG_PATH,
            varname='temp',
            depth_range=(0, 50),
            years=years,
            log=LOG
        )
    
    if DO_DEPTH_AVG_TEMP_50_500:
        print("\n>>> Computing depth-averaged temperature (50-500m)...")
        so.fesom_depth_average(
            src_path=SRC_PATH_FESOM,
            dest_path=DEST_PATH_DEPTH_AVG,
            mesh_diag_path=MESH_DIAG_PATH,
            varname='temp',
            depth_range=(50, 500),
            years=years,
            log=LOG
        )
    
    if DO_DEPTH_AVG_SALT_0_50:
        print("\n>>> Computing depth-averaged salinity (0-50m)...")
        so.fesom_depth_average(
            src_path=SRC_PATH_FESOM,
            dest_path=DEST_PATH_DEPTH_AVG,
            mesh_diag_path=MESH_DIAG_PATH,
            varname='salt',
            depth_range=(0, 50),
            years=years,
            log=LOG
        )
    
    if DO_DEPTH_AVG_SALT_50_500:
        print("\n>>> Computing depth-averaged salinity (50-500m)...")
        so.fesom_depth_average(
            src_path=SRC_PATH_FESOM,
            dest_path=DEST_PATH_DEPTH_AVG,
            mesh_diag_path=MESH_DIAG_PATH,
            varname='salt',
            depth_range=(50, 500),
            years=years,
            log=LOG
        )
    
    if DO_TOTAL_KINETIC_ENERGY:
        print("\n>>> Computing total kinetic energy...")
        so.fesom_total_kinetic_energy(
            src_path=SRC_PATH_FESOM,
            mesh_diag_path=MESH_DIAG_PATH,
            meshpath=MESH_PATH,
            years=years,
            mask=KE_MASK,
            log=LOG,
            savepath=DEST_PATH_KE
        )
    
    if DO_DEPTH_EXTRACT:
        for extract_config in DO_DEPTH_EXTRACT:
            varname = extract_config['varname']
            depth = extract_config['depth']
            print(f"\n>>> Extracting {varname} at {depth}m...")
            so.fesom_extract_at_depth(
                src_path=SRC_PATH_FESOM,
                dest_path=DEST_PATH_DEPTH_EXTRACT,
                mesh_diag_path=MESH_DIAG_PATH,
                varname=varname,
                depth=depth,
                years=years,
                log=LOG
            )
    
    # -------------------------------------------------------------------------
    # Cavity/ice shelf processing (loop over all regional masks)
    # -------------------------------------------------------------------------
    
    for cavity_mask in CAVITY_MASKS:
        mask_name = cavity_mask['name']
        print(f"\n{'='*40}")
        print(f"Processing cavity region: {mask_name}")
        print(f"{'='*40}")
        
        if DO_SUBSHELF_FRESHWATERFLUX:
            print(f"\n>>> Computing subshelf freshwater flux ({mask_name})...")
            so.fesom_subshelf_freshwaterflux(
                src_path=SRC_PATH_FESOM,
                mesh_diag_path=MESH_DIAG_PATH,
                mesh_path=MESH_PATH,
                mask=cavity_mask,
                years=years,
                log=LOG,
                savepath=DEST_PATH_SUBSHELF
            )
        
        if DO_SUBSHELF_HEATFLUX:
            print(f"\n>>> Computing subshelf heat flux ({mask_name})...")
            so.fesom_subshelf_heatflux(
                src_path=SRC_PATH_FESOM,
                mesh_diag_path=MESH_DIAG_PATH,
                mesh_path=MESH_PATH,
                mask=cavity_mask,
                years=years,
                log=LOG,
                savepath=DEST_PATH_SUBSHELF
            )
        
        if DO_SUBSHELF_HYDROGRAPHY_TEMP:
            print(f"\n>>> Extracting subshelf temperature ({mask_name})...")
            so.fesom_subshelf_hydrography(
                src_path=SRC_PATH_FESOM,
                mesh_diag_path=MESH_DIAG_PATH,
                mesh_path=MESH_PATH,
                mask=cavity_mask,
                years=years,
                variable='temp',
                mean=True,
                log=LOG,
                savepath=DEST_PATH_SUBSHELF_HYDRO
            )
        
        if DO_SUBSHELF_HYDROGRAPHY_SALT:
            print(f"\n>>> Extracting subshelf salinity ({mask_name})...")
            so.fesom_subshelf_hydrography(
                src_path=SRC_PATH_FESOM,
                mesh_diag_path=MESH_DIAG_PATH,
                mesh_path=MESH_PATH,
                mask=cavity_mask,
                years=years,
                variable='salt',
                mean=True,
                log=LOG,
                savepath=DEST_PATH_SUBSHELF_HYDRO
            )
    
    # Unit conversion (only for current year)
    if DO_FRESHWATERFLUX_TO_GTY:
        print("\n>>> Converting freshwater flux to Gt/yr...")
        so.freshwaterflux_to_massflux_Gty(
            src_path=DEST_PATH_SUBSHELF,
            dst_path=DEST_PATH_SUBSHELF_GTY,
            rho_fw=1000,
            year=year,
            log=LOG
        )
    
    if DO_FRESHWATERFLUX_TO_GTM:
        print("\n>>> Converting freshwater flux to Gt/month...")
        so.freshwaterflux_to_massflux_Gtm(
            src_path=DEST_PATH_SUBSHELF,
            dst_path=DEST_PATH_SUBSHELF_GTM,
            rho_fw=1000,
            year=year,
            log=LOG
        )
    
    # -------------------------------------------------------------------------
    # Sea ice processing
    # -------------------------------------------------------------------------
    
    if DO_SEA_ICE_AREA:
        print("\n>>> Computing sea ice area...")
        for box in [SEA_ICE_BOX_SO, SEA_ICE_BOX_ARCTIC, SEA_ICE_BOX_ARCMED]:
            so.fesom_sea_ice_area(
                src_path=SRC_PATH_FESOM_PREPROCESSED,
                mesh_diag_path=MESH_DIAG_PATH,
                years=years,
                box=box,
                siconc_threshold=SICONC_THRESHOLD,
                savepath=DEST_PATH_SEA_ICE,
                log=LOG
            )
    
    if DO_SEA_ICE_VOLUME:
        print("\n>>> Computing sea ice volume...")
        for box in [SEA_ICE_BOX_SO, SEA_ICE_BOX_ARCTIC, SEA_ICE_BOX_ARCMED]:
            so.fesom_sea_ice_volume(
                src_path=SRC_PATH_FESOM_PREPROCESSED,
                mesh_diag_path=MESH_DIAG_PATH,
                years=years,
                box=box,
                savepath=DEST_PATH_SEA_ICE,
                log=LOG
            )
    
    if DO_SEA_ICE_EXTENT:
        print("\n>>> Computing sea ice extent...")
        for box in [SEA_ICE_BOX_SO, SEA_ICE_BOX_ARCTIC, SEA_ICE_BOX_ARCMED]:
            so.fesom_sea_ice_extent(
                src_path=SRC_PATH_FESOM_PREPROCESSED,
                mesh_diag_path=MESH_DIAG_PATH,
                years=years,
                box=box,
                siconc_threshold=SICONC_THRESHOLD,
                savepath=DEST_PATH_SEA_ICE,
                log=LOG
            )
    
    # -------------------------------------------------------------------------
    # Region mean vertical profiles
    # -------------------------------------------------------------------------
    
    if DO_REGION_MEAN_PROFILE:
        for mask_file in MASK_FILES_REGION:
            for profile_config in DO_REGION_MEAN_PROFILE:
                varname = profile_config['varname']
                print(f"\n>>> Computing region mean vertical profile for {varname} ({mask_file})...")
                so.fesom_region_mean_vertical_profile(
                    src_path=SRC_PATH_FESOM,
                    dest_path=DEST_PATH_REGION_AVG,
                    mesh_diag_path=MESH_DIAG_PATH,
                    mask_file=mask_file,
                    varname=varname,
                    years=years,
                    log=LOG
                )
        
    print("\n" + "=" * 78)
    print(f"Postprocessing complete for year {year}")
    print("=" * 78)


# ============================================================================
#                         ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fesom_postprocessing.py <YEAR>")
        print("Example: python fesom_postprocessing.py 1979")
        sys.exit(1)
    
    year = int(sys.argv[1])
    run_postprocessing(year)
