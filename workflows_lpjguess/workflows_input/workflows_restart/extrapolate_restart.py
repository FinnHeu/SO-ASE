#!/usr/bin/env python3
"""
# =============================================================================
# DARS2 to DARS2cav: LPJ-GUESS Restart File Extrapolation
# =============================================================================
# Author: Finn Heukamp, 2025
#
# This script converts LPJ-GUESS restart files from a DARS2 simulation to 
# be compatible with a DARS2cav simulation using a template-first approach.
#
# Conversion strategy:
#   1. For cells that exist in both DARS2 and target: use DARS2 data
#   2. For new cells: use template data from the DARS2cav template files
#   3. Only use nearest neighbor as last resort
#
# This approach avoids propagating invalid soil data from DARS2 restarts
# and ensures compatibility with the cavity-enabled DARS2cav mesh.
#
# The script reads LPJ-GUESS state files (.state) which contain binary
# vegetation and soil data for each grid cell, along with spatial indexing
# information for efficient cell lookup.
#
# =============================================================================
"""

# =============================================================================
# ============================== READ MODULES =================================
# =============================================================================
import sys
import logging

# Import LPJ-GUESS restart helper functions
import so_ase.lpjguess.helpers_restart as lpj_helpers


# =============================================================================
# ============================== USER SETTINGS ================================
# =============================================================================

# Source restart files (DARS2 simulation)
source_dir = '/work/bb1469/a270089/runtime/awiesm3-v3.4.2/AWI-ESM3-VEG-HR-CMIP7-Spinup_cont3/restart/lpj_guess/lpjg_state_1680'

# Template restart files (DARS2cav/CAV-ICB structure)
template_dir = '/work/ba1550/a270186/simulations/awiesm3-v3.4.1-CAV-ICB/preproduction/CAV-ICB-PICTRL-TEST-with_icb-7/restart/lpj_guess/lpjg_state_1608'

# Output directory for converted restart files
output_dir = '/home/a/a270186/python_modules/SO-ASE/test/'

# =============================================================================
# ============================ SET LOG FILES ==================================
# =============================================================================

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


if __name__ == '__main__':
    print('=============================================================================')
    print('====================== LPJ-GUESS RESTART EXTRAPOLATION =====================')
    print('=============================================================================')
    print(' ')
    
    try:
        # Build indices for source and template
        source_index, source_coords = lpj_helpers.build_cell_index(source_dir, "source")
        template_index, template_coords = lpj_helpers.build_cell_index(template_dir, "template")
        
        # Convert restart files using template-first approach
        stats, nn_distances = lpj_helpers.convert_restart_files_template_first(
            source_dir, template_dir, output_dir
        )
        
        # Verify conversion results
        lpj_helpers.verify_conversion(output_dir, template_dir)
        
        # Print detailed summary
        lpj_helpers.print_conversion_summary(
            stats, nn_distances, output_dir, source_index, template_index
        )
        
        logger.info("")
        logger.info("Conversion successful!")
        sys.exit(0)
        
    except Exception as e:
        logger.error(f"Conversion failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
