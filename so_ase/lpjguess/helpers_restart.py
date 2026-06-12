# so_ase/lpjguess/helpers_restart.py

import os
import struct
import numpy as np
from glob import glob
from pathlib import Path
from scipy.spatial import cKDTree
from typing import Dict, List, Tuple, Set
import logging


def read_state_index(state_file: str) -> List[Tuple[float, float, int, int]]:
    """
    Read the index from a LPJ-GUESS state file.
    
    Returns list of (lon, lat, file_position, data_size) tuples.
    The index is stored at the end of the file, sorted by coordinates.
    Data is stored in the order cells were serialized (not sorted).
    
    Parameters
    ----------
    state_file : str
        Path to the LPJ-GUESS state file.
        
    Returns
    -------
    List[Tuple[float, float, int, int]]
        List of tuples containing (longitude, latitude, file_position, data_size).
    """
    with open(state_file, 'rb') as f:
        # Get file size
        f.seek(0, 2)
        file_size = f.tell()
        
        # Read number of elements (last 8 bytes)
        f.seek(-8, 2)
        num_elements = struct.unpack('<Q', f.read(8))[0]
        
        if num_elements == 0:
            return []
        
        # Calculate index size (each entry: 16 bytes key + 8 bytes position)
        index_entry_size = 16 + 8  # 2 doubles + 1 streamsize
        index_size = num_elements * index_entry_size
        
        # Seek to start of index
        index_start = file_size - 8 - index_size
        f.seek(index_start)
        
        entries = []
        for _ in range(num_elements):
            lon, lat = struct.unpack('<dd', f.read(16))
            pos = struct.unpack('<q', f.read(8))[0]
            entries.append((lon, lat, pos))
        
        # Sort entries by position to calculate sizes correctly
        # (data is written in order, but index is sorted by coordinates)
        entries_by_pos = sorted(entries, key=lambda x: x[2])
        
        # Calculate data sizes based on position order
        pos_to_size = {}
        for i, (lon, lat, pos) in enumerate(entries_by_pos):
            if i < len(entries_by_pos) - 1:
                next_pos = entries_by_pos[i + 1][2]
            else:
                next_pos = index_start
            pos_to_size[pos] = next_pos - pos
        
        # Return entries with sizes
        result = []
        for lon, lat, pos in entries:
            size = pos_to_size[pos]
            result.append((lon, lat, pos, size))
        
        return result


def read_cell_data(state_file: str, position: int, size: int) -> bytes:
    """
    Read raw cell data from a state file.
    
    Parameters
    ----------
    state_file : str
        Path to the LPJ-GUESS state file.
    position : int
        File position where the cell data starts.
    size : int
        Size of the cell data in bytes.
        
    Returns
    -------
    bytes
        Raw cell data.
    """
    with open(state_file, 'rb') as f:
        f.seek(position)
        return f.read(size)


def build_cell_index(directory: str, name: str) -> Tuple[Dict[Tuple[float, float], Tuple[str, int, int]], np.ndarray]:
    """
    Build an index of all cells in state files.
    
    Parameters
    ----------
    directory : str
        Directory containing the state files.
    name : str
        Name for logging purposes (e.g., "source", "template").
        
    Returns
    -------
    Tuple[Dict[Tuple[float, float], Tuple[str, int, int]], np.ndarray]
        - Dictionary mapping (lon, lat) to (filename, position, size)
        - Numpy array of all coordinates for KDTree
    """
    logging.info(f"Building {name} index from {directory}")
    
    cell_index = {}
    all_coords = []
    
    state_files = sorted(glob(os.path.join(directory, '*.state')))
    
    for state_file in state_files:
        entries = read_state_index(state_file)
        for lon, lat, pos, size in entries:
            coord = (lon, lat)
            if coord not in cell_index:
                cell_index[coord] = (state_file, pos, size)
                all_coords.append([lon, lat])
    
    logging.info(f"  Found {len(cell_index)} unique cells in {len(state_files)} files")
    
    return cell_index, np.array(all_coords)


def build_target_structure(template_dir: str) -> Dict[int, List[Tuple[float, float]]]:
    """
    Read the target grid structure from template state files.
    
    Parameters
    ----------
    template_dir : str
        Directory containing template state files.
        
    Returns
    -------
    Dict[int, List[Tuple[float, float]]]
        Dictionary mapping rank -> list of (lon, lat) coordinates.
    """
    logging.info(f"Reading target structure from {template_dir}")
    
    target_structure = {}
    
    state_files = sorted(glob(os.path.join(template_dir, '*.state')))
    
    for state_file in state_files:
        rank = int(os.path.basename(state_file).replace('.state', ''))
        entries = read_state_index(state_file)
        target_structure[rank] = [(lon, lat) for lon, lat, _, _ in entries]
    
    total_cells = sum(len(cells) for cells in target_structure.values())
    logging.info(f"  Target has {total_cells} cells across {len(target_structure)} ranks")
    
    return target_structure


def find_nearest_neighbor(coord: Tuple[float, float], 
                          tree: cKDTree, 
                          all_coords: np.ndarray) -> Tuple[float, float]:
    """
    Find the nearest neighbor coordinate using spherical approximation.
    
    Parameters
    ----------
    coord : Tuple[float, float]
        Target coordinate (longitude, latitude).
    tree : cKDTree
        KDTree built on transformed coordinates.
    all_coords : np.ndarray
        Array of all source coordinates.
        
    Returns
    -------
    Tuple[float, float]
        Nearest neighbor coordinate (longitude, latitude).
    """
    lon, lat = coord
    
    # Convert to approximate Cartesian for better distance on sphere
    lat_rad = np.radians(lat)
    x = lon * np.cos(lat_rad)
    y = lat
    
    # Transform all coordinates similarly
    all_lat_rad = np.radians(all_coords[:, 1])
    all_x = all_coords[:, 0] * np.cos(all_lat_rad)
    all_y = all_coords[:, 1]
    
    # Build KDTree on transformed coordinates
    transformed = np.column_stack([all_x, all_y])
    
    # Find nearest
    dist, idx = cKDTree(transformed).query([x, y])
    
    return tuple(all_coords[idx])


def write_state_file(output_file: str, 
                     cells_data: List[Tuple[Tuple[float, float], bytes]]):
    """
    Write a new state file with the given cell data.
    
    Parameters
    ----------
    output_file : str
        Path to the output state file.
    cells_data : List[Tuple[Tuple[float, float], bytes]]
        List of ((lon, lat), raw_data) tuples.
    """
    with open(output_file, 'wb') as f:
        # Write cell data and track positions
        index_entries = []
        
        for (lon, lat), data in cells_data:
            pos = f.tell()
            f.write(data)
            index_entries.append((lon, lat, pos))
        
        # Sort index by coordinates (as LPJ-GUESS expects)
        index_entries.sort(key=lambda x: (x[0], x[1]))
        
        # Write index
        for lon, lat, pos in index_entries:
            f.write(struct.pack('<dd', lon, lat))
            f.write(struct.pack('<q', pos))
        
        # Write number of elements
        f.write(struct.pack('<Q', len(index_entries)))


def convert_restart_files_template_first(source_dir: str, 
                                        template_dir: str, 
                                        output_dir: str) -> Dict[str, int]:
    """
    Convert LPJ-GUESS restart files using template-first approach.
    
    Priority strategy:
    1. For cells that exist in both source and target: use source data
    2. For new cells: use template data  
    3. Only use nearest neighbor as last resort
    
    Parameters
    ----------
    source_dir : str
        Directory containing source restart files.
    template_dir : str
        Directory containing template restart files.
    output_dir : str
        Directory where converted files will be written.
        
    Returns
    -------
    Dict[str, int]
        Statistics dictionary with conversion counts.
    """
    import shutil
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"Output directory: {output_dir}")
    
    # Build indices for both source and template
    source_index, source_coords = build_cell_index(source_dir, "source")
    template_index, template_coords = build_cell_index(template_dir, "template")
    
    # Build KDTree for nearest neighbor lookup (from source only)
    logging.info("Building spatial index for nearest neighbor lookup...")
    lat_rad = np.radians(source_coords[:, 1])
    transformed_coords = np.column_stack([
        source_coords[:, 0] * np.cos(lat_rad),
        source_coords[:, 1]
    ])
    tree = cKDTree(transformed_coords)
    
    # Read target structure from template
    target_structure = build_target_structure(template_dir)
    
    # Statistics
    stats = {
        'source_exact_match': 0,
        'template_match': 0,
        'nearest_neighbor': 0,
        'files_written': 0
    }
    nn_distances = []
    
    # Process each target rank
    for rank, target_coords in sorted(target_structure.items()):
        logging.info(f"Processing rank {rank} ({len(target_coords)} cells)...")
        
        cells_data = []
        
        for lon, lat in target_coords:
            coord = (lon, lat)
            
            # Priority 1: Try exact match in source
            if coord in source_index:
                state_file, pos, size = source_index[coord]
                data = read_cell_data(state_file, pos, size)
                cells_data.append((coord, data))
                stats['source_exact_match'] += 1
                
            # Priority 2: Try exact match in template
            elif coord in template_index:
                state_file, pos, size = template_index[coord]
                data = read_cell_data(state_file, pos, size)
                cells_data.append((coord, data))
                stats['template_match'] += 1
                
            # Priority 3: Use nearest neighbor from source (last resort)
            else:
                # Find nearest neighbor
                query_point = np.array([lon * np.cos(np.radians(lat)), lat])
                dist, idx = tree.query(query_point)
                nn_coord = tuple(source_coords[idx])
                
                state_file, pos, size = source_index[nn_coord]
                data = read_cell_data(state_file, pos, size)
                cells_data.append((coord, data))
                stats['nearest_neighbor'] += 1
                
                # Calculate actual spherical distance (approximate)
                nn_lon, nn_lat = nn_coord
                dlat = lat - nn_lat
                dlon = (lon - nn_lon) * np.cos(np.radians((lat + nn_lat) / 2))
                dist_deg = np.sqrt(dlat**2 + dlon**2)
                nn_distances.append(dist_deg)
        
        # Write output file
        output_file = os.path.join(output_dir, f'{rank}.state')
        write_state_file(output_file, cells_data)
        stats['files_written'] += 1
    
    # Copy meta.bin from template (since we're using template structure)
    src_meta = os.path.join(template_dir, 'meta.bin')
    dst_meta = os.path.join(output_dir, 'meta.bin')
    if os.path.exists(src_meta):
        shutil.copy2(src_meta, dst_meta)
        logging.info("Copied meta.bin from template")
    
    return stats, nn_distances


def verify_conversion(output_dir: str, template_dir: str):
    """
    Verify the converted files match expected structure.
    
    Parameters
    ----------
    output_dir : str
        Directory containing converted files.
    template_dir : str
        Directory containing template files for comparison.
    """
    logging.info("")
    logging.info("Verifying converted files...")
    
    # Read output structure
    output_files = sorted(glob(os.path.join(output_dir, '*.state')))
    
    total_cells = 0
    for state_file in output_files:
        entries = read_state_index(state_file)
        total_cells += len(entries)
        
    logging.info(f"Total cells in output: {total_cells}")
    
    # Compare with template
    template_files = sorted(glob(os.path.join(template_dir, '*.state')))
    template_cells = 0
    for state_file in template_files:
        entries = read_state_index(state_file)
        template_cells += len(entries)
    
    logging.info(f"Total cells in template: {template_cells}")
    
    if total_cells == template_cells:
        logging.info("✓ Cell count matches template")
    else:
        logging.warning(f"✗ Cell count mismatch: {total_cells} vs {template_cells}")
    
    # Check meta.bin
    if os.path.exists(os.path.join(output_dir, 'meta.bin')):
        logging.info("✓ meta.bin present")
    else:
        logging.warning("✗ meta.bin missing")


def print_conversion_summary(stats: Dict[str, int], 
                           nn_distances: List[float],
                           output_dir: str,
                           source_index: Dict,
                           template_index: Dict):
    """
    Print detailed summary of conversion results.
    
    Parameters
    ----------
    stats : Dict[str, int]
        Conversion statistics.
    nn_distances : List[float]
        List of nearest neighbor distances.
    output_dir : str
        Directory containing converted files.
    source_index : Dict
        Source cell index for problematic cell analysis.
    template_index : Dict
        Template cell index for problematic cell analysis.
    """
    import numpy as np
    
    logging.info("")
    logging.info("=" * 60)
    logging.info("TEMPLATE-FIRST CONVERSION COMPLETE")
    logging.info("=" * 60)
    logging.info(f"Output directory: {output_dir}")
    logging.info(f"State files written: {stats['files_written']}")
    logging.info(f"Cells from source (DARS2): {stats['source_exact_match']}")
    logging.info(f"Cells from template (CAV-ICB): {stats['template_match']}")
    logging.info(f"Cells from nearest neighbor: {stats['nearest_neighbor']}")
    
    total_cells = stats['source_exact_match'] + stats['template_match'] + stats['nearest_neighbor']
    logging.info(f"Total cells: {total_cells}")
    
    if nn_distances:
        logging.info(f"Nearest neighbor statistics:")
        logging.info(f"  Min distance: {min(nn_distances):.4f} degrees")
        logging.info(f"  Max distance: {max(nn_distances):.4f} degrees")
        logging.info(f"  Mean distance: {np.mean(nn_distances):.4f} degrees")
        logging.info(f"  Median distance: {np.median(nn_distances):.4f} degrees")
    
    # Check for problematic cells
    problematic_cells = [
        (15.8824, 76.8618),
        (49.7368, 80.5151),
        (32.3077, 80.2341),
        (48.4615, 80.2341),
        (53.0769, 80.2341),
        (57.6923, 80.2341),
    ]
    
    logging.info(f"\nProblematic cell analysis:")
    for lon, lat in problematic_cells:
        coord = (lon, lat)
        if coord in source_index:
            source = "SOURCE (would crash)"
        elif coord in template_index:
            source = "TEMPLATE (safe)"
        else:
            source = "NEAREST NEIGHBOR (risky)"
        logging.info(f"  ({lon:.4f}, {lat:.4f}): {source}")