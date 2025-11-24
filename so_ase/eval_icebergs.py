# so_ase.eval_icebergs.py

import numpy as np

def read_iceberg_initial_files(icebergpath):

    def read_dat(filepath):
        data = []
        with open(filepath) as f:
            for line in f:
                data.append(f.readline())
        return np.array(data, dtype=float)
    
    icb_lon = read_dat(icebergpath + "/icb_lon.dat")
    icb_lat = read_dat(icebergpath + "/icb_lat.dat")
    icb_height = read_dat(icebergpath + "/icb_height.dat")
    icb_felem = read_dat(icebergpath + "/icb_felem.dat").astype(int)
    icb_length = read_dat(icebergpath + "/icb_length.dat")
    icb_scaling = read_dat(icebergpath + "/icb_scaling.dat")

    return icb_lon, icb_lat, icb_height, icb_felem, icb_length, icb_scaling