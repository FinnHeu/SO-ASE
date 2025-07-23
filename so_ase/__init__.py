# so_ase/__init__.py

__version__ = "0.1.0"

from .plotting_maps import create_map
from .eval_sea_ice import fesom_ice_area, fesom_ice_volume
from .eval_ocean import fesom_ocean_heat_transport_as_residual, fesom_timeseries_of_mean_vertical_profile_in_region, regression2D_fesom
from .helpers_mesh import read_nodes, read_elements, read_aux3d, find_nodes_in_box
from .helpers_regression import dataset_regression_on_time_1D
