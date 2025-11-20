# so_ase/__init__.py

__version__ = "0.1.0"

from .plotting_maps import create_map, plot_on_elements
from .eval_sea_ice import fesom_ice_area, fesom_ice_volume, nsidc_ice_area, hadlsst_ice_area
from .eval_ocean import fesom_ocean_heat_transport_as_residual, fesom_timeseries_of_mean_vertical_profile_in_region
from .helpers_mesh import read_nodes, read_elements, read_aux3d, read_cavity_depth_at_node, read_element_levels, find_nodes_in_box, reproject_to_latlon
from .helpers_regression import regression2D_fesom, regression2D_nsidc, regression2D_hadlsst, anomaly2D_fesom, dataset_regression_on_time_1D
from .helpers_plots import add_notebook_path_to_fig
