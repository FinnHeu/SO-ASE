# so_ase/__init__.py

__version__ = "0.1.0"

from .plotting_maps import create_map
from .eval_sea_ice import fesom_ice_area, fesom_ice_volume
from .eval_ocean import fesom_ocean_heat_transport_as_residual