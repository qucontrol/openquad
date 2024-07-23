
__version__ = "alpha+dev"

# Expose top-level interface
from .geometries import Rn, S2, SO3

# Expose submodules
from . import (
    geometries,
    newton_cotes,
    gauss,
    lebedev,
    monte_carlo,
    grid,
    graef,
    womersley,
    karney,
    fibonacci,
)
