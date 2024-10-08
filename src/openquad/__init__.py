# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

__version__ = '0.3.0'

# Expose top-level interface
from .geometries import Rn, S2, SO3

__all__ = ['Rn', 'S2', 'SO3']

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
