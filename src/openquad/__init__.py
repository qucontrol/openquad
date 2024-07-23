# OpenQuad - open database for multi-dimensional numerical integration
# Copyright (C) 2024  Alexander Blech
#
# This library is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this library. If not, see <https://www.gnu.org/licenses/>.

__version__ = "0.1.0-dev"

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
