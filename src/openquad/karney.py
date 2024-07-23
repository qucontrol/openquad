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

import numpy as np

from .base import QuadratureWithoutDegree, AtomicSO3Quadrature
from .data.karney import SO3_COVERINGS


class KarneySO3(AtomicSO3Quadrature, QuadratureWithoutDegree):

    _TABLE = np.array(SO3_COVERINGS, dtype='object')
    _available_sizes = _TABLE[:, 0].astype(int)

    def _points_weights(self):
        """Nearly optimal SO3 coverings from Charles F.F. Karney.

        These coverings are based on 4d polytopes.

        Coverings taken from:
        https://github.com/cffk/orientation

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.

        """
        points, weights = self._load_points_weights('karney', 'so3_covering')
        weights = weights * (8*np.pi**2)
        return points, weights
