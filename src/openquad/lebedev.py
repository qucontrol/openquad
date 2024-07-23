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

from .base import AtomicS2Quadrature, QuadratureWithDegree
from .data.lebedev import LEBEDEV_LAIKOV_QUADRATURES


class LebedevLaikov(AtomicS2Quadrature, QuadratureWithDegree):

    _TABLE = np.asarray(LEBEDEV_LAIKOV_QUADRATURES, dtype=int)
    _available_sizes = _TABLE[:, 1]
    _available_degrees = _TABLE[:, 0]

    def _points_weights(self):
        """Lebedev-Laikov quadrature.

        We simply wrap sphericalquadpy's Lebedev class for now.

        #TODO
        Parameters
        ----------
        size : int
            Number of equidistant sampling points.
        order : int
            Quadrature order.

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Array of ones. This is only returned to retain a common interface.

        """
        points, weights = self._load_points_weights('lebedev', 's2_gauss')
        weights = weights * (4*np.pi)
        return points, weights
