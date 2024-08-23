# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

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
