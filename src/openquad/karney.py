# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

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
