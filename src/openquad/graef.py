import numpy as np

from .base import (AtomicSO3Quadrature, AtomicS2Quadrature, QuadratureWithPolyAcc)
from .data.graef.table import (
    S2_GAUSS_METHODS,
    S2_SPHERICAL_DESIGNS,
    SO3_GAUSS_METHODS,
    SO3_EQUAL_WEIGHT_METHODS,
)


class GraefS2Gauss(AtomicS2Quadrature, QuadratureWithPolyAcc):

    _TABLE = np.asarray(S2_GAUSS_METHODS, dtype=int)
    _available_sizes = _TABLE[:, 1]
    _available_p_accs = _TABLE[:, 0]

    def _points_weights(self):
        """FFT optimized Gauss-type S2 quadrature schemes from Manuel Gräf.

        Quadrature points and weights taken from:
        https://www-user.tu-chemnitz.de/~potts/workgroup/graef/quadrature/

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.

        """
        points, weights = self._load_points_weights('graef', 's2_gauss')
        weights = weights * (4*np.pi)
        return points, weights


class GraefS2Design(AtomicS2Quadrature, QuadratureWithPolyAcc):

    _TABLE = np.asarray(S2_SPHERICAL_DESIGNS, dtype=int)
    _available_sizes = _TABLE[:, 1]
    _available_p_accs = _TABLE[:, 0]

    def _points_weights(self):
        """FFT optimized S2 spherical designs from Manuel Gräf.

        Quadrature points and weights taken from:
        https://www-user.tu-chemnitz.de/~potts/workgroup/graef/quadrature/

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.

        """
        points, weights = self._load_points_weights('graef', 's2_design')
        weights = weights * (4*np.pi)
        return points, weights


class GraefSO3Gauss(AtomicSO3Quadrature, QuadratureWithPolyAcc):

    _TABLE = np.asarray(SO3_GAUSS_METHODS, dtype=int)
    _available_sizes = _TABLE[:, 1]
    _available_p_accs = _TABLE[:, 0]

    def _points_weights(self):
        """FFT optimized Gauss-type SO3 quadrature schemes from Manuel Gräf.

        Quadrature points and weights taken from:
        https://www-user.tu-chemnitz.de/~potts/workgroup/graef/quadrature/

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.

        """
        points, weights = self._load_points_weights('graef', 'so3_gauss')
        weights = weights * (8*np.pi**2)
        return points, weights


class GraefSO3EqualWeight(AtomicSO3Quadrature, QuadratureWithPolyAcc):

    _TABLE = np.asarray(SO3_EQUAL_WEIGHT_METHODS, dtype=int)
    _available_sizes = _TABLE[:, 1]
    _available_p_accs = _TABLE[:, 0]

    def _points_weights(self):
        """FFT optimized equal-weight SO3 quadrature schemes from Manuel Gräf.

        Quadrature points and weights taken from:
        https://www-user.tu-chemnitz.de/~potts/workgroup/graef/quadrature/

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.

        """
        points, weights = self._load_points_weights('graef', 'so3_equalweight')
        weights = weights * (8*np.pi**2)
        return points, weights
