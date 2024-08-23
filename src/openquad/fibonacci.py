# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

from functools import lru_cache
import numpy as np

from .base import (AtomicSO3Quadrature, AtomicS2Quadrature, QuadratureWithoutDegree)


@lru_cache(maxsize=100)
def fibonacci(i):
    """Return the `i`th member of the Fibonacci series.

    Parameters
    ----------
    i : int
        Index of the Fibonacci series. Must be larger than zero.

    Returns
    -------
    f_i : float
        Fibonacci number at index `i`.

    """
    if i < 0:
        raise ValueError("Index i must not be negative")
    elif i < 2:
        return i
    else:
        return fibonacci(i-1) + fibonacci(i-2)


def _get_available_sizes_for_ZCW2(max_size):
    """Return the available sizes for the ZCW2 method up to `max_size`.

    Note:
    According to :math:`g_M` from Eq. (58) in [1], the ZCW method uses
    for given `M` the `M+6`th Fibonacci number.
    [1]: M. Eden, M.H. Levitt, J. Magn. Reson. 132, 220-239 (1998).

    """
    available_sizes = []
    M_offset = 6
    M = -1
    current_size = -1
    while max_size > current_size:
        M += 1
        current_size = fibonacci(M_offset+M+2)
        available_sizes.append(current_size)
    return np.asarray(available_sizes)


class FibonacciSphere(AtomicS2Quadrature, QuadratureWithoutDegree):

    #TODO: is there a restriction for the available sizes?

    def _points_weights(self):
        """Spherical Fibonacci point sets for S2 quadrature.

        Algorithm taken from Eq. (4) in
        Marques et al., Computer Graphics Forum 32, 8 (2013)
        with an additional modulo operation in order to map the angles to
        the interval $[0, 2\pi]$.

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.

        """
        N = self.size
        j = np.arange(self.size)
        theta = np.arccos(1 - (2*j + 1) / N)
        phi = (4*np.pi * j / (1+np.sqrt(5))) % (2*np.pi)
        points = np.stack((theta, phi), axis=0)
        weights = np.ones(self.size) / self.size * (4*np.pi)
        return points, weights


class ZCW2(AtomicS2Quadrature, QuadratureWithoutDegree):

    _TABLE = _get_available_sizes_for_ZCW2(max_size=1e+7)
    _available_sizes = _TABLE

    def _points_weights(self):
        """ZCW method for S2 quadrature.

        Algorithm taken from
        M. Eden, M.H. Levitt, J. Magn. Reson. 132, 220-239 (1998)

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.

        """
        M_offset = 6
        c = [1, 2, 1] # for integration of full sphere
        N = self.size
        M = self._get_M(N)
        g_M = fibonacci(M + M_offset)
        j = np.arange(self.size)
        theta = np.arccos(c[0] * (c[1] * np.mod(j/N, 1) - 1))
        phi = 2*np.pi / c[2] * np.mod(j*g_M/N, 1)
        points = np.stack((theta, phi), axis=0)
        weights = np.ones(self.size) / self.size * (4*np.pi)
        return points, weights

    def _get_M(self, N):
        """
        Utility function returning the integer :math:`M` corresponding a given
        number of sampling points :math:`N`. See Eq. (59) in
        M. Eden, M.H. Levitt, J. Magn. Reson. 132, 220-239 (1998).
        """
        M_offset = 6
        M = 0
        while fibonacci(M + 2 + M_offset) < N:
            M += 1
        assert fibonacci(M + 2 + M_offset) == N, "Something went wrong..."
        return M


class ZCW3(AtomicSO3Quadrature, QuadratureWithoutDegree):

    def _points_weights(self):
        """ZCW method for SO3 quadrature.

        #TODO: not yet implemented

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.

        """
        raise NotImplementedError("ZCW3 not yet implemented")
        #TODO
