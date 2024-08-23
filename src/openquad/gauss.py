# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

import numpy as np
from scipy import special

from .base import AtomicR1Quadrature, QuadratureWithDegree


class GaussLegendre(AtomicR1Quadrature, QuadratureWithDegree):

    _has_endpoints = False

    def _get_size(self, degree):
        return int(np.ceil((degree + 1) / 2))

    def _get_degree(self, size):
        return 2 * size - 1

    def _points_weights(self, a, b):
        """Calculate Gauss-Legendre quadrature points and weights.
    
        This function returns the quadrature points and weights suitable for
        integrating a function on the interval :math:`[a, b]`. The Jacobian for
        the transformation from the interval :math:`[-1, 1]` to :math:`[a, b]`
        is already included in the weights.
    
        Parameters
        ----------
        size : int
            Number of sampling points.
        order : int
            Quadrature order. Equivalent to the number of sampling points.
        a : float
            Lower limit of integration.
        b : float
            Upper limit of integration.
    
        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.
    
        """
        x, w = special.roots_legendre(self.size)
        return (b-a)/2 * x + (a+b)/2, w * (b-a)/2


class GaussLobattoLegendre(AtomicR1Quadrature, QuadratureWithDegree):

    _has_endpoints = True

    def _get_size(self, degree):
        size = int(np.ceil((degree + 3) / 2))
        if self._periodic:
            size = size - 1
        return size

    def _get_degree(self, size):
        if self._periodic:
            size = size + 1
        return 2 * size - 3

    def _points_weights(self, a, b):
        """Calculate Gauss-Lobatto-Legendre quadrature points and weights.

        This function returns the quadrature points and weights suitable for
        integrating a function on the interval :math:`[a, b]`. The Jacobian for
        the transformation from the interval :math:`[-1, 1]` to :math:`[a, b]`
        is already included in the weights.

        Parameters
        ----------
        size : int
            Number of sampling points.
        order : int
            Quadrature order. Equivalent to the number of sampling points.
        a : float
            Lower limit of integration.
        b : float
            Upper limit of integration.

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.

        """
        N = self.size
        if N < 2:
            msg = ("Gauss-Lobatto rule requires at least two sample points")
            raise ValueError(msg)
        coeffs = np.append(np.zeros(N-1), 1)
        x = np.concatenate(([-1.0],
                            np.polynomial.legendre.Legendre(coeffs).deriv().roots(),
                            [1.0]))
        w = 2 / (N*(N-1) * special.eval_legendre(N-1, x)**2)
        return (b-a)/2 * x + (a+b)/2, w * (b-a)/2
