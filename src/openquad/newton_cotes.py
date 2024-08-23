# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

import numpy as np
from scipy.integrate import trapezoid, simpson, romb

from .base import AtomicR1Quadrature, QuadratureWithFixedDegree


class CompositeTrapezoid(AtomicR1Quadrature, QuadratureWithFixedDegree):

    _has_endpoints = True

    def _get_degree(self):
        return 1

    def _points_weights(self, a, b):
        """Composite trapezoid rule.
    
        Composite trapezoid rule for integrating a function on the interval
        :math:`[a,b]`.

        If `period=True` assume peridoic boundary conditions. In this case, the
        function value of point `a` is reused for the point `b`, effectively
        increasing the number of subintervals of :math:`[a,b]` by one.
    
        Parameters
        ----------
        size : int
            Number of equidistant sampling points.
        order : int
            Quadrature order. Equivalent to the number of sampling points.
        a : float
            Lower limit of integration.
        b : float
            Upper limit of integration.
        periodic : logical, optional
            If True, assume periodic boundary conditions.
    
        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Array of ones. This is only returned to retain a common interface.
        """
        N = self.size
        if self._periodic: N += 1
        if N < 2:
            msg = "Trapezoid method requires at least two sample points"
            raise ValueError(msg)
        x, dx = np.linspace(a, b, N, retstep=True)
        w = np.zeros(N) + dx
        w[0] = 0.5 * w[0]
        w[-1] = 0.5 * w[-1]
        return x, w


class CompositeSimpson(AtomicR1Quadrature, QuadratureWithFixedDegree):

    _has_endpoints = True

    def _get_degree(self):
        return 3

    def _points_weights(self, a, b):
        """Composite Simpson's rule.
    
        Composite Simpson's rule for integrating a function on the interval
        :math:`[a,b]`.
    
        Parameters
        ----------
        size : int
            Number of equidistant sampling points.
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
            Array of ones. This is only returned to retain a common interface.
        """
        N = self.size
        if self._periodic: N += 1
        if N < 3:
            msg = "Simpson method requires at least three sample points"
            raise ValueError(msg)
        x, dx = np.linspace(a, b, N, retstep=True)
        i = np.arange(N)
        w = np.zeros(N)
        w[i % 2 == 0] = 2
        w[i % 2 == 1] = 4
        w[0] = 1
        if N % 2 == 1:
            # Even number of subintervals: just do normal Simpson
            w[-1] = 1
        else:
            # Odd number of subintervals: do Simpson until the second to last
            # interval:
            w[-2] = 1
            # and then add the following weights:
            # see https://en.wikipedia.org/wiki/Simpson%27s_rule#Composite_Simpson's_rule_for_irregularly_spaced_data
            w[-3] -= 1 / 4
            w[-2] += 2
            w[-1] = 5 / 4
        w = w * dx / 3
        return x, w


class Romberg(AtomicR1Quadrature, QuadratureWithFixedDegree):

    _has_endpoints = True

    def _get_degree(self):
        # Romberg is in some situations equivalent to Trapz, Simps and Boole,
        # i.e. it can have degree at least up to five. I don't know, if there
        # is a general formula for the degree for higher order of Romb.
        return None #TODO?

    def _points_weights(self, a, b):
        """Romberg integration.

        The number of sample points must be one plus a non-negative power of 2.
    
        This function currently wraps SciPy's `romb` implementation for
        integrating a function on the interval :math:`[a, b]`.
    
        Parameters
        ----------
        size : int
            Number of equidistant sampling points.
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
            Array of ones. This is only returned to retain a common interface.
        """
        # Check that we have a valid number of sample points
        N = self.size
        if self._periodic: N += 1
        n = 1
        while n < N - 1:
            n <<= 1
        if n != N - 1:
            raise ValueError("Number of samples for Romberg integration"
                             "must be one plus a non-negative power of 2")
        x = np.linspace(a, b, N)
        w = np.ones(N)
        return x, w

    def _integrate_sample(self, f_samples, axis=-1):
        dx = self.points[1] - self.points[0]
        f_samples = np.moveaxis(f_samples, axis, -1)
        if self._periodic: # add entry for endpoint of interval
            f_samples = np.concatenate((f_samples, f_samples[...,:1]),
                                       axis=-1)
        return romb(f_samples, dx=dx, axis=-1)
