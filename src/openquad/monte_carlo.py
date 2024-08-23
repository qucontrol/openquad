# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

import numpy as np
import quaternionic

from .base import QuadratureWithoutDegree, AtomicR1Quadrature, AtomicS2Quadrature, AtomicSO3Quadrature


class MonteCarloR1(AtomicR1Quadrature, QuadratureWithoutDegree):

    _has_endpoints = False

    def __init__(self, size, a, b, jacobian=None, seed=None, **kwargs):
        super().__init__(size=size, a=a, b=b, jacobian=jacobian, seed=seed)

    def _points_weights(self, a, b, seed):
        """One-dimensional Monte-Carlo integration.
    
        Parameters
        ----------
        size : int
            Number of random sampling points.
    
        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Array of ones. This is only returned to retain a common interface.

        """
        rng = np.random.default_rng(seed=seed)
        x = rng.uniform(a, b, self.size)
        w = np.ones(self.size, dtype=float) / self.size * (b-a)/2
        return x, w


class MonteCarloS2(AtomicS2Quadrature, QuadratureWithoutDegree):

    def __init__(self, size, seed=None):
        super().__init__(size=size, seed=seed)

    def _points_weights(self, seed):
        """Monte-Carlo integration on the unit sphere.
    
        Parameters
        ----------
        size : int
            Number of random sampling points.
    
        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Array of ones. This is only returned to retain a common interface.

        """
        x = self._sample_spherical_polar_angles(self.size, seed)
        w = np.ones(self.size, dtype=float) / self.size * (4*np.pi)
        return x, w

    def _sample_spherical_polar_angles(self, n, seed=None):
        """Generate a uniform distribution of points on the unit sphere.

        The unit sphere is sampled in spherical polar coodinates.
        
        Parameters
        ----------
        n : int
            Number of points to generate.
        seed : int, optional
            Seed for the pseudo-random number generator. If None, it will
            be determined by the OS.
    
        Returns
        -------
        angles : ndarray
            Array of spherical polar angles of shape ``(2, n)``.
            The angles are stored in the first axis, with the polar angle
            coming first.

        """
        rng = np.random.default_rng(seed=seed)
        phi = rng.uniform(0, 2*np.pi, n)
        u = rng.uniform(-1, 1, n)
        theta = np.arccos(u)
        return np.array([theta, phi])


class MonteCarloSO3(AtomicSO3Quadrature, QuadratureWithoutDegree):

    def __init__(self, size, seed=None, method='quaternions'):
        super().__init__(size=size, seed=seed, method=method)

    def _points_weights(self, seed, method):
        """Three-angle Monte-Carlo integration.
    
        Parameters
        ----------
        size : int
            Number of random sampling points.
        method : str
            Sampling method. Possible options are `'quaternions'` and
            `'angles'`.
    
        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Array of ones. This is only returned to retain a common interface.

        """
        x = self.random_rotations(self.size, method, seed)
        w = np.ones(self.size, dtype=float) / self.size * (8*np.pi**2)
        return x, w

    def random_rotations(self, n, method, seed=None):
        """Generate uniform random 3d rotations.
        
        The parameter `method` controls which method is used to sample
        orientation space: For ``method='quaternions'``, orientation space is
        sampled in quaternion representation with Showmake's method.  For
        ``method='angles'``, orientation space is sampled in Euler angle
        representation. Both methods yield uniform sampling without clustering
        at the poles.
    
        Parameters
        ----------
        n : int
            Number of rotations to generate.
        method : str
            Sampling method.
        seed : int, optional
            Seed for the pseudo-random number generator. If None, it will
            be determined by the OS.
    
        Returns
        -------
        angles : ndarray
            Array of Euler angles of shape ``(3, n)``.
            The Euler angles are stored in the first axis.

        """
        if method == 'quaternions':
            quaternions = self._sample_quaternions(n, seed)
            x = quaternions.to_euler_angles.T
        elif method == 'angles':
            x = self._sample_euler_angles(n, seed)
        else:
            raise ValueError(f"Invalid sampling method: {method}")
        return x

    def _sample_quaternions(self, n, seed=None):
        """Generate uniform random 3d rotations in quaternion representation.
        
        This function uses Shoemake's method.
    
        Parameters
        ----------
        n : int
            Number of rotations to generate.
        seed : int, optional
            Seed for the pseudo-random number generator. If None, it will
            be determined by the OS.
    
        Returns
        -------
        quaternsions : quaternionic array
            Array of length `n` with 3d rotations in quaternion representation.
            The quaternion elements are stored in the last axis.

        """
        rng = np.random.default_rng(seed=seed)
        x = rng.random((3, n))
        s1 = np.sin(2*np.pi*x[0])
        s2 = np.sin(2*np.pi*x[1])
        c1 = np.cos(2*np.pi*x[0])
        c2 = np.cos(2*np.pi*x[1])
        r1 = np.sqrt(1 - x[2])
        r2 = np.sqrt(x[2])
        return quaternionic.array(np.transpose([c2*r2, s1*r1, c1*r1, s2*r2]))

    def _sample_euler_angles(self, n, seed=None):
        """Generate uniform random 3d rotations in Euler angle representation.
        
        Parameters
        ----------
        n : int
            Number of rotations to generate.
        seed : int, optional
            Seed for the pseudo-random number generator. If None, it will
            be determined by the OS.
    
        Returns
        -------
        angles : ndarray
            Array of Euler angles of shape ``(3, n)``.
            The Euler angles are stored in the first axis.

        """
        rng = np.random.default_rng(seed=seed)
        alpha = rng.uniform(0, 2*np.pi, n)
        gamma = rng.uniform(0, 2*np.pi, n)
        u = rng.uniform(-1, 1, n)
        beta = np.arccos(u)
        return np.array([alpha, beta, gamma])

    @staticmethod
    def random_rotation(seed=None):
        """Generate a random 3d rotation in quaterion representation.

        Equivalent to ``np.squeeze(random_rotations(1, seed))``.

        Parameters
        ----------
        seed : int, optional
            Seed for the pseudo-random number generator. If None, it will
            be determined by the OS.

        Returns
        -------
        q : quaternion
            A single unit quaternion.

        """
        return np.squeeze(random_rotations(1, seed))
