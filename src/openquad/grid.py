# OpenQuad - open database for multi-dimensional numerical integration
# SPDX-FileCopyrightText: 2024  Alexander Blech
# SPDX-License-Identifier: MPL-2.0

from abc import ABC, abstractmethod
import numpy as np
import quaternionic


def angles_from_xyz(xyz, axis=0):
    """Calculate the spherical polar angles of a set of cartesian
    coordinates.

    $\theta$ denotes the polar angle, while $\phi$ denotes the azimuthal angle.
    
    Parameters
    ----------
    xyz : array_like
        Array with cartesian coordinates in the order $(x, y, z)$.
    axis : int
        Specifies the axis, along which the coordinates are stored.
        By default, the first axis.

    Returns
    -------
    theta_phi : ndarray
        Array of polar and azimuthal angles in the order $(\theta, \phi)$ stored in the axis
        given by `axis`.
    """
    # We assume that we have points on a unit sphere:
    assert np.allclose(np.linalg.norm(xyz, axis=axis), 1.0)
    xyz = np.moveaxis(xyz, axis, 0)
    x, y, z = xyz
    θ = np.arccos(z / 1.0)
    ϕ = np.arctan2(y, x) % (2*np.pi)
    return np.moveaxis([θ, ϕ], 0, axis)


def xyz_from_angles(theta_phi, axis=0):
    """Calculate the cartesian coordinates from a set of spherical
    polar angles.

    $\theta$ denotes the polar angle, while $\phi$ denotes the azimuthal angle.
    
    Parameters
    ----------
    theta_phi : array_like
        Array with polar and azimuthal angles in the order $(\theta, \phi)$.
    axis : int
        Specifies the axis, along which the values for $\theta$ and $\phi$ are stored.
        By default, the first axis.

    Returns
    -------
    xyz : ndarray
        Array of cartesian coordinates in the order $(x, y, z)$ stored in the axis
        given by `axis`.
    """
    theta_phi = np.moveaxis(theta_phi, axis, 0)
    θ, ϕ = theta_phi
    x = np.sin(θ) * np.cos(ϕ)
    y = np.sin(θ) * np.sin(ϕ)
    z = np.cos(θ)
    return np.moveaxis([x, y, z], 0, axis)


def euler_angles_from_quaternions(quaternions, axis=0):
    """Convert quaternions to Euler angles.

    Paramters
    ---------
    angles : quaternion array
        Array of quaternions
    axis : int, optional
        Axis that contains the quaternion elements. By default, use the first
        axis.

    Returns
    -------
    angles : ndarray
        Array of Euler angles. The Euler angles are stored in axis
        `axis`.

    """
    quaternions = np.moveaxis(quaternions, axis, -1)
    angles = quaternions.as_euler_angles % (2*np.pi)
    return np.moveaxis(angles, -1, axis)


def quaternions_from_euler_angles(angles, axis=0):
    """Convert Euler angles to quaternions.

    Paramters
    ---------
    angles : array_like
        Array of Euler angles.
    axis : int, optional
        Axis that contains the Euler angles. By default, use the first axis.

    Returns
    -------
    quaternions : quaternionic array
        Array of quaternions. The quaternion elements are stored in axis
        `axis`.

    """
    angles = np.moveaxis(angles, axis, 0)
    quaternions = quaternionic.array.from_euler_angles(*angles)
    return np.moveaxis(quaternions, -1, axis)


class Grid(ABC):
    """TODO"""
    
    @abstractmethod
    def _get_file_data_and_header(self, representation):
        """TODO"""
        data = None
        header = None
        return data, header

    # TODO: add writing a binary format
    def save(self, filename, representation, weights=None):
        """TODO"""
        data, header = self._get_file_data_and_header(representation)
        if weights is not None:
            data = np.vstack((*data, weights))
            header = header + f"{'weights':>25}"
        index = np.arange(1, data.shape[-1]+1)
        fmt = '%10i' + data.shape[0] * '%25.16e'
        np.savetxt(filename, np.vstack((index, *data)).T,
                   fmt=fmt, header=header)


class S2Grid(Grid):
    """Grid of points on a sphere.

    Grid points are stored either as spherical polar angles `theta` and `phi`
    or in terms of cartesian coordinates `x`, `y` and `z`.
    """
    def __init__(self, representation, points, axis=0):
        """
        Parameters
        ----------
        representation : str
            Representation of the points in `points`.
            Currently supported: `'angles'`, `'xyz'`.
        points : array-like
            Array of point coordinates in a given representation.
        axis : int
            The axis of `points` that contains the coordinates on input.
            Eventually, the coordinates will be stored in the first axis.
        """
        self.original_representation = representation
        self.original_shape = points.shape
        self.original_axis = axis
        self.N = len(points)
        self._angles = None
        self._xyz = None

        points = np.moveaxis(points, axis, 0)
        if representation == 'angles':
            self.angles = points.reshape(2, -1)
        elif representation == 'xyz':
            self.xyz = points.reshape(3, -1)
        else:
            raise ValueError(f"Unknown representation: {representation}")

    @property
    def angles(self):
        """Array of spherical polar angles."""
        if self._angles is None:
            self._angles = angles_from_xyz(self._xyz, axis=0)
        return self._angles

    @angles.setter
    def angles(self, angles):
        self._angles = angles

    @property
    def xyz(self):
        """Array of cartesian coordinates."""
        if self._xyz is None:
            self._xyz = xyz_from_angles(self._angles, axis=0)
        return self._xyz

    @xyz.setter
    def xyz(self, xyz):
        self._xyz = xyz

    def _get_file_data_and_header(self, representation):
        header = (
            "Spherical grid"
            + f"\nRepresentation: {representation}\n"
            + f"Number of grid points: {self.N}\n"
        )
        if representation == 'angles':
            data = self.angles
            header = header + f"{'index':>8}{'theta (rad)':>25}{'phi (rad)':>25}"
        elif representation == 'xyz':
            data = self.xyz
            header = header + f"{'index':>8}{'x':>25}{'y':>25}{'z':>25}"
        else:
            raise ValueError(f"Unknown representation: {representation}")
        return data, header


class SO3Grid(Grid):
    """Grid of orientation or rotations.

    Orientations are stored in several representations:
    Euler angles, quaternions.
    """
    def __init__(self, representation, points, axis=0):
        """
        Parameters
        ----------
        representation : str
            Representation of the orientations in `points`.
            Currently supported: `'Euler angles'`, `'quaternions'`.
        points : array-like
            Array of orientations in a given representation.
        axis : int
            The axis of `points` that contains the orientations on input.
            Eventually, the coordinates will be stored in the first axis.
        """
        points = np.moveaxis(points, axis, 0)
        self.original_representation = representation
        self.original_shape = points.shape
        self.original_axis = axis
        self.N = points.shape[-1]
        self._angles = None
        self._quaternions = None

        if representation == 'Euler angles':
            self.angles = points.reshape(3, -1)
        elif representation == 'quaternions':
            self.quaternions = points.reshape(-1, 4)
        else:
            raise ValueError(f"Unknown representation: {representation}")

    @property
    def angles(self):
        """Array of Euler angles."""
        if self._angles is None:
            self._angles = euler_angles_from_quaternions(self._quaternions)
        return self._angles

    @angles.setter
    def angles(self, angles):
        self._angles = angles

    @property
    def quaternions(self):
        """Array of quaternions."""
        if self._quaternions is None:
            self._quaternions = quaternions_from_euler_angles(self._angles)
        return self._quaternions

    @quaternions.setter
    def quaternions(self, quaternions):
        self._quaternions = quaternions

    def _get_file_data_and_header(self, representation):
        header = (
            "Orientation grid"
            + f"\nRepresentation: {representation}\n"
            + f"Number of grid points: {self.N}\n"
        )
        if representation == 'Euler angles':
            data = self.angles
            header = header + f"{'index':>8}{'alpha (rad)':>25}{'beta (rad)':>25}{'gamma (rad)':>25}"
        elif representation == 'quaternions':
            data = self.quaternions
            header = header + f"{'index':>8}{'q_0':>25}{'q_1':>25}{'q_2':>25}{'q_3':>25}"
        else:
            raise ValueError(f"Unknown representation: {representation}")
        return data, header
