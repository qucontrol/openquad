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

"""
Collection of quadrature grids from various sources.

All data is stored as NumPy arrays in the `.npy` binary format, which can be
loaded by `numpy.save()`.
The data from the original sources have been processed in order to obtain a
unified output format. Details about the format are given below.
For conversion from quaternions to Euler angles, the `quaternionic` package [1]
has been used. A self-written implementation has been used for conversion
between cartesian coordinates and spherical polar angles.

For information on the individual original sources, check out the file `table.py`
in the corresponding subfolder.

All files are formatted in the same way:

- S2 spherical designs:
  NumPy array of dimension `(2, N)` containing the spherical polar angles in
  the order `(polar angle, azimuthal angle)` in the first axis and with `N` the
  number of quadrature points.
  All angles are given in the interval from [0,pi] or [0,2pi], respectively.
- S2 Gauss-type quadratures:
  NumPy array of dimension `(3, N)` containing the spherical polar angles and the weights in
  the order `(polar angle, azimuthal angle, weights)` in the first axis and with `N` the
  number of quadrature points.
- SO3 Chebyshev-type quadratures:
  NumPy array of dimension `(3, N)` containing the Euler angles in the z-y-z convention
  the order `(alpha, beta, gamma)` in the first axis and with `N` the
  number of quadrature points.
- SO3 Gauss-type quadratures:
  NumPy array of dimension `(4, N)` containing the Euler angles in the z-y-z convention and the weights
  the order `(alpha, beta, gamma, weights)` in the first axis and with `N` the
  number of quadrature points.

All angles are given in the interval from [0,pi] or [0,2pi], respectively.


[1]: https://github.com/moble/quaternionic
"""
