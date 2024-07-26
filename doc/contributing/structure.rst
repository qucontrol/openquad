.. _package-structure:

Structure of this package
-------------------------

The top level interface, located in :mod:`openquad.geometries`, provides
classes for quadratures on different geometries, e.g. a unit sphere, the
rotation group, n-dimensional cartesian space and other. They are derived from
the abstract class :class:`GeometryQuadrature` defining their main attributes
and methods.

On a lower level, individual quadrature methods are organized in modules
corresponding to their respective source. They are implemented in terms of
classes derived from the abstract base classes in :mod:`base`. On the one hand,
these base classes define the basic properties for quadratures on different
geometries. On the other hand, they distinguish between quadratures with and
without non-zero degree of exactness.

Quadrature points and weights are located in the ``data`` directory, organized in
separate folders for the different sources. Each subfolder contains a file
``__init__.py`` giving information about the source and listing available
degrees and sizes. Actual coordinates and weights are stored in NumPy's
``.npy`` `binary file format`_ in the following convention, ``n`` denoting the
number of sample points:

- S2 spherical designs:
  Array of shape ``(2, n)``, along the first axis containing the spherical polar
  angles in the order ``[polar angle, azimuthal angle]``.
- S2 Gauss-type quadratures:
  Array of shape ``(3, n)`` along the first axis containing the spherical polar
  angles and the weights in the order ``[polar angle, azimuthal angle,
  weights]``.
- SO3 equal-weight quadratures:
  Array of shape ``(3, n)`` along the first axis containing the Euler angles in
  the z-y-z convention the order ``[alpha, beta, gamma]``.
- SO3 Gauss-type quadratures:
  Array of shape ``(4, n)`` along the first axis containing the Euler angles in
  the z-y-z convention and the weights the order ``[alpha, beta, gamma,
  weights]``.

All angles are given in the interval from :math:`[0,\pi]` or :math:`[0,2\pi]`,
respectively.

.. _binary file format: https://numpy.org/doc/stable/reference/generated/numpy.lib.format.html
