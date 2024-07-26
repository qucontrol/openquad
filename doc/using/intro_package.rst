.. _intro-structure:

Structure of |openquad|
-----------------------

Quadrature methods in this package are organized according to the geometry of
the domain of integration. |openquad| provides a dedicated class for each
geometry, which form the top-level interface of this data base. All top-level
classes share the same key functionality, amended by geometry-specific
features.

Currently, the following geometries are implemented:

- :math:`n`-dimensional Euclidian spaces, :math:`\mathrm{R}^n`: :class:`openquad.Rn`
- the surface of the 2d unit sphere, :math:`\mathrm{S}^2`: :class:`openquad.S2`
- the space of all orientations in 3d-space, i.e., the rotation group
  :math:`\mathrm{SO}(3)`: :class:`openquad.SO3`

.. todo: link to API reference instead of manually generating this list

