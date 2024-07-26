.. _howtos:

How-tos
-------

Recipies to achieve common tasks using |openquad|, targeted towards users familiar with the basic concepts.

For detailed description of the classes and methods contained in |openquad|, see the :ref:`API reference <api>`.

.. todo: move the testsetup to conftest.py
.. testsetup::

   >>> import numpy as np


How to access quadrature points and weights
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Obtain Gauss-Legendre sample points and weights for degree `71` on the interval
`[-10, 5]`:

.. testcode::

    >>> from openquad import Rn
    >>> quad = Rn([
    ...     ('GaussLegendre', dict(degree=21, a=-10, b=5)),
    ... ])

Quadrature points and weights are stored in :attr:`~openquad.Rn.points` and :attr:`~openquad.Rn.weights`,

.. doctest::

    >>> quad.points.shape
    (1, 11)
    >>> quad.weights.shape
    (11,)


How to integrate a python function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Given you have a python function defined on the unit sphere in terms of
spherical polar coordinates.

.. testcode::

    >>> def func(theta, phi):
    ...     """Spherical harmonic |Y_{2, 1}|^2."""
    ...     return ((np.sin(2*theta) * np.cos(phi))*np.sqrt(15/(16*np.pi)))**2
 
For example, initiate a :math:`\mathrm{S}^2` `spherical design`_ of degree `7`:

.. _spherical design: https://en.wikipedia.org/wiki/Spherical_design

.. testcode::

    >>> from openquad import S2
    >>> quad = S2([
    ...     ('S2-Design-Graef', dict(degree=7)),
    ... ])

To evaluate the integral of ``func`` over :math:`\mathrm{S}^2` use the :meth:`~openquad.S2.integrate` method:
    
.. doctest::

    >>> quad.integrate(func)
    np.float64(0.9999999999975318)


How to export quadrature points and weights
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For example, create a quadrature method for an integral over the three `Euler angles`_,
using `Lebedev-Laikov quadrature`_ of degree `5` for the first two angles
combined with the composite trapezoid rule with `6` sample points for the third
angle. 

.. testcode::

    >>> from openquad import SO3
    >>> quad = SO3([
    ...     ('LebedevLaikov', dict(degree=5)),
    ...     ('Trapezoid', dict(size=6)),
    ... ])
    
Save sample points and weights to a text with the :meth:`~openquad.SO3.savetxt` method.

.. testcode::

    >>> quad.savetxt('points_and_weights.dat')

.. _Lebedev-Laikov quadrature: https://en.wikipedia.org/wiki/Lebedev_quadrature
.. _Euler angles: https://en.wikipedia.org/wiki/Euler_angles
