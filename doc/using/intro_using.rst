.. _usage-basics:

Using |openquad|
----------------

This guide describes the basic steps for using the key functionality of
|openquad|, illustrated with a concrete example. If you are not sure how to
use this package, this guide is a good starting point.


Initialization
^^^^^^^^^^^^^^

As a first step, import the class for the desired geometry. For example,
:class:`openquad.S2` for the integral over the surface of the 2d unit sphere,
:math:`\mathrm{S}^2`:

.. testcode::

   >>> import numpy as np
   >>> from openquad import S2

We will also use NumPy later on.

Initiate the class with the method specification (see :class:`~openquad.S2`) in
the form of a list of tuples, containing the a string with the name of the method
and a dictionary with required parameters.

Either use a method tailored to this geometry

.. testcode::

   >>> quad = S2([('S2-Design-Womersley', degree=7)])

or combine lower-dimensional methods with compatible geometries

.. testcode::

   >>> quad = S2([
   ...     ('GaussLegendre', dict(degree=7))]),
   ...     ('Trapezoid', dict(size=8)),
   ... ], polar_sampling='cos')

Some classes take additional parameters, like ``polar_sampling`` for
:attr:`openquad.S2`.


Common functionality
^^^^^^^^^^^^^^^^^^^^

All geometry classes share a set of common attributes and functions,
supplemented by functionality specific to certain geometries.


Sample points and weights
"""""""""""""""""""""""""

Sample points and weights are stored as NumPy arrays in :attr:`~openquad.S2.points` and
:attr:`~openquad.S2.weights`, respecitvely.

.. doctest::

    >>> quad.weights.shape
    (32,)
    >>> quad.points.shape
    (2, 32)

Weights are normalized such that their sum equals the volume of the integration
domain. For :math:`\mathrm{S}^2`, this is the area :math:`4\pi`:

.. doctest::

   >>> np.sum(quad.weights) / (4 * np.pi)
   np.float(1.0)

Sample points in :attr:`~openquad.S2.points` are given in the default
coordinates of the selected integration domain. For :class:`~openquad.S2`,
these are spherical polar angles.  Other coordinates might be available, e.g.
:attr:`~openquad.S2.angles` or :attr:`~openquad.S2.xyz`.

.. doctest::

    >>> np.array_equal(quad.angles, quad.points)
    True
    >>> quad.xyz.shape
    (3, 32)


Exporting quadratures
"""""""""""""""""""""

You can save quadrature points and weights as a textfile with
:meth:`~openquad.S2.savetxt`.

.. testcode::

   >>> quad.savetxt('points_and_weights.txt')


Integration
"""""""""""

Each class is equipped with the :meth:`~openquad.S2.integrate` function, which can handle
arrays and Python callables.

Suppose the integrand :math:`f(x)` is a Python function, e.g.

.. testcode::

    >>> def f(theta, phi):
    ...     """Spherical harmonic |Y_{2, 1}|^2."""
    ...     return ((np.sin(2*theta) * np.cos(phi))*np.sqrt(15/(16*np.pi)))**2

To perform the integral of this function directly

.. doctest::

    >>> quad.integrate(f)
    np.float64(1.0000000000000002)

In some situations it may be desirable or necessary to access the function
values available on the quadrature grid.

.. doctest::

    >>> f_values = f(*quad.angles)
    >>> f_values.shape
    (32,)

You can perform the integration on the array data at a later point with

.. doctest::

    >>> quad.integrate(f_values)
    np.float64(1.0000000000000002)


Other parameters
""""""""""""""""

Other attributes that are available for all top-level classes include:

- :attr:`~openquad.S2.dim`: the dimension of the domain :math:`\mathcal{D}`.
- :attr:`~openquad.S2.size`: the number of sample points.
- :attr:`~openquad.S2.shape`: the shape of :attr:`~openquad.S2.points`.

.. - :attr:`~openquad.S2.source`: original sources of the comprising quadrature methods.

See the :ref:`API reference <api>` for details.

.. todo: give a more precise target of that link
