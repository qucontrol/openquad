.. _usage-basics:

Using |openquad|
----------------

This guide describes the basic steps for using the key functionality of
|openquad|, illustrated with a concrete example. If you are not sure how to
use this package, this guide is a good starting point.


Initialization
^^^^^^^^^^^^^^

As a first step, import the class for the desired geometry. For example,
:class:`S2` for the integral over the surface of the 2d unit sphere,
:math:`\mathrm{S}^2`:

.. testcode::

   >>> import numpy as np
   >>> from openquad import S2

We will also use NumPy later on.

Initiate the class with the method specification (see :class:`S2`) in
the form of a list of tuples, containing the name of the method and
required parameters.

Either use a method tailored to this geometry

.. testcode::

   >>> quad = S2([('S2-Design-Womersley', degree=7)])

or combine lower-dimensional methods with compatible geometries

.. testcode::

   >>> quad = S2([
   ...     ('GaussLegendre', degree=7)]),
   ...     ('Trapezoid', size=8),
   ... ], polar_sampling='cos')

Some classes take additional parameters, like :attr:`~S2.polar_sampling`.


Common functionality
^^^^^^^^^^^^^^^^^^^^

All geometry classes share a set of common attributes and functions,
supplemented by functionality specific to certain geometries.


Sample points and weights
"""""""""""""""""""""""""

Sample points and weights are stored as NumPy arrays in :attr:`~S2.points` and
:attr:`~S2.weights`, respecitvely.

.. doctest::

    >>> quad.weights.shape
    (28,)
    >>> quad.points.shape
    (2, 28)

Weights are normalized such that their sum equals the volume of the integration
domain. For :math:`\mathrm{S}^2`, this is the area :math:`4\pi`:

.. doctest::

   >>> np.sum(quad.weights) / (4 * np.pi)
   1.0

Sample points in :attr:`~S2.points` are given in the default coordinates of the
selected integration domain. For :class:`S2`, these are spherical polar angles.
Other coordinates might be available, e.g. :attr:`~S2.angles` or :attr:`~S2.xyz`.

.. doctest::

    >>> np.array_equal(quad.angles, quad.points)
    True
    >>> quad.xyz.shape
    (3, 28)


Exporting quadratures
"""""""""""""""""""""

You can save quadrature points and weights as a textfile with
:meth:`savetxt`.

.. testcode::

   >>> quad.savetxt('points_and_weights.txt')


Integration
"""""""""""

Each class is equipped with the :meth:`integrate` function, which can handle
arrays and Python callables.

Suppose the integrand :math:`f(x)` is a Python function, e.g.

.. testcode::

    >>> def f(theta, phi):
    ...     return np.cos(theta) * np.sin(phi)

To perform the integral of this function directly

.. doctest::

    >>> quad.integrate(f)
    1.0

In some situations it may be desirable or necessary to access the function
values available on the quadrature grid.

.. doctest::

    >>> f_values = f(*quad.angles)
    >>> f_values.shape
    (2, 28)

You can perform the integration on the array data at a later point with

.. doctest::

    >>> quad.integrate(f_values)
    1.0


Other parameters
""""""""""""""""

Other attributes that are available for all top-level classes include:

- :attr:`~S2.dim`: the dimension of the domain :math:`\mathcal{D}`.
- :attr:`~S2.size`: the number of sample points.
- :attr:`~S2.shape`: the shape of :attr:`~S2.points`.

.. - :attr:`~S2.source`: original sources of the comprising quadrature methods.

See the :ref:`API reference <api>` for details.

.. todo: give a more precise target of that link
