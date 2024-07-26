.. _quickstart:

Quickstart guide
----------------


What is this package?
^^^^^^^^^^^^^^^^^^^^^

Combining a comprehenseive collection of quadrature methods in one single
interface, |openquad| helps you to efficiently compute integrals over various
domains with ease, including 1d intervals, spherical surfaces, Euler angles and
more. The object-oriented interface provides Python classes for each geometry,
for example the :class:`openquad.S2` class for integration over the 2d unit
sphere.  Python functions and array data can be integrated directly using the
:meth:`~openquad.S2.integrate` method. Alternatively, sample points and
quadrature weights can be exported with the :meth:`~openquad.S2.savetxt` method
for easy integration with other software.

The best part: combine individual quadrature methods to create custom
multi-dimensional quadratures.


How do I use it?
^^^^^^^^^^^^^^^^

.. todo: testsetup into conftest.py
.. testsetup::

   >>> import numpy as np
   >>> def your_awesome_function(arg0, arg1, arg2=None):
   ...     return arg0*arg1

Want to integrate a Python function over the interval :math:`[-10,5]` with
Gauss-Legendre quadrature?

.. testcode::

    >>> from openquad import Rn
    >>> quad = Rn([
    ...     ('GaussLegendre', dict(degree=21, a=-10, b=5)),
    ... ])
    >>> quad.integrate(your_awesome_function)

Have a Python function that you need to evaluate on a uniform grid over the
surface of the unit sphere?

.. testcode::
   :hide:

   >>> def your_awesome_function(theta, phi):
   ...  return np.sin(theta) * np.cos(phi)

.. testcode::

    >>> from openquad import S2
    >>> quad = S2([
    ...     ('S2-Covering-Fibonacci', dict(size=142)),
    ... ])
    >>> grid = your_awesome_function(*quad.angles)

Need quadrature points and weights for performing an efficient orientation
average in your molecular dynamics simulation?

.. testcode::

    >>> from openquad import SO3
    >>> quad = SO3([
    ...     ('LebedevLaikov', dict(degree=15)),
    ...     ('Trapezoid', dict(size=16)),
    ... ])
    >>> quad.savetxt('points_and_weights.dat')


How do I get it?
^^^^^^^^^^^^^^^^

Install the latest stable release with

.. code-block:: bash

    pip install openquad


Tell me more!
^^^^^^^^^^^^^

Have a look at the :ref:`user guide <using>` and the :ref:`example gallery
<examples>`.

.. tip::

   Don't know, which method to choose? Read our `paper`_.

.. _paper: https://arxiv.org/abs/2407.17434
