.. _quickstart:

Quickstart guide
----------------


What is this package?
^^^^^^^^^^^^^^^^^^^^^

Combining a comprehenseive collection of quadrature methods in one single
interface, :mod:`orientation_average` helps you to efficiently compute
integrals over various domains with ease, including 1d intervals, spherical
surfaces, Euler angles and more. The object-oriented interface provides Python
classes for each geometry, for example the :class:`S2` class for integration
over the 2d unit sphere. Python functions and array data can be integrated
directly using the :method:`integrate` method. Alternatively, sample points and
quadrature weights can be exported with the :method:`savetxt` method for easy
integration with other software.

The best part: combine individual quadrature methods to create custom
multi-dimensional quadratures.


How do I use it?
^^^^^^^^^^^^^^^^

.. todo: testsetup into conftest.py
.. testsetup::

   >>> import numpy as np

Want to integrate a Python function over the interval :math:`[-10,5]` with
Gauss-Legendre quadrature?

.. testcode::

    >>> from orientation_average import Rn
    >>> quad = Rn([
    >>>     ('GaussLegendre', degree=71, a=-10, b=5),
    >>> ])
    >>> quad.integrate(your_awesome_function)

Have a Python function that you need to evaluate on a uniform grid over the surface of the unit sphere?

.. testcode::
   :hide:

   >>> def your_awesome_function(theta, phi):
   >>>  return np.sin(theta) * np.cos(phi)

.. testcode::

    >>> from orientation_average import S2
    >>> quad = S2([
    >>>     ('S2-Covering', size=142),
    >>> ])
    >>> grid = your_awesome_function(*quad.angles)

Need quadrature points and weights for performing an efficient orientation average in your molecular dynamics simulation?

.. testcode::

    >>> from orientation_average import SO3
    >>> quad = SO3([
    >>>     ('LebedevLaikov', degree=35),
    >>>     ('Trapezoid', size=36),
    >>> ])
    >>> quad.savetxt('points_and_weights.dat')


How do I get it?
^^^^^^^^^^^^^^^^

Install the latest stable release with

.. code-block:: bash

    pip install orientation_average


Tell me more!
^^^^^^^^^^^^^

Have a look at the :ref:`user guide <using>` and the :ref:`example gallery <examples>`.

Don't know, which method to choose? Read our `paper`_.

.. _paper: https://arxiv.org/abs/2407.17434
