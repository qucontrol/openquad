.. _implemented-methods:

Quadrature methods
------------------

Currently, the following quadrature methods are implemented, sorted according
to their geometry.

The first entry denotes the input string for the method specification, cf. the
``method_spec`` parameter in :class:`openquad.Rn`.


:math:`\mathbb{R}^1` methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``GaussLegendre``: `Gauss-Legendre quadrature <https://en.wikipedia.org/wiki/Gauss-Legendre_quadrature>`_
* ``GaussLegendreLobatto``: `Gauss-Lobatto-Legendre quadrature <https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss-Lobatto_rules>`_
* ``Trapezoid``: `Composite trapezoidal rule <https://en.wikipedia.org/wiki/Trapezoidal_rule>`_, see :func:`scipy.integrate.trapezoid`
* ``Simpson``: `Composite Simpson's rule <https://en.wikipedia.org/wiki/Simpson's_rule>`_, see :func:`scipy.integrate.simpson`
* ``Romberg``: `Romberg's method <https://en.wikipedia.org/wiki/Romberg's_method>`_, see :func:`scipy.integrate.romberg`
* ``1d-MonteCarlo``:`Monte Carlo integration <https://en.wikipedia.org/wiki/Monte_Carlo_integration>`_

These methods need the following additional keyword paramters:

* ``size``: float: Number of sampling points.
* ``degree``: float: Degree of exactness. For Gauss methods alternative to ``size``.
* ``a`` : float: Lower boundary of the integral.
* ``b`` : float: Upper boundary of the integral.
* ``jacobian`` : callable, optional: Jacobian to apply to the coordinate.
* ``periodic`` : logical, optinal: If ``True``, assume periodic boundary conditions.


:math:`\mathrm{S}^2` methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``S2-Gauss-LebedevLaikov``: spherical Gauss quadrature from `Lebedev-Laikov`_
* ``S2-Gauss-Graef``: spherical Gauss qudarature from `Manuel Gräf`_
* ``S2-Design-Graef``: spherical designs from `Manuel Gräf`_
* ``S2-Design-Womersley``: spherical designs from `Robert Womersley`_
* ``S2-Covering-Fibonacci``: near-uniform coverings based on `Fibonacci numbers <https://en.wikipedia.org/wiki/Fibonacci_sequence>`_
* ``S2-Covering-ZCW``: near-uniform coverings with the `ZCW`_ method
* ``S2-MonteCarlo``: `Monte Carlo integration <https://en.wikipedia.org/wiki/Monte_Carlo_integration>`_

These methods need the following additional keyword paramters:

* ``size``: float: Number of sampling points.
* ``degree``: float: Degree of exactness. For Gauss methods and spherical desings alternative to ``size``.


:math:`\mathrm{SO}(3)` methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``SO3-Gauss-Graef``: spherical Gauss qudarature from `Manuel Gräf`_
* ``SO3-Chebyshev-Graef``: spherical Chebyshev quadrature from `Manuel Gräf`_
* ``SO3-Covering-Karney``: near-uniform coverings from `Charles Karney`_
* ``SO3-MonteCarlo``: `Monte Carlo integration <https://en.wikipedia.org/wiki/Monte_Carlo_integration>`_

These methods need the following additional keyword paramters:

* ``size``: float: Number of sampling points.
* ``degree``: float: Degree of exactness. For Gauss and Chebyshev methods alternative to ``size``.


.. note::

   The following extensions are planned for the next release:
   
   * Implementing the three-dimensional unit sphere, :math:`\mathrm{S}^3`.
   * Adding spherical desings and coverings from `Neil Sloane and Ronald Hardin`_.
   * Providing a way to use :math:`\mathrm{S}^3` quadratures for :math:`\mathrm{SO}(3)`
     and vice versa.
   * Implementing a class for the unit circle, :math:`\mathrm{S}^1`.
   * and more...


.. _Manuel Gräf: https://www-user.tu-chemnitz.de/~potts/workgroup/graef/quadrature/
.. _Robert Womersley: https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/
.. _Neil Sloane and Ronald Hardin: http://www.neilsloane.com/
.. _Charles Karney: https://github.com/cffk/orientation
.. _Lebedev-Laikov: https://doi.org/10.1016/0041-5553(75)90133-0
.. _ZCW: https://doi.org/10.1006/jmre.1998.1427
