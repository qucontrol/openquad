.. _implemented-methods:

Quadrature methods
------------------

Currently, the following quadrature methods are implemented, sorted according
to their geometry.


:math:`\mathbb{R}^1` methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `Gauss-Legendre quadrature <https://en.wikipedia.org/wiki/Gauss-Legendre_quadrature>`_
* `Gauss-Lobatto-Legendre quadrature <https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss-Lobatto_rules>`_
* `Composite trapezoidal rule <https://en.wikipedia.org/wiki/Trapezoidal_rule>`_, see :scipy:`trapezoid`
* `Composite Simpson's rule <https://en.wikipedia.org/wiki/Simpson's_rule>`_, see :scipy:`simpson`
* `Romberg's method <https://en.wikipedia.org/wiki/Romberg's_method>`_, see :scipy:`romberg`
* `Monte Carlo integration <https://en.wikipedia.org/wiki/Monte_Carlo_integration>`_


:math:`\mathrm{S}^2` methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* spherical Gauss quadrature from `Lebedev-Laikov`_
* spherical Gauss qudarature from `Manuel Gräf`_
* spherical designs from `Manuel Gräf`_
* spherical designs from `Robert Womersley`_
* near-uniform coverings based on `Fibonacci numbers <https://en.wikipedia.org/wiki/Fibonacci_sequence>`_
* near-uniform coverings with the `ZCW`_ method
* `Monte Carlo integration <https://en.wikipedia.org/wiki/Monte_Carlo_integration>`_


:math:`\mathrm{SO}(3)` methods
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* spherical Gauss qudarature from `Manuel Gräf`_
* spherical Chebyshev quadrature from `Manuel Gräf`_
* near-uniform coverings from `Charles Karney`_
* `Monte Carlo integration <https://en.wikipedia.org/wiki/Monte_Carlo_integration>`_


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
