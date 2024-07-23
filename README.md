OpenQuad
========

Your open database for multi-dimensional numerical integration

TODO: add badges

---

OpenQuad offers a collection of highly efficient quadrature methods for
evaluating integrals on different domains and geometries, including 1d
intervals, spherical surfaces, Euler angles and more.  These methods serve as
building blocks for constructing custom quadratures, making this package
versatile for integrals over arbitrary multi-dimensional domains and tensor
spaces.

**Key features**:

- Access a large collection of state-of-the-art quadrature methods (Gauss
  quadratures, spherical designs, uniform coverings, ...).
- Combine methods to create customized product quadratures.
- Integrate Python functions and array data with a single command.
- Export quadrature points and weights for use in other software.

Found a bug? [Open an issue](https://github.com/qucontrol/openquad/issues).  
Missing a feature? [Open a pull request](https://github.com/qucontrol/openquad/pulls).

We appreciate and welcome your [contribution]()!

---
> [!NOTE]
> The first release will be published soon!  
> Give this repo a star to stay tuned!
---

Installation
------------

This package is available on
[PyPi](https://pypi.org/project/openquad). Install it with `pip`
into your active environment:

```bash
python -m pip install orientation_average
```

Usage
-----

Create a quadrature method for an integral over the three [Euler angles](),
using [Lebedev-Laikov quadrature]() of degree `5` for the first two angles
combined with the composite trapezoid rule with `6` sample points for the third
angle, and export sample points and weights in a text file:

```python
from orientation_average import SO3

quad = SO3([
    ('Lebedev-Laikov', degree=5),
    ('Trapezoid', size=6),
])
quad.savetxt('points_and_weights.dat')
```

Integrate a Python function `func(theta, phi)` over the surface of the
two-dimensional unit sphere using a [spherical design]() of degree `7`:

```python
from orientation_average import S2

quad = S2([
    ('S2-Design-Graef', degree=7),
])
quad.integrate(func)
```

Obtain Gauss-Legendre sample points and weights for degree `71` on the interval
`[-10, 5]`:

```python
from orientation_average import Rn

quad = Rn([
    ('Gauss-Legendre', degree=71, a=-10, b=5),
])
quad.points
quad.weights
```

For further information, including advanced examples, background information,
and details on the implementation, see the [documentation]().


Citation
--------

If this package was useful for your research, please [cite it]().
