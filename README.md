<h1 align="center">
<img src="docs/_static/images/logo.svg" width="300">
<br>
Your open database for multi-dimensional numerical integration
</h1><br>

<!-- TODO: add badges -->

---
> [!NOTE]
> The first stable release is under active development.  
> Watch this repo to stay tuned for any updates!
---

<!-- start including on doc landing page -->
OpenQuad offers a collection of highly efficient quadrature methods for
evaluating integrals on different domains and geometries, including 1d
intervals, spherical surfaces, Euler angles and more. These methods serve as
building blocks for constructing custom quadratures, making this package
versatile for integrals over arbitrary multi-dimensional domains and tensor
spaces.
<!-- end including on doc landing page -->

**Key features**:

- Access a large collection of state-of-the-art quadrature methods (Gauss
  quadratures, spherical designs, uniform coverings, ...).
- Combine methods to create customized product quadratures.
- Integrate Python functions and array data with a single command.
- Export quadrature points and weights for use in other software.

Found a bug? [Open an issue](https://github.com/qucontrol/openquad/issues).  
Missing a feature? [Start a discussion](https://github.com/qucontrol/openquad/discussions).

We appreciate and welcome your [contribution][contribute]!

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
python -m pip install openquad
```

Usage
-----

Create a quadrature method for an integral over the three [Euler angles][angles],
using [Lebedev-Laikov quadrature][lebedev] of degree `5` for the first two angles
combined with the composite trapezoid rule with `6` sample points for the third
angle, and export sample points and weights in a text file:

```python
from openquad import SO3

quad = SO3([
    ('LebedevLaikov', dict(degree=5)),
    ('Trapezoid', dict(size=6)),
])
quad.savetxt('points_and_weights.dat')
```

Integrate a Python function `func(theta, phi)` over the surface of the
two-dimensional unit sphere using a [spherical design][designs] of degree `7`:

```python
from openquad import S2

quad = S2([
    ('S2-Design-Graef', dict(degree=7)),
])
quad.integrate(func)
```

Obtain Gauss-Legendre sample points and weights for degree `71` on the interval
`[-10, 5]`:

```python
from openquad import Rn

quad = Rn([
    ('GaussLegendre', dict(degree=71, a=-10, b=5)),
])
quad.points
quad.weights
```

For further information, including advanced examples, background information,
and details on the implementation, see the [documentation][docs].


Citation
--------

If this package was useful for your research, please [cite it][cite].


[angles]: https://en.wikipedia.org/wiki/Euler_angles
[designs]: https://en.wikipedia.org/wiki/Spherical_design
[lebedev]: https://en.wikipedia.org/wiki/Lebedev_quadrature
[docs]: https://qucontrol.github.io/openquad
[cite]: https://qucontrol.github.io/openquad/using/cite.html
[contribute]: https://qucontrol/github.io/openquad/contributing/index.html
