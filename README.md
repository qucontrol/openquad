<p align="center">
<img alt="openquad" src="https://raw.githubusercontent.com/qucontrol/openquad/main/doc/_static/images/logo.svg" width="300">
<br>
Open database for multi-dimensional numerical integration
</p>

[![Source code on Github](https://img.shields.io/badge/github-qucontrol/openquad-blue.svg)](https://github.com/qucontrol/openquad)
[![Documentation](https://img.shields.io/badge/docs-gh--pages-blue.svg)](https://qucontrol.github.io/openquad)
[![Openquad on the Python Package Index](https://img.shields.io/pypi/v/openquad.svg)](https://pypi.python.org/pypi/openquad)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![arXiv](https://img.shields.io/badge/arXiv-2407.17434-b31b1b.svg)](https://arxiv.org/abs/2407.17434)

---
> [!NOTE]
> The first stable release is under active development.  
> Stay tuned for any updates!
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

[**Get started!**][docs]

Found a bug? [Open an issue](https://github.com/qucontrol/openquad/issues).  
Missing a feature? [Start a discussion](https://github.com/qucontrol/openquad/discussions).

We appreciate and welcome your [contribution][contribute]!


Installation
------------

This package is available on
[PyPi](https://pypi.org/project/openquad). Install it with `pip`
into your active environment:

```bash
python -m pip install openquad
```

Basic usage
-----------

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

Integrate a function `func(theta, phi)` over the surface of the
two-dimensional unit sphere using a [spherical design][designs] of degree `7`:

```python
from openquad import S2

quad = S2([
    ('S2-Design-Graef', dict(degree=7)),
])
quad.integrate(func)
```

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

For further information, including advanced examples, background information,
and details on the implementation, see the [documentation][docs].


Citation
--------

If this package was useful for your research, please [cite it][cite].


[angles]: https://en.wikipedia.org/wiki/Euler_angles
[designs]: https://en.wikipedia.org/wiki/Spherical_design
[lebedev]: https://en.wikipedia.org/wiki/Lebedev_quadrature
[docs]: https://qucontrol.github.io/openquad
[cite]: https://qucontrol.github.io/openquad/latest/using/cite.html
[contribute]: https://qucontrol/github.io/openquad/latest/contributing/index.html
