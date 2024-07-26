import numpy as np
from scipy import special, integrate

import pytest
from contextlib import nullcontext as does_not_raise

import spherical
import quaternionic

from openquad.lebedev import LebedevLaikov
from openquad.graef import GraefS2Design, GraefS2Gauss
from openquad.womersley import WomersleyS2Design
from openquad.fibonacci import FibonacciSphere, ZCW2
from openquad.monte_carlo import MonteCarloS2


def f(theta, phi):
    """Test integrand. Integral should yield 1."""
    assert np.shape(theta) == np.shape(phi)
    poly = (
        1
        + 3 * special.sph_harm(3, 3, phi, theta)
        + 4 * special.sph_harm(2, 4, phi, theta)
        + 5 * special.sph_harm(4, 5, phi, theta)
    )
    return poly / (4*np.pi)


def test_testfuction():
    """Make sure that the integral over the test function yields 1."""
    assert integrate.dblquad(lambda t, p: np.sin(t)*f(t, p).real, 0, 2*np.pi, 0, np.pi)[0] == pytest.approx(1)
    assert integrate.dblquad(lambda t, p: np.sin(t)*f(t, p).imag, 0, 2*np.pi, 0, np.pi)[0] == pytest.approx(0)


@pytest.mark.parametrize(
    "method, size, degree, has_degree",
    [
        (LebedevLaikov, 6, 3, True),
        (GraefS2Gauss, 6, 3, True),
        (GraefS2Design, 6, 3, True),
        (WomersleyS2Design, 6, 3, True),
        (FibonacciSphere, 13, None, False),
        (ZCW2, 21, None, False),
    ]
)
def test_initialization(method, size, degree, has_degree):
    # missing required arguments
    with pytest.raises(TypeError):
        quad = method()

    # duplicate arguments
    if has_degree:
        with pytest.raises(TypeError, match="not both"):
            quad = method(degree=degree, size=size)

    # unexpected arguments
    with pytest.raises(TypeError, match="unexpected keyword"):
        quad = method(size=size, unexpected=True)

    # invalid number of points or accuracy
    with pytest.raises(ValueError, match="larger than zero"):
        quad = method(size=0)
    if has_degree:
        with pytest.raises(ValueError, match="not available"):
            quad = method(size=1)
        with pytest.raises(ValueError, match="not available"):
            quad = method(degree=0)

    # valid initializations
    with does_not_raise():
        quad = method(size=size)
    if has_degree:
        with does_not_raise():
            quad = method(degree=degree)


@pytest.mark.parametrize(
    "method, size, degree, kwargs",
    [
        (LebedevLaikov, 14, 5, {}),
        (GraefS2Gauss, 12, 5, {}),
        (GraefS2Design, 12, 5, {}),
        (WomersleyS2Design, 12, 5, {'symmetric':True}),
        (WomersleyS2Design, 18, 5, {'symmetric':False}),
        (FibonacciSphere, 21, None, {}),
        (ZCW2, 21, None, {}),
    ]
)
def test_fields_and_integration(method, size, degree, kwargs):
    tol = 1e-6
    quad = method(size=size, **kwargs)
    
    # test fields
    assert quad.dim == 2
    assert quad.degree == degree
    assert quad.size == size

    # test weights
    assert quad.weights.shape == (size,)
    assert np.sum(quad.weights) == pytest.approx(4*np.pi)

    # test grid dimension
    assert quad.points.shape == (2, size)
    assert np.all(quad.angles == quad.points)

    # test for callable
    assert quad.integrate(f) == pytest.approx(1, rel=tol)
    # test for samples
    theta, phi = quad.points
    f_samples = f(theta, phi)
    assert quad.integrate(f_samples) == pytest.approx(1, rel=tol)



@pytest.mark.parametrize(
    "method, tol, kwargs",
    [
        (LebedevLaikov, 1e-11, {}),
        (GraefS2Gauss, 1e-9, {}),
        (GraefS2Design, 1e-8, {}),
        (WomersleyS2Design, 1e-12, {}),
    ]
)
def test_poly_acc(method, tol, kwargs):
    max_degree = 50
    wigner = spherical.Wigner(max_degree)
    available_degree = method._available_degrees
    tested_degree = available_degree[available_degree <= max_degree]
    for degree in tested_degree:
        quad = method(degree=degree, **kwargs)
        # First, test that the weights are normalized:
        assert np.sum(quad.weights) == pytest.approx(4*np.pi)
        # Then, test the integration of all polynomials:
        q = quaternionic.array.from_spherical_coordinates(*quad.angles)
        Ylm = wigner.sYlm(0, q)[:, :wigner.Yindex(degree, degree)+1]
        results = quad.integrate(Ylm, axis=0)
        truevals = np.concatenate(
            ([np.sqrt(4*np.pi)], np.zeros(Ylm.shape[-1]-1))
        )
        assert np.abs(results - truevals).max() == pytest.approx(0, abs=tol)
