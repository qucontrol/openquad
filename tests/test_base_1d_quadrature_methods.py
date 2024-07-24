import numpy as np
from scipy import integrate, special

import pytest
from contextlib import nullcontext as does_not_raise

from openquad.newton_cotes import CompositeTrapezoid, CompositeSimpson, Romberg
from openquad.gauss import GaussLegendre, GaussLobattoLegendre

#@pytest.fixture()
def f(x):
    """Test integrand. Integral should yield 1."""
    return (5*x + 4*x**2 + 3*x**3 + 2*x**4 + 1*x**5) / 215750


def test_test_function():
    """Make sure that the integral over the test function yields 1."""
    a = 0
    b = 10
    assert integrate.quad(f, a, b)[0] == pytest.approx(1)


@pytest.mark.parametrize(
    "method, has_poly_acc",
    [
        (CompositeTrapezoid, False),
        (CompositeSimpson, False),
        (Romberg, False),
        (GaussLegendre, True),
        (GaussLobattoLegendre, True),
    ]
)
def test_initialization(method, has_poly_acc):
    # missing required arguments
    with pytest.raises(TypeError):
        quad = method()
    with pytest.raises(TypeError):
        quad = method(size=5)
    with pytest.raises(TypeError):
        quad = method(size=5, a=0)
    with pytest.raises(TypeError):
        quad = method(size=5, b=1)
    with pytest.raises(TypeError):
        quad = method(a=0, b=1)
    
    # duplicate arguments
    if has_poly_acc:
        with pytest.raises(TypeError, match="not both"):
            quad = method(degree=5, size=5, a=0, b=1)

    # unexpected arguments
    with pytest.raises(TypeError, match="unexpected keyword"):
        quad = method(size=5, a=0, b=1, unexpected=True)

    # invalid integration interval
    with pytest.raises(ValueError):
        quad = method(size=5, a=0, b=0)

    # invalid jacobian
    with pytest.raises(TypeError):
        quad = method(size=5, a=0, b=1, jacobian=np.ones(5))

    # invalid number of points
    with pytest.raises(ValueError, match="larger than zero"):
        quad = method(size=0, a=0, b=1)

    # valid initializations
    with does_not_raise():
        quad = method(size=5, a=0, b=1)
        quad = method(size=5, a=0, b=-1)
        quad = method(size=5, a=0, b=1, jacobian=np.sin)
        quad = method(size=5, a=0, b=1, jacobian=None)
    if has_poly_acc:
        with does_not_raise():
            quad = method(degree=5, a=0, b=1)


@pytest.mark.parametrize(
    "method, size, has_normalized_weights, poly_acc",
    [
        (CompositeTrapezoid, 50, True, lambda size: 1),
        (CompositeSimpson, 21, True, lambda size: 3),
        (CompositeSimpson, 22, True, lambda size: 3),
        (Romberg, 9, False, lambda size: None),
        (GaussLegendre, 5, True, lambda size: 2*size - 1),
        (GaussLobattoLegendre, 5, True, lambda size: 2*size - 3),
    ]
)
def test_fields_and_integration(method, size, has_normalized_weights, poly_acc):
    a = 0
    b = 10
    tol = 1e-3
    quad = method(size=size, a=a, b=b)

    # test fields
    assert quad.dim == 1
    assert quad.size == size
    assert quad.degree == poly_acc(size)

    # test weights
    if has_normalized_weights:
        assert np.sum(quad.weights) == pytest.approx(b-a)
    else:
        assert np.all(quad.weights == np.ones(size))

    # test grid dimension
    assert quad.points.shape == (size,)

    # test for callable
    assert quad.integrate(f) == pytest.approx(1, rel=tol)
    # test for samples
    x = quad.points
    f_samples = f(x)
    assert quad.integrate(f_samples) == pytest.approx(1, rel=tol)


@pytest.mark.parametrize(
    "method, tol",
    [
        (GaussLegendre, 1e-13),
        (GaussLobattoLegendre, 1e-13),
    ]
)
def test_poly_acc(method, tol):
    """Test the polynomial accuracy of the given method by checking that all
    polynomials up to the required degree are integrated exactly.

    Here, Legendre polynomials integrated from -1 to 1.
    """
    max_degree  = 100
    for degree in np.arange(max_degree +1):
        quad = method(degree=degree, a=-1, b=1)
        # First, test that the weights are normalized:
        assert np.sum(quad.weights) == pytest.approx(2)
        # Then, test the integration of all polynomials:
        orders = np.arange(degree+1)[:, np.newaxis]
        Pls = special.eval_legendre(orders, quad.points)
        results = quad.integrate(Pls, axis=-1)
        truevals = np.concatenate(([2], np.zeros(Pls.shape[0]-1)))
        assert np.abs(results - truevals).max() == pytest.approx(0, abs=tol)

def none_type_quadrature():
    points = np.linspace(0, 1, 10)
    quad = NoQuadrature(points=points)
    assert quad.size == points.size
    assert np.all(quad.points == points)

    func = lambda x: x
    assert np.all(quad.integrate(func) == func(points))
