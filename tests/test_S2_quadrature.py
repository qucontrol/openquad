import numpy as np
from scipy import special, integrate

import pytest
from contextlib import nullcontext as does_not_raise

from openquad import S2


def f(theta, phi):
    """Test integrand. Integral should yield 1."""
    assert np.shape(theta) == np.shape(phi)
    poly = (
        1
        + 3 * special.sph_harm(0, 1, phi, theta)
        + 3 * special.sph_harm(3, 3, phi, theta)
        + 4 * special.sph_harm(2, 4, phi, theta)
        + 5 * special.sph_harm(4, 5, phi, theta)
    )
    return poly / (4*np.pi)


def test_test_function():
    """Make sure that the integral over the test function yields 1."""
    assert integrate.dblquad(lambda t, p: np.sin(t)*f(t, p).real, 0, 2*np.pi, 0, np.pi)[0] == pytest.approx(1)
    assert integrate.dblquad(lambda t, p: np.sin(t)*f(t, p).imag, 0, 2*np.pi, 0, np.pi)[0] == pytest.approx(0)


## see answer in this post:
## https://stackoverflow.com/a/68012715
#@pytest.mark.parametrize(
#    "methods_spec, expectation, match",
#    [
#        ([('trapezoid',)], pytest.raises(ValueError), "invalid choice of methods"),
#        ([('gausslegendre',)], pytest.raises(ValueError), "invalid choice of methods"),
#        ([('gausslegendre',), ('gausslegendre',), ('gausslegendre',)], pytest.raises(ValueError), "invalid choice of methods"),
#        ([('gausslegendre',), ('s2-gauss-lebedevlaikov',)], pytest.raises(ValueError), "invalid choice of methods"),
#        ([('ZCW',)], pytest.raises(ValueError), "invalid choice of methods"),
#        
def test_initialization():
    # Wrong initializations
    # Before checking the method specification for consistency,
    # the methods shall raise an error, if a wrong type is given.
    with pytest.raises(TypeError): # missing method specifier
        quad = S2()
    with pytest.raises(ValueError, match="empty"): # empty method specifier
        quad = S2([])

    # Invalid methods
    # Given the type of the specifier is correct, first the validity of the
    # method shall be checked before checking the remaining options.
    with pytest.raises(ValueError, match='Invalid'):
        quad = S2([('invalid', {})])
    with pytest.raises(ValueError, match="Invalid"):
        quad = S2([('trapezoid', {})])
    with pytest.raises(ValueError, match="Invalid"):
        quad = S2([('gausslegendre', {})])
    with pytest.raises(ValueError, match="Invalid"):
        quad = S2([
            ('gausslegendre', {}),
            ('gausslegendre', {}),
            ('gausslegendre', {})])
    with pytest.raises(ValueError, match="Invalid"):
        quad = S2([('gausslegendre', {}), ('lebedevlaikov', {})])
    with pytest.raises(ValueError, match="Invalid"):
        quad = S2([('so3-montecarlo', {})])

    # Invalid method options
    # Now given that the type of the specifier is correct and that a valid
    # combination of methods has been choosen, check the remaining options
    # for the method.
    with pytest.raises(TypeError, match='not both'): # either degree or size
        quad = S2([
            ('s2-gauss-lebedevlaikov', {'size':6, 'degree':3}),
        ])
    with pytest.raises(TypeError, match='unexpected'): # unexpected argument
        quad = S2([
            ('s2-gauss-lebedevlaikov', {'size':6, 'unexpected':True}),
        ])

    # Valid initializations
    with does_not_raise():
        quad = S2([
            ('gAuSslObAtTolEgEnDrE', {'size':5}),
            ('trapezoid', {'size':5}),
        ])
        quad = S2([
            ('simpson', {'size':5}),
            ('Trapezoid', {'size':5}),
        ])
        quad = S2([('s2-gauss-lebedevlaikov', {'size':6})])
        quad = S2([('s2-gauss-lebedevlaikov', {'degree':3})])
        quad = S2([
            ('none', {'points':[0]}),
            ('trapezoid', {'size':5}),
        ])


def test_fields_product_method():
    """Test the various fields of the quadrature class.
    Here we test a product method.
    """
    # Test a product quadrature:
    n1 = 5
    n2 = 6
    x1 = np.linspace(-1, 1, n1)
    quad = S2([
        ('trapezoid', {'size':n1}),
        ('gausslegendre', {'size':n2}),
    ])
    # test submethod properties:
    assert quad._methods[0].degree == 1
    assert quad._methods[1].degree == 2*n2-1
    assert quad._methods[0].size == n1
    assert quad._methods[1].size == n2
    assert len(quad._method_weights) == 2
    assert len(quad._method_points) == 2
    assert quad._method_points[0].shape == (1, n1)
    assert quad._method_points[1].shape == (1, n2)
    assert np.all(quad._method_points[0] == x1)

    # test hidden properties
    assert quad._ndims == [1, 1]
    assert quad._dims == [(0,), (1,)] # do I really need this?
    assert quad._sizes == [n1, n2]
    assert quad._degrees == [1, 2*n2-1]
    assert quad._points.size == (2*n1*n2)
    assert quad._points.shape == (2, n1*n2)

    # test public properties:
    assert quad.dim == 2
    assert quad.size == n1*n2
    assert quad.shape == (n1, n2)
    assert quad.weights.shape == (n1*n2,)
    assert quad.angles.shape == (2, n1*n2)
    assert quad.xyz.shape == (3, n1*n2)
    assert np.all(quad.angles == quad._points)
    #assert quad.grid.angles.shape == (2, n1*n2)
    #assert quad.grid.xyz.shape == (3, n1*n2)
    assert quad.weights.reshape(quad.shape).shape == (n1, n2)
    assert quad.angles.reshape(2, *quad.shape).shape == (2, n1, n2)
    assert quad.xyz.reshape(3, *quad.shape).shape == (3, n1, n2)


def test_fields_non_product_method():
    """Test the various fields of the quadrature class.
    Here we test a non-product method.
    """
    degree = 3 # degree of scheme
    n = 6 # number of sample points
    quad = S2([
        ('s2-gauss-lebedevlaikov', {'degree':degree}),
    ])
    # test submethod properties:
    assert quad._methods[0].degree == degree
    assert quad._methods[0].size == n
    assert len(quad._method_weights) == 1
    assert len(quad._method_points) == 1
    assert quad._method_points[0].shape == (2, n)

    # test hidden properties:
    assert quad._ndims == [2]
    assert quad._dims == [(0, 1)] # do I really need this?
    assert quad._sizes == [n]
    assert quad._degrees == [degree]
    assert quad._points.shape == (2, n)

    # test public properties:
    assert quad.dim == 2
    assert quad.size == n
    assert quad.shape == (n,)
    assert quad.weights.shape == (n,)
    assert quad.angles.shape == (2, n)
    assert quad.xyz.shape == (3, n)
    print(quad.angles)
    assert np.all(quad.angles == quad._points)
    #assert quad.grid.angles.shape == (2, n)
    #assert quad.grid.xyz.shape == (3, n)
    assert quad.weights.reshape(quad.shape).shape == (n,)
    assert quad.angles.reshape(2, *quad.shape).shape == (2, n)
    assert quad.xyz.reshape(3, *quad.shape).shape == (3, n)
    


def test_integration():
    quad = S2([
        ('gausslegendre', dict(size=6)),
        ('trapezoid', dict(size=5)),
    ], polar_sampling='cos')
    assert quad.integrate(f) == pytest.approx(1)
    f_samples = f(*quad.angles)
    assert quad.integrate(f_samples) == pytest.approx(1)

    quad = S2([
        ('gausslegendre', dict(size=6)),
        ('trapezoid', dict(size=5)),
    ], polar_sampling='angle')
    assert quad.integrate(f) == pytest.approx(1)
    f_samples = f(*quad.angles)
    assert quad.integrate(f_samples) == pytest.approx(1)

    quad = S2([
        ('s2-gauss-lebedevlaikov', {'degree':9}),
    ])
    assert quad.integrate(f) == pytest.approx(1)
    f_samples = f(*quad.angles)
    assert quad.integrate(f_samples) == pytest.approx(1)
