import numpy as np
from scipy import special, integrate

import pytest
from contextlib import nullcontext as does_not_raise

from openquad import Rn

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


def test_test_function():
    """Make sure that the integral over the test function yields 1."""
    assert integrate.dblquad(lambda t, p: np.sin(t)*f(t, p).real, 0, 2*np.pi, 0, np.pi)[0] == pytest.approx(1)
    assert integrate.dblquad(lambda t, p: np.sin(t)*f(t, p).imag, 0, 2*np.pi, 0, np.pi)[0] == pytest.approx(0)


def test_initialization():
    # Wrong initializations
    # Before checking the method specification for consistency,
    # the methods shall raise an error, if a wrong type is given.
    with pytest.raises(TypeError): # missing method specifier
        quad = Rn()
    with pytest.raises(ValueError, match="empty"): # empty method specifier
        quad = Rn([])
    # unsopproted types for method specifier
    # (we allow for ducktyping)
    #with pytest.raises(TypeError): #TODO remove test
    #    quad = Rn(())
    #with pytest.raises(TypeError):
    #    quad = Rn('str')
    #with pytest.raises(TypeError): #TODO: remove?
    #    quad = Rn([('str', [])])
    #with pytest.raises(TypeError): #TODO: remove?
    #    quad = Rn([('str', ())])

    # Invalid methods
    # Given the type of the specifier is correct, first the validity of the
    # method shall be checked before checking the remaining options.
    with pytest.raises(ValueError, match='Invalid'):
        quad = Rn([('invalid', {})])
    with pytest.raises(ValueError, match='Invalid'):
        quad = Rn([('s2-gauss-lebedevlaikov', {})]) # example of 2d method
    with pytest.raises(ValueError, match='Invalid'):
        quad = Rn([('so3-montecarlo', {})]) # example of 3d method
    
    # Invalid method options
    # Now given that the type of the specifier is correct and that a valid
    # combination of methods has been choosen, check the remaining options
    # for the method.
    with pytest.raises(TypeError, match='not both'): # either degree or size
        quad = Rn([
            ('gausslegendre', {'size':5, 'degree':4, 'a':0, 'b':1}),
        ])
    with pytest.raises(TypeError, match='`size` or `degree` must be given'):
        quad = Rn([
            ('gausslegendre', {'a':0, 'b':np.pi}),
        ])
    with pytest.raises(ValueError, match='`a`'): # no interval specified
        quad = Rn([
            ('trapezoid', {'size':33, 'b':np.pi, 'jacobian':np.sin}),
        ])
    with pytest.raises(TypeError, match='unexpected'): # unexpected argument
        quad = Rn([
            ('trapezoid', {'size':5, 'a':0, 'b':1, 'unexpected':True}),
        ])

    # Test optional parameters
    with pytest.raises(TypeError, match='jacobian'):
        quad = Rn([
            ('trapezoid', {'size':33, 'a':0, 'b':np.pi, 'jacobian':np.ones(11)}),
        ])

    # Valid initializations
    with does_not_raise():
        quad = Rn([('trapezoid', {'size':33, 'a':0, 'b':np.pi})])
        quad = Rn([('tRaPeZoId', {'size':33, 'a':0, 'b':np.pi})])
        quad = Rn([('trapezoid', {'size':33, 'a':0, 'b':np.pi})])
        quad = Rn([
            ('trapezoid', {'size':33, 'a':0, 'b':np.pi, 'jacobian':np.sin}),
        ])


def test_only_1d_methods_allowed():
    with pytest.raises(ValueError, match='Invalid'):
        quad = Rn([('s2-gauss-lebedevlaikov', {})]) # example of 2d method
    with pytest.raises(ValueError, match='Invalid'):
        quad = Rn([('so3-montecarlo', {})]) # example of 3d method


def test_fields():
    """Test the various fields of the quadrature class."""
    n1 = 5
    n2 = 6
    x1 = np.linspace(0, 1, n1)
    quad = Rn([
        ('trapezoid', {'size':n1, 'a':0, 'b':1}),
        ('gausslegendre', {'size':n2, 'a':0, 'b':2}),
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

    # test hidden properties:
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
    assert quad.points.shape == (2, n1*n2)
    assert np.all(quad.points == quad._points)
    assert np.all(quad.r == quad.points) # TODO: remove?
    #assert quad.grid.points.shape == (2, n1*n2)
    assert quad.weights.reshape(quad.shape).shape == (n1, n2)
    assert quad.points.reshape(quad.dim, *quad.shape).shape == (2, n1, n2)


def test_integration():
    quad = Rn([
        ('gausslegendre', {'size':11,  'a':0, 'b':np.pi, 'jacobian':np.sin}),
        ('simpson', {'size':13, 'a':0, 'b':2*np.pi}),
    ])
    # Complete integration
    # callable
    assert quad.integrate(f) == pytest.approx(1)
    # flat array
    f_samples = f(*quad.points)
    assert quad.integrate(f_samples) == pytest.approx(1)
    # TODO: shaped array?

    #TODO: below?
    ## Partial integration (of a sample)
    #f_samples = f(*quad.points)
    #print(quad.integrate(f_samples, dim=0))
    #print(quad.integrate(f_samples, dim=(0, 1)))
