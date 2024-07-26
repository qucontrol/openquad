import numpy as np
from scipy import integrate

import spherical
import quaternionic

import pytest
from contextlib import nullcontext as does_not_raise

from openquad import SO3


#@pytest.fixture()
def f_q(quaternions):
    """Test integrand for quaternions.
    Integral should yield 1.
    """
    l_max = 2
    wigner = spherical.Wigner(l_max)
    poly = (
        1
        + 1 * wigner.D(quaternions)[..., wigner.Dindex(2, 0, 0)]
        + 1 * wigner.D(quaternions)[..., wigner.Dindex(2, 1, -2)]
    )
    if np.ndim(quaternions) > 1:
        assert poly.shape == (len(quaternions),)
    else:
        assert poly.shape == ()
    return poly / (8*np.pi**2)


#@pytest.fixture()
def f_abg(alpha, beta, gamma):
    """Test integrand for Euler angles.
    Integral should yield 1.
    """
    assert np.shape(alpha) == np.shape(beta) == np.shape(gamma)
    quaternions = quaternionic.array.from_euler_angles(alpha, beta, gamma)
    return f_q(quaternions)


#@pytest.fixture()
def f_vector_q(quaternions):
    """Test vector integrand, containing the orientations in the first axis.
    Integral should yield 1.
    """
    l_max = 1
    wigner = spherical.Wigner(l_max)
    idx_min = wigner.Dindex(1, -1, -1)
    idx_max = wigner.Dindex(1+1, -1+1, -2+1)
    poly = (
        1
        + 1 * wigner.D(quaternions)[..., idx_min:idx_max].reshape(-1, 3, 3)
    )
    assert poly.shape == (len(quaternions), 3, 3)
    return poly / (8*np.pi**2)


def test_test_function():
    """Make sure that the integral over the test function yields 1."""
    assert integrate.tplquad(lambda a, b, g: np.sin(b)*f_abg(a, b, g).real, 0, 2*np.pi, 0, np.pi, 0, 2*np.pi)[0] == pytest.approx(1)
    #assert integrate.tplquad(lambda a, b, g: np.sin(b)*f_abg(a, b, g).imag, 0, 2*np.pi, 0, np.pi, 0, 2*np.pi, epsabs=1e-1)[0] == pytest.approx(0)


def test_initialization():
    # Wrong initializations
    # Before checking the method specification for consistency,
    # the methods shall raise an error, if a wrong type is given.
    with pytest.raises(TypeError): # missing method specifier
        quad = SO3()
    with pytest.raises(ValueError, match="empty"): # empty method specifier
        quad = SO3([])

    # Invalid methods
    # Given the type of the specifier is correct, first the validity of the
    # method shall be checked before checking the remaining options.
    with pytest.raises(ValueError, match="Invalid"):
        quad = SO3([('trapezoid', {})])
    with pytest.raises(ValueError, match="Invalid"):
        quad = SO3([('s2-gauss-lebedevlaikov', {})])
    with pytest.raises(ValueError, match="Invalid"):
        quad = SO3([
            ('simpson', {}),
            ('simpson', {}),
            ('simpson', {}),
            ('simpson', {}),
         ])
    with pytest.raises(ValueError, match="Invalid"):
        quad = SO3([
            ('s2-gauss-lebedevlaikov', {}), 
            ('gausslegendre', {}), 
            ('gausslegendre', {})
        ])
    with pytest.raises(ValueError, match="Invalid"):
        quad = SO3([('s2-gauss-lebedevlaikov', {}), ('s2-gauss-lebedevlaikov', {})])

    # Invalid method options
    # Now given that the type of the specifier is correct and that a valid
    # combination of methods has been choosen, check the remaining options
    # for the method.
    with pytest.raises(TypeError, match='not both'): # either degree or size
        quad = SO3([
            ('s2-gauss-lebedevlaikov', {'size':6, 'degree':3}),
            ('gausslegendre', {'size':6, 'degree':11}),
        ])
    with pytest.raises(TypeError, match='unexpected'): # unexpected argument
        quad = SO3([
            ('so3-montecarlo', {'size':6, 'unexpected':True}),
        ])

    # Valid initializations
    with does_not_raise():
        quad = SO3([
            ('s2-gauss-LebeDevlaikov', {'degree':5}),
            ('Trapezoid', {'size':5}),
        ])
        quad = SO3([
            ('simpson', {'size':5}),
            ('Trapezoid', {'size':5}),
            ('Gausslegendre', {'size':5}),
        ])
        quad = SO3([('so3-montecarlo', {'size':5})])
        quad = SO3([('so3-montecarlo', {'size':5, 'seed':1})])


def test_fields():
    """Test the various fields of the quadrature class."""
    degree_12 = 3
    n12 = 6
    n3 = 4
    quad = SO3([
        ('s2-gauss-lebedevlaikov', {'degree': degree_12}),
        ('trapezoid', {'size': n3}),
    ])
    # test submethod properties:
    assert quad._methods[0].degree == degree_12
    assert quad._methods[1].degree == 1
    assert quad._methods[0].size == n12
    assert quad._methods[1].size == n3
    assert len(quad._method_weights) == 2
    assert len(quad._method_points) == 2
    assert quad._method_points[0].shape == (2, n12)
    assert quad._method_points[1].shape == (1, n3)

    # test private properties
    assert quad._ndims == [2, 1]
    assert quad._dims == [(0, 1), (2,)] # do I really need this?
    assert quad._sizes == [n12, n3]
    assert quad._degrees == [degree_12, 1]
    assert quad._points.size == (3*n12*n3)
    assert quad._points.shape == (3, n12*n3)

    # test public properties:
    assert quad.dim == 3
    assert quad.size == n12*n3
    assert quad.shape == (n12, n3)
    assert quad.weights.shape == (n12*n3,)
    assert quad.angles.shape == (3, n12*n3)
    assert quad.quaternions.shape == (n12*n3, 4)
    assert np.all(quad.angles == quad._points)
    #assert quad.grid.angles.shape == (3, n12*n3)
    #assert quad.grid.quaternions.shape == (n12*n3, 4) ???
    assert quad.weights.reshape(quad.shape).shape == (n12, n3)
    assert quad.angles.reshape(3, *quad.shape).shape == (3, n12, n3)
    assert quad.quaternions.reshape(*quad.shape, 4).shape == (n12, n3, 4)


def test_integration():
    """Test integration for different methods.
    Euler angles are used in this test.
    """
    quad = SO3([
        ('trapezoid', dict(size=7)),
        ('gausslegendre', dict(size=9)),
        ('trapezoid', dict(size=5)),
    ], polar_sampling='cos')
    assert quad.integrate(f_abg) == pytest.approx(1)
    f_samples = f_abg(*quad.angles)
    assert quad.integrate(f_samples) == pytest.approx(1)

    quad = SO3([
        ('trapezoid', dict(size=7)),
        ('gausslegendre', dict(size=9)),
        ('trapezoid', dict(size=5)),
    ], polar_sampling='angle')
    assert quad.integrate(f_abg) == pytest.approx(1)
    f_samples = f_abg(*quad.angles)
    assert quad.integrate(f_samples) == pytest.approx(1)

    quad = SO3([
        ('trapezoid', {'size':11}),
        ('s2-gauss-lebedevlaikov', {'degree':5}),
    ])
    assert quad.integrate(f_abg) == pytest.approx(1)
    f_samples = f_abg(*quad.angles)
    assert quad.integrate(f_samples) == pytest.approx(1)

    quad = SO3([
        ('s2-gauss-lebedevlaikov', {'degree':5}),
        ('trapezoid', {'size':11}),
    ])
    assert quad.integrate(f_abg) == pytest.approx(1)
    f_samples = f_abg(*quad.angles)
    assert quad.integrate(f_samples) == pytest.approx(1)

    quad = SO3([
        ('so3-montecarlo', {'size':500, 'seed':0}),
    ])
    assert quad.integrate(f_abg) == pytest.approx(1)
    f_samples = f_abg(*quad.angles)
    assert quad.integrate(f_samples) == pytest.approx(1)


def test_integration_quaternions():
    """Test integration with quaternions."""
    quad = SO3([
        ('s2-gauss-lebedevlaikov', {'degree':5}),
        ('trapezoid', {'size':11}),
    ])
    f_samples = f_q(quad.quaternions)
    assert quad.integrate(f_samples) == pytest.approx(1)


def test_integration_axis():
    """Test axis keyword in integration method."""
    quad = SO3([
        ('s2-gauss-lebedevlaikov', {'degree':7}),
        ('trapezoid', {'size':11}),
    ])
    f_vector = f_vector_q(quad.quaternions)
    result = quad.integrate(f_vector, axis=0)
    assert result.shape == (3, 3)
    assert result == pytest.approx(1)
