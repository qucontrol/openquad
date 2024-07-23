import numpy as np
from scipy import special, integrate

import spherical
import quaternionic

import pytest
from contextlib import nullcontext as does_not_raise

from openquad.monte_carlo import MonteCarloSO3
from openquad.graef import GraefSO3EqualWeight, GraefSO3Gauss
from openquad.womersley import WomersleySO3EqualWeight
from openquad.fibonacci import ZCW3
from openquad.karney import KarneySO3


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


def f_abg(alpha, beta, gamma):
    """Test integrand for Euler angles.
    Integral should yield 1.
    """
    assert np.shape(alpha) == np.shape(beta) == np.shape(gamma)
    quaternions = quaternionic.array.from_euler_angles(alpha, beta, gamma)
    return f_q(quaternions)


def test_test_function():
    """Make sure that the integral over the test function yields 1."""
    assert integrate.tplquad(lambda a, b, g: np.sin(b)*f_abg(a, b, g).real, 0, 2*np.pi, 0, np.pi, 0, 2*np.pi)[0] == pytest.approx(1)
    #assert integrate.tplquad(lambda a, b, g: np.sin(b)*f_abg(a, b, g).imag, 0, 2*np.pi, 0, np.pi, 0, 2*np.pi, epsabs=1e-1)[0] == pytest.approx(0)


@pytest.mark.parametrize(
    "method, size, p_acc",
    [
        (MonteCarloSO3, 10, None),
        (GraefSO3Gauss, 4, 1),
        (GraefSO3EqualWeight, 4, 1),
        #(WomersleySO3EqualWeight, 4, 1),
        #(ZCW3, 4, None),
        (KarneySO3, 24, None),
    ]
)
def test_initialization(method, size, p_acc):
    # missing required arguments
    with pytest.raises(TypeError):
        quad = method()

    # duplicate arguments
    if p_acc is not None:
        with pytest.raises(TypeError, match="not both"):
            quad = method(p_acc=3, size=6)

    # unexpected arguments
    with pytest.raises(TypeError, match="unexpected keyword"):
        quad = method(size=size, unexpected=True)

    # invalid number of points or accuracy
    with pytest.raises(ValueError, match="larger than zero"):
        quad = method(size=0)
    if p_acc is not None:
        with pytest.raises(ValueError, match="not available"):
            quad = method(size=1)
        with pytest.raises(ValueError, match="not available"):
            quad = method(p_acc=0)

    # valid initializations
    with does_not_raise():
        quad = method(size=size)
    if p_acc is not None:
        with does_not_raise():
            quad = method(p_acc=p_acc)


@pytest.mark.parametrize(
    "method, size, p_acc, tol, kwargs",
    [
        (MonteCarloSO3, 1000, None, 1e-2, {'seed':0, 'method':'quaternions'}),
        (MonteCarloSO3, 5000, None, 1e-2, {'seed':0, 'method':'angles'}),
        (KarneySO3, 24, None, 1e-8, {}),
        (GraefSO3Gauss, 23, 3, 1e-8, {}),
        (GraefSO3EqualWeight, 24, 3, 1e-8, {}),
        #(WomersleySO3EqualWeight, 10, 10, 1e-8, {}),
    ]
)
def test_fields_and_integration(method, size, p_acc, tol, kwargs):
    quad = method(size=size, **kwargs)
    
    # test fields
    assert quad.dim == 3
    assert quad.p_acc == p_acc
    assert quad.size == size

    # test weights
    assert quad.weights.shape == (size,)
    assert np.sum(quad.weights) == pytest.approx(8*np.pi**2)

    # test grid dimension
    assert quad.points.shape == (3, size)
    assert np.all(quad.angles == quad.points)

    # test for callable
    assert quad.integrate(f_abg) == pytest.approx(1, rel=tol)
    # test for samples
    alpha, beta, gamma = quad.points
    f_samples = f_abg(alpha, beta, gamma)
    assert quad.integrate(f_samples) == pytest.approx(1, rel=tol)


@pytest.mark.parametrize(
    "method, tol, kwargs",
    [
        #TODO: set accuracy
        (GraefSO3Gauss, 1e-9, {}),
        (GraefSO3EqualWeight, 1e-6, {}),
        #(WomersleySO3EqualWeight, 1e-12, {}),
    ]
)
def test_poly_acc(method, tol, kwargs):
    max_p_acc = 30
    wigner = spherical.Wigner(max_p_acc)
    available_p_accs = method._available_p_accs
    tested_p_accs = available_p_accs[available_p_accs <= max_p_acc]
    for p_acc in tested_p_accs:
        quad = method(p_acc=p_acc, **kwargs)
        # First, test that the weights are normalized:
        assert np.sum(quad.weights) == pytest.approx(8*np.pi**2)
        # Then, test the integration of all polynomials:
        q = quaternionic.array.from_euler_angles(*quad.angles)
        Dlmm = wigner.D(q)[:, :wigner.Dindex(p_acc, p_acc, p_acc)+1]
        results = quad.integrate(Dlmm, axis=0)
        truevals = np.concatenate(
            ([8*np.pi**2], np.zeros(Dlmm.shape[-1]-1))
        )
        assert np.abs(results - truevals).max() == pytest.approx(0, abs=tol)
