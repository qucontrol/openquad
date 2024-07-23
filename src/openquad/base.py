from pathlib import Path
from abc import ABC, abstractmethod
import numpy as np


QUAD_DATA_DIR = (Path.cwd() / __file__).parent / 'data'
"""Directory where quadrature grids are stored."""


class AtomicQuadrature(ABC):

    @property
    @abstractmethod
    def dim(self):
        """Dimension of the quadrature method.
        This should be an integer.
        This should be a class attribute.
        """
        raise NotImplementedError("Missing class attribute `dim`")

    @property
    def size(self):
        """Total number of sample points.
        This must be an integer.
        """
        try:
            return self._size
        except AttributeError:
            msg = f"{type(self)} missing attribute `size`"
            raise NotImplementedError(msg)

    @size.setter
    def size(self, value):
        self._check_size_allowed(value)
        self._size = value

    def _check_size_allowed(self, value):
        if not np.issubdtype(type(value), np.integer):
            raise TypeError(f"`size` must be an integer, not {type(value)}")
        if value < 1:
            raise ValueError("`size` must be larger than zero")
        if hasattr(self, "_available_sizes"):
            if value not in self._available_sizes:
                msg = (f"`size` {value} not available. Possible values are:"
                       f"{self._available_sizes}")
                raise ValueError(msg)

    @property
    def p_acc(self):
        """Polynomial accuracy of the quadrature method.
        If the method does not have a polynomial accuracy, this should be
        `None`. Otherwise this should be an integer.
        """
        try:
            return self._p_acc
        except AttributeError:
            msg = f"{type(self)} missing attribute `p_acc`"
            raise NotImplementedError(msg)

    @p_acc.setter
    def p_acc(self, value):
        self._check_poly_acc_allowed(value)
        self._p_acc = value

    def _check_poly_acc_allowed(self, value):
        if value is not None and value < 0:
            raise ValueError("`p_acc` must not be negative")
        if hasattr(self, "_available_p_accs"):
            if value not in self._available_p_accs:
                msg = (f"`p_acc` {value} not available. Possible values are:"
                       f"{self._available_p_accs}")
                raise ValueError(msg)

    def __init__(self, *pw_args, jacobian=None, **pw_kwargs):
        """Base initialization of the quadrature method.

        Here, the quadrature weights and sample points are initialized.
        """
        x, w = self._points_weights(*pw_args, **pw_kwargs)
        if jacobian is not None:
            w = self._apply_jacobian(x, w, jacobian)
        self.points = x
        self.weights = w
        #TODO: ensure that points are contiguous?
        # Maybe necessary for reliable reshaping?

    @abstractmethod
    def _points_weights(self, **kwargs):
        """TODO!

        Returns
        -------
        points : ndarray
            Sample points.
        weights : ndarray
            Quadrature weights.
        """
        msg = f"{self.__class__} missing method _points_weights"
        raise NotImplementedError(msg)

    def _load_points_weights(self, folder, type, suffix=''):
        if self.p_acc is not None:
            filename = f"{type}_size{self.size}_acc{self.p_acc}{suffix}.npy"
        else:
            filename = f"{type}_size{self.size}{suffix}.npy"
        data = np.load(QUAD_DATA_DIR / folder / filename,
                       allow_pickle=False, fix_imports=False)
        if data.shape[0] not in [self.dim, self.dim+1]:
            raise TypeError(
                f"Read data has wrong shape {data.shape} for dim {self.dim}"
            )
        points = data[0:self.dim]
        if data.shape[0] == self.dim + 1:
            weights = data[-1]
        else:
            weights = np.ones(self.size) / self.size
        return points, weights

    def _apply_jacobian(self, x, w, jacobian):
        if self.dim == 1: x = np.expand_dims(x, 0)
        try:
            w = w * jacobian(*x)
        except:
            raise TypeError("Cannot apply jacobian")
        return w

    def integrate(self, f, f_kwargs={}, axis=-1):
        """Integrate `f`.

        `f` can be a callable, in which case `f(*points, **f_kwargs)` is integrated.
        """
        if callable(f):
            if self.dim == 1:
                f_samples = f(self.points, **f_kwargs)
            else:
                f_samples = f(*self.points, **f_kwargs)
        else:
            f_samples = f
        return self._integrate_sample(f_samples, axis=axis)

    def _integrate_sample(self, f_samples, axis=-1):
        """Integrate an array of sample points.

        This is the default implementation as a sum of sample points with their
        corresponding weights.
        Subclasses may overwrite this method, if necessary.
        """
        f_samples = np.moveaxis(f_samples, axis, -1)
        return np.sum(f_samples * self.weights, axis=-1)


class QuadratureWithoutPolyAcc(AtomicQuadrature):

    def __init__(self, size, *args, **kwargs):
        self.size = size
        self.p_acc = None
        super().__init__(*args, **kwargs)


class QuadratureWithFixedPolyAcc(AtomicQuadrature):

    def __init__(self, size, *args, **kwargs):
        self.size = size
        self.p_acc = self._get_p_acc()
        super().__init__(*args, **kwargs)

    @abstractmethod
    def _get_p_acc(self):
        """Return the polynomial accuracy of this quadrature method.

        For quadrature methods where the polynomial accuracy does not depend on
        the size of the quadrature grid, this method needs to be defined and
        return a single integer.

        """
        msg = f"{self.__class__} missing method _get_p_acc"
        raise NotImplementedError(msg)


class QuadratureWithPolyAcc(AtomicQuadrature):

    def __init__(self, *args, size=None, p_acc=None, **kwargs):
        if size is not None and p_acc is None:
            self.size = size
            self.p_acc = self._get_p_acc(size)
        elif p_acc is not None and size is None:
            self.p_acc = p_acc
            self.size = self._get_size(p_acc)
        else:
            raise TypeError("Either `size` or `p_acc` must be given, not both")
        super().__init__(*args, **kwargs)

    def _get_size(self, p_acc):
        """Return the size corresponding to the given polynomial accuracy.

        If the quadrature class has the attribute `_available_sizes`, search
        the size in this list. If a functional relation between the size and
        the polynomial accuracy of the method exits, this method needs to be
        overwritten in the corresponding child class.

        """
        if hasattr(self, '_available_sizes'):
            selection = np.nonzero(self._available_p_accs == p_acc)
            return self._available_sizes[selection].min()
        else:
            msg = "f{self.__class__} missing method _get_size"
            raise NotImplementedError(msg)

    def _get_p_acc(self, size):
        """Return the polynomial accuracy corresponding to the given size.

        If the quadrature class has the attribute `_available_p_accs`, search
        the accuracy in this list. If a functional relation between the size
        and the polynomial accuracy of the method exits, this method needs to
        be overwritten in the corresponding child class.

        """
        if hasattr(self, '_available_p_accs'):
            selection = np.nonzero(self._available_sizes == size)
            return self._available_p_accs[selection].max()
        else:
            msg = f"{self.__class__} missing method _get_p_acc"
            raise NotImplementedError(msg)


class AtomicR1Quadrature(AtomicQuadrature):

    dim = 1

    @property
    @abstractmethod
    def _has_endpoints(self):
        """Flag indicating, if the quadrature grid contains the endpoints of
        the integration interval.
        """
        pass

    def __init__(self, a, b, *args, periodic=False, **kwargs):
        if a == b:
            raise ValueError("`a` must not be equal to `b`")
        self._periodic = periodic
        super().__init__(*args, a=a, b=b, **kwargs)
        if self._periodic and self._has_endpoints: # drop endpoint
            self.weights[0] += self.weights[-1]
            self.weights = self.weights[:-1]
            self.points = self.points[:-1]


class AtomicS2Quadrature(AtomicQuadrature):

    dim = 2

    def __init__(self, *args, swap_angles=False, **kwargs):
        super().__init__(*args, **kwargs)
        if swap_angles:
            self.points = np.stack((self.points[1], self.points[0]))
            self._swapped_angles = True
        else:
            self._swapped_angles = False
        self.angles = self.points


class AtomicSO3Quadrature(AtomicQuadrature):

    dim = 3

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.angles = self.points


class NoQuadrature(ABC):

    dim = 1

    def __init__(self, points, *args, **kwargs):
        """Dummy method that does not perform quadrature.

        Can be used for not integrating certain dimensions and instead only
        evaluate the intgerand at certain values.

        Parameters
        ----------
        values : array-like
            Grid points on which to evaluate the integrand.

        """
        points = np.asarray(points)
        self.size = points.size
        self.p_acc = None
        self.points = points
        self.weights = np.ones(points.size)

    def integrate(self, f, f_kwargs={}, axis=-1):
        """Dummy integration routine that does not perform integration.

        `f` can be a callable, in which case `f(points, **f_kwargs)` is evaluated.
        """
        if callable(f):
            f_samples = f(self.points, **f_kwargs)
        else:
            f_samples = f
        return f_samples
