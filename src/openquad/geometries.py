"""Module for integrating angular functions.

This module aims for usability and flexibility...
#TODO

"""
from abc import ABC, abstractmethod
from copy import deepcopy
import numpy as np

import quaternionic

from .base import NoQuadrature
from .grid import xyz_from_angles, S2Grid, SO3Grid
from .newton_cotes import CompositeTrapezoid, CompositeSimpson, Romberg
from .gauss import GaussLegendre, GaussLobattoLegendre
from .lebedev import LebedevLaikov
from .graef import GraefS2Gauss, GraefS2Design, GraefSO3Gauss, GraefSO3EqualWeight
from .karney import KarneySO3
from .womersley import WomersleyS2Design, WomersleySO3EqualWeight
from .fibonacci import FibonacciSphere, ZCW2, ZCW3
from .monte_carlo import MonteCarloR1, MonteCarloS2, MonteCarloSO3

π = np.pi


# TODO: move to submodule?
class GeometryQuadrature(ABC):

    _available_methods = {
        'none': NoQuadrature,
        # 1d methods:
        'composite trapezoid': CompositeTrapezoid,
        'composite simpson': CompositeSimpson,
        'romberg': Romberg,
        'gauss-legendre': GaussLegendre,
        'gauss-lobatto-legendre': GaussLobattoLegendre,
        'monte-carlo-1d': MonteCarloR1,
        # two-angle methods:
        'lebedev-laikov': LebedevLaikov,
        's2-gauss-graef': GraefS2Gauss,
        's2-design-graef': GraefS2Design,
        's2-design-womersley': WomersleyS2Design,
        'fibonacci-sphere': FibonacciSphere,
        'zcw2': ZCW2,
        'monte-carlo-s2': MonteCarloS2,
        # three-angle methods:
        'so3-covering-karney': KarneySO3,
        'so3-gauss-graef': GraefSO3Gauss,
        'so3-equalweight-graef': GraefSO3EqualWeight,
        'so3-equalweight-womersley': WomersleySO3EqualWeight,
        'zcw3': ZCW3,
        'monte-carlo-so3': MonteCarloSO3,
    }
    _method_aliases = {
        # dummy quadrature
        'none': 'none',
        # composite trapezoid rule
        'composite trapezoid': 'composite trapezoid',
        'trapezoid': 'composite trapezoid',
        'trapz': 'composite trapezoid',
        # composite simpson rule
        'composite simpson': 'composite simpson',
        'simpson': 'composite simpson',
        'simps': 'composite simpson',
        # Romberg integration
        'romberg': 'romberg',
        'romb': 'romberg',
        # Gauss-Legendre quadrature
        'gauss-legendre': 'gauss-legendre',
        'gl': 'gauss-legendre',
        # Gauss-Lobatto-Legendre
        'gauss-lobatto-legendre': 'gauss-lobatto-legendre',
        'gll': 'gauss-lobatto-legendre',
        # Lebedev Laikov
        'lebedev-laikov': 'lebedev-laikov',
        'lebedev': 'lebedev-laikov',
        # 1d Monte-Carlo
        'monte-carlo-1d': 'monte-carlo-1d',
        'mc1': 'monte-carlo-1d',
        # two-angle Monte-Carlo
        'monte-carlo-s2': 'monte-carlo-s2',
        'mcs2': 'monte-carlo-s2',
        # three-angle Monte-Carlo
        'monte-carlo-so3': 'monte-carlo-so3',
        'mcso3': 'monte-carlo-so3',
        # other methods without aliases
        's2-gauss-graef': 's2-gauss-graef',
        's2-design-graef': 's2-design-graef',
        's2-design-womersley': 's2-design-womersley',
        'fibonacci-sphere': 'fibonacci-sphere',
        'zcw2': 'zcw2',
        'so3-covering-karney': 'so3-covering-karney',
        'so3-gauss-graef': 'so3-gauss-graef',
        'so3-equalweight-graef': 'so3-equalweight-graef',
        'so3-equalweight-womersley': 'so3-equalweight-womersley',
        'zcw3': 'zcw3',
    }

    @property
    @abstractmethod
    def dim(self):
        """Dimension of the space that this quadrature method is working on."""
        pass

    def __init__(self, method_specs):
        self._method_points = []
        self._method_weights = []
        self._methods = []
        self._ndims = []
        self._dims = []
        self._sizes = []
        self._degrees = []

        # Resolve and check selected quadrature methods:
        method_objects = []
        method_options = []
        if len(method_specs) < 1:
            raise ValueError("Method speficication must not be empty")
        for i, method_spec in enumerate(method_specs):
            method_str, options = deepcopy(method_spec)
            method_objects.append(self._select_method_object(method_str))
            method_options.append(options)
        self._check_method_objects(method_objects)

        # Finalize data structures for weights and quadrature grids and define
        # remaining attributes:
        self._initialize_quadrature_methods(method_objects, method_options)
        self._construct_weights_and_meshgrids()
        self.shape = tuple(self._sizes)
        self.size = self.weights.size

    def _resolve_method_name(self, method_str):
        """Map the given `method_str` to the full method name."""
        if method_str.lower() in self._method_aliases:
            method_name = self._method_aliases[method_str.lower()]
        else:
            raise ValueError(f"Invalid method name: {method_str}")
        return method_name

    def _select_method_object(self, method_str):
        method_name = self._resolve_method_name(method_str)
        if method_name in self._available_methods:
            method_object = self._available_methods[method_name]
        else:
            raise ValueError(f"Method not available: {method_name}")
        return method_object

    def _check_method_objects(self, method_objects):
        total_dim = np.sum([method.dim for method in method_objects])
        if total_dim != self.dim:
            raise ValueError("Invalid combination of methods")

    def _initialize_quadrature_methods(self, method_objects, method_options):
        """Fill options and initialize quadrature method objects."""
        collected_dims = 0
        for method_object, options in zip(method_objects, method_options):

            # figure out, on which dimensions the current method is acting:
            self._dims.append(tuple(
                np.arange(method_object.dim) + collected_dims
            ))
            collected_dims += method_object.dim

            # fill remaining options
            options = self._check_and_fill_options(options, self._dims[-1])
            method = method_object(**options)
            if method.points.ndim == 1:
                self._method_points.append(np.expand_dims(method.points, 0))
            else:
                self._method_points.append(method.points)
            self._method_weights.append(method.weights)
            self._methods.append(method)
            self._ndims.append(method.dim)
            self._sizes.append(method.size)
            self._degrees.append(method.degree)

    def _check_and_fill_options(self, options, dims):
        """Check the integrations options and fill remaining fields."""
        a, b, jacobian = self._get_interval(dims)
        for key, value in zip(
                ['a', 'b', 'jacobian'],
                [a, b, jacobian]
                ):
            if key not in options and value == 'user':
                raise ValueError(f"Option `{key}` needs to be specified by the user")
            elif key in options and value == 'auto':
                raise ValueError(f"Option `{key}` is set automatically and must not be given")
            elif key not in options and value != 'auto':
                options[key] = value
        # use periodic boundaries if possible
        if (
            (isinstance(self, S2) and dims == (1,)) or
            (isinstance(self, SO3) and dims in [(0,), (2,)])
        ):
            options['periodic'] = True
        # determine ordering of angles for spherical quadrature schemes:
        if self.dim == 3 and len(dims) == 2 and dims == (0, 1):
            options['swap_angles'] = True
            # assuming swap_angles=False by default in 2d quadrature classes
        return options

    def _construct_weights_and_meshgrids(self):
        """Construct combined quadrature weights and grid."""
        self.weights = np.prod(
            np.meshgrid(*self._method_weights, indexing='ij', copy=False),
            axis=0,
        ).flatten()
        # For creating the quadrature meshgrid, I need to handle the cases
        # separately depending on the dimension of the subgrids:
        if self._ndims[0] == self.dim:
            # the points form already the meshgrid
            meshgrid = self._method_points[0]
            # TODO: make sure it is contigous?
            # TODO: could be merged with the case below?
        elif np.all(np.asarray(self._ndims) == 1):
            # If all subgrids have dimension one, it is easy:
            meshgrid = np.meshgrid(*self._method_points, indexing='ij', copy=False),
        elif self._ndims[0] == 1 and self._ndims[1] == 2:
            # If some subgrids have more than one dimension, I need to build
            # the meshgrid gradually:
            p0, p1 = np.meshgrid(self._method_points[0][0], self._method_points[1][0],
                                 indexing='ij', copy=False)
            p0, p2 = np.meshgrid(self._method_points[0][0], self._method_points[1][1],
                                 indexing='ij', copy=False)
            meshgrid = np.stack((p0, p1, p2), axis=0)
        elif self._ndims[0] == 2 and self._ndims[1] == 1:
            p0, p2 = np.meshgrid(self._method_points[0][0], self._method_points[1][0],
                                 indexing='ij', copy=False)
            p1, p2 = np.meshgrid(self._method_points[0][1], self._method_points[1][0],
                                 indexing='ij', copy=False)
            meshgrid = np.stack((p0, p1, p2), axis=0)
        else:
            raise ValueError("Error in constructing quadrature grid")
        #TODO: is there a more elegant solution than the one above?
        self._points = np.reshape(meshgrid, (self.dim, -1))
        #TODO: do I need to set copy=True?
        # This could be necessary to make sure that reshaping is reliable,
        # independent from the data's memory layout

    @abstractmethod
    def _get_interval(self, dims):
        """Determine the integration boundaries and the Jacobian.

        This function needs to be defined for every quadrature subclass.
        For all parameters that are to be set by the user in the `options`
        dictionary, this function needs to return ``'user'``.
        For all parameters that are determined by quadpy automatically,
        this function needs to return ``'auto'``.
        For all optional parameters, this function needs to return `None`.

        Parameters
        ----------
        dims : tuple
            Tuple containing the dimensions that a quadrature scheme is working
            on.

        Returns
        -------
        a : float
            Lower bound of the integration interval.
        b : float
            Upper bound of the integration interval.
        jacobian : callable
            Jacobian factor for the integration interval.
            This needs to be the sample points as single argument.

        """
        #TODO: raise NotImplementedError ?
        # see https://docs.python.org/3/library/exceptions.html#NotImplementedError
        a = None
        b = None
        jacobian = None
        return a, b, jacobian

    def integrate(self, f, f_kwargs={}, axis=-1, dim=None):#, squeeze=True):#, order=None, only=None):
        """Integrate using the given samples of a function.

        This method integrates the data with the integration methods provided
        during initialization of the quadrature object.

        For each integration method, `y` is integrated along each 1d slice on
        the given axis.

        Parameters
        ----------
        y : callable or array_like
            Representation of the integrand. Either an array-like or an callable.
        axis : sequence of ints, optional
            Axes along which to perform the integration. The number of axes provided must match the number
            of integration methods hold by the quadrature object.
            By default, start with the last axis of `y`.
        squeeze : bool, optional
            Some SciPy integration method return zero for integration of an
            axis with length one. If True, axes of length one are squeezed,
            i.e. the result of the integration is the single value of the
            integrand along this axes. Axes that are not integrated are not
            affected by the squeeze.
        dim : int or sequence of ints, optional
            Specifies which dimensions are to be integrated.
            If None, integrate all dimensions.

        Returns
        -------
        result : float or ndarray
            The result of the integration. If the dimension of `y` equals the
            number of integration methods, the result is a float.

        """
        if callable(f):
            if self.dim == 1:
                f_array = f(self._points, **f_kwargs)
            else:
                f_array = f(*self._points, **f_kwargs)
        else:
            f_array = f

        # options for dim:
        # 0
        # (0, 1)
        # ((0, 1), 2) #no!
        if dim is not None: # integrate speficic dimensions
            # This blew up in complexity; I won't allow it for now.
            raise NotImplementedError
            if np.ndim(dim) == 0: dim = (dim,)
            assert len(dim) > 0 and len(dim) <= self.dim
            dim = sorted(dim) # the implementation below depends on it!
            method_index = []
            for d in dim:
                # figure out, which method to use for dimension d; after the
                # loop, it will be the method at index m
                for m in range(len(self._ndims)):
                    if d in self._dims[m]:
                        break
                method_index.append(m)
            method_index, ndims = np.unique(method_index, return_counts=True)
            # check that the user wants to integrate over all necessary
            # dimensions
            for i in range(len(method_index)):
                if ndims[i] != self._ndims[i]:
                    raise ValueError(f"Invalid selection of dimensions")
            # Now determine the axis for each dimension to integrate
            axes = []
            reduced_axes = 0
            for m in np.flip(method_index):
                axes.append(-len(self._methods)+m+reduced_axes)
                reduced_axes += 1

        # the following was an older try:
        #if dim is not None:
        #    if np.ndim(dim) == 0: dim = (dim,)
        #    assert len(dim) > 0 and len(dim) <= self.dim
        #    # reshape the sample points
        #    y = np.reshape(y, y.shape[:axis]+self.shape+y.shape[axis:])
        #    result = y
        #    for d in dim:
        #        result = self._methods[d].integrate(y, axis=axis+d)

        else: # integrate all dimensions of a sample
            method_index = list(range(len(self._methods)))
            axes = np.repeat(-1, len(self._methods))

        f_array = np.moveaxis(f_array, axis, -1)
        shape = f_array.shape[:-1]+self.shape
        result = f_array.reshape(shape)
        # iterate over angles/methods, starting with the last:
        for i_m, m in enumerate(np.flip(method_index)):
            ndim_before = result.ndim
            result = self._methods[m].integrate(result, axis=axes[i_m])
            if result.ndim != ndim_before - 1:
                # No reduction has been performed. Hence, we need to skip this
                # axis for the remaining iterations
                axes = axes - 1
        return result

    @abstractmethod
    def write_grid(self, filename, weights=False):
        pass


class SO3(GeometryQuadrature):

    dim = 3

    def __init__(self, method_specs, polar_sampling='cos'):
        allowed_samplings = ['angle', 'cos']
        if polar_sampling not in allowed_samplings:
            raise ValueError(
                f"`polar_sampling needs to be one of {allowed_samplings}"
            )
        self._polar_sampling = polar_sampling
        super().__init__(method_specs)
        # in case we sample the polar integral in terms of cos(angle),
        # transform the sample points back to the angle:
        if self._polar_sampling == 'cos' and len(self._methods) > 2:
            self._points[1] = (
                np.flip(np.arccos(self._points[1]), axis=0)
            )
        self.angles = self._points
        self.quaternions = quaternionic.array.from_euler_angles(*self.angles)

    def _get_interval(self, dims):
        """Determine the integration boundaries and the Jacobian."""
        if len(dims) == 3: # SO3 quadrature scheme
            return 'auto', 'auto', 'auto'
        elif len(dims) == 2: # spherical quadrature scheme
            return 'auto', 'auto', 'auto'
        elif len(dims) == 1 and dims[0] == 1: # 1d method for beta
            if self._polar_sampling == 'angle':
                return 0, π, np.sin
            elif self._polar_sampling == 'cos':
                return -1, 1, None
        elif len(dims) == 1 and (dims[0] == 0 or dims[0] == 2) : # 1d method for alpha or gamma
            return 0, 2*π, None
        else:
            raise ValueError("Cannot determine integration interval")

    def write_grid(self, filename, weights=False):
        grid = SO3Grid('Euler angles', self.angles)
        if weights:
            grid.save(filename, 'Euler angles', weights=self.weights)
        else:
            grid.save(filename, 'Euler angles')


class S2(GeometryQuadrature):

    dim = 2

    def __init__(self, method_specs, polar_sampling='cos'):
        allowed_samplings = ['angle', 'cos']
        if polar_sampling not in allowed_samplings:
            raise ValueError(
                f"`polar_sampling needs to be one of {allowed_samplings}"
            )
        self._polar_sampling = polar_sampling
        super().__init__(method_specs)
        # in case we sample the polar integral in terms of cos(angle),
        # transform the sample points back to the angle:
        if self._polar_sampling == 'cos' and len(self._methods) > 1:
            self._points[0] = (
                np.flip(np.arccos(self._points[0]), axis=0)
            )
        self.angles = self._points
        self.xyz = xyz_from_angles(self.angles, axis=0)

    def _get_interval(self, dims):
        """Determine the integration boundaries and the Jacobian."""
        if len(dims) == 2: # Spherical quadrature scheme
            return 'auto', 'auto', 'auto'
        elif len(dims) == 1 and dims[0] == 0: # Polar angle
            if self._polar_sampling == 'angle':
                return 0, π, np.sin
            elif self._polar_sampling == 'cos':
                return -1, 1, None
        elif len(dims) == 1 and dims[0] == 1: # Azimuthal angle
            return 0, 2*π, None
        else:
            raise ValueError("Cannot determine integration interval")

    def write_grid(self, filename, weights=False):
        grid = S2Grid('angles', self.angles)
        if weights:
            grid.save(filename, 'angles', weights=self.weights)
        else:
            grid.save(filename, 'angles')


class Rn(GeometryQuadrature):

    dim = None #TODO: can this be more elegant?
    # Maybe, make dim an instance attribute for consistency

    def __init__(self, method_specs):
        self.dim = len(method_specs)
        # Note: determining the dimension this way, effectively
        # allows to use only 1d methods for this class.
        super().__init__(method_specs)
        self.points = self._points
        self.r = self._points #TODO: deprecate self.r?

    def _get_interval(self, dims):
        """Determine the integration boundaries and the Jacobian."""
        # All integration bounds need to be set by the user in this quadrature class.
        # Specification of Jacobians is optional.
        return 'user', 'user', None

    def write_grid(self, filename, weights=False):
        msg = "Writing grids not yet implemented for this class"
        raise NotImplementedError(msg)
