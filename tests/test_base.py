import pytest
from contextlib import nullcontext as does_not_raise

from openquad.geometries import GeometryQuadrature


def test_all_method_aliases_are_valid():
    """Test that all method aliases point to valid method names."""
    method_aliases = GeometryQuadrature._method_aliases
    method_names = GeometryQuadrature._available_methods
    for method in method_aliases.values():
        assert method in method_names.keys()
