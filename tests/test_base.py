import numpy as np

import pytest
from contextlib import nullcontext as does_not_raise

from openquad.quad import TopologicalQuadrature


def test_all_method_aliases_are_valid():
    """Test that all method aliases point to valid method names."""
    method_aliases = TopologicalQuadrature._method_aliases
    method_names = TopologicalQuadrature._available_methods
    for method in method_aliases.values():
        assert method in method_names.keys()
