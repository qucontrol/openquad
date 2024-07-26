"""Set up the environment for doctests.

This file is automatically evaluated by pytest. It ensures that we can write
doctests without distracting import statements in the doctest.
"""
import pytest
import numpy as np


@pytest.ficture(autouse=True)
def set_doctest_env(doctest_namespace):
    doctest_namespace['np'] = np
