.. _install:

Installation instructions
-------------------------


Prerequisites and compatibility
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This package uses `NumPy`_ and `SciPy`_ for data handling and numerical operations.

.. todo: do I use scipy?

Some modules make use of Mike Boyle's `quaternionic`_ package.
This allows easy interoperability with `spherical`_.

.. todo: mention spherical only in further resources

.. _NumPy: https://numpy.org/
.. _SciPy: https://scipy.org/
.. _quaternionic: https://github.com/moble/quaternionic
.. _spherical: https://github.com/moble/spherical


Installation
^^^^^^^^^^^^

This package is available on `PyPi`_. You can install it with

.. code-block:: bash

    python -m pip install openquad

This will install the most recent stable release into your active
`environment`_ along with any dependencies.

.. _PyPi: https://pypi.org/...
.. _environment: https://docs.python.org/3/glossary.html#term-virtual-environment

If you intend to contribute to this package, follow the instructions for a
:ref:`development installation <dev_install>`.
