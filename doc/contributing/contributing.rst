Ready to contribute?
--------------------


Development workflow
^^^^^^^^^^^^^^^^^^^^

Follow `Aaron Meurer's description of the git workflow`_.

.. _Aaron Meurer's description of the git workflow: https://www.asmeurer.com/git-workflow/

Briefly:

1. Create a personal fork of this repository.
2. Clone this repo from your fork.
3. Create a branch for your contribution.
4. Make your changes and commit them (testing locally).
5. Push changes to the branch on your remote.
6. Open a pull request via the GitHub website of your fork.

If you are a member of the `qucontrol organization`_, you can omit the first
step and work in the parent GitHub repository directly.

.. _qucontrol organization: https://github.com/qucontrol


Commit messages should be clear and concise. Use the following template:

.. code-block:: none

   Short (72 chars or less) summary

   More detailed explanatory text. Wrap it to 72 characters. The blank
   line separating the summary from the body is critical (unless you omit
   the body entirely).

   Write your commit message in the imperative: "Fix bug" and not "Fixed
   bug" or "Fixes bug." The git commit subject line should always be able
   to complete the sentence
   "If applied, this commit will <your subject line here>".

   Reference any issue that is being addressed in the commit, as e.g. "#1"
   for issue #1. If the commit closes an issue, state this on the last line
   of the message (see below). This will automatically close the issue on
   Github as soon as the commit is pushed there.

   Closes #1

See `Linking a pull request to an issue`_ for details on how to reference
issues in commits.

.. _Linking a pull request to an issue: https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue


Relases follow `semantic versioning`_ and need to be compatible with `PyPA
version specifiers`_ (former `PEP440`_).

.. _semantic versioning: https://semver.org/
.. _PyPA version specifiers: https://packaging.python.org/en/latest/specifications/version-specifiers/#version-specifiers
.. _PEP440: https://peps.python.org/pep-0440/


Styles and conventions
^^^^^^^^^^^^^^^^^^^^^^

Code needs to be compatible with `PEP 8`_. It is formatted with `Ruff`_, which
follows the Black code style, with Ruff's default line length of 88
characters.

Docstrings are written according to the `NumPy Style Guide`_, e.g. see this
`comprehensive example`_.

.. _PEP 8: https://peps.python.org/pep-0008/
.. _Ruff: https://docs.astral.sh/ruff/formatter/
.. _NumPy Style Guide: https://numpydoc.readthedocs.io/en/latest/format.html
.. _comprehensive example: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html


.. _dev_install:

Development installation
^^^^^^^^^^^^^^^^^^^^^^^^

We recommend to use a `virtual environment`_ for developing your contribution.
See for example `this guide`_ for details on how to create one.

.. _virtual environment: https://docs.python.org/3/glossary.html#term-virtual-environment
.. _this guide: https://docs.python.org/3/library/venv.html#module-venv

To install the latest commit of the package along with the tools for testing
and linting, clone the package and run the following command in the package
root directory:

.. code-block:: bash

   python -m pip install -e . --group dev

This needs pip 25.1 or newer. With `uv`_, run ``uv pip install -e . --group
dev`` instead. Add ``--group docs`` to also install the tools for building
the documentation, and ``.[examples]`` instead of ``.`` for the packages used
in the examples.

If you don’t have `pip`_ installed, the `Python installation guide`_,
respectively the `Python Packaging User Guide`_ can guide you through the
process.

.. _uv: https://docs.astral.sh/uv/
.. _pip: https://pip.pypa.io/en/stable/
.. _Python installation guide: https://docs.python-guide.org/starting/installation/
.. _Python Packaging User Guide: https://packaging.python.org/en/latest/tutorials/installing-packages/


Running the tests
^^^^^^^^^^^^^^^^^

Run the test suite from the package root directory:

.. code-block:: bash

   python -m pytest


Building the documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Build the documentation from the package root directory, and serve it
locally to view it in a browser:

.. code-block:: bash

   sphinx-build -b html doc doc/_build/html -d doc/_build/doctree
   python -m http.server -d doc/_build/html
