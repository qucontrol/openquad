.. OpenQuad - open database for multi-dimensional numerical integration
.. SPDX-FileCopyrightText: 2026  Alexander Blech
.. SPDX-License-Identifier: MPL-2.0

Making a release
----------------

Releases are made by the maintainers. Version numbers follow `semantic
versioning`_ and the `PyPA version specifiers`_. Between releases, the version
in ``src/openquad/__init__.py`` carries a ``+dev`` suffix, as in
``0.3.0+dev``. Below, ``X.Y.Z`` stands for the new version.

.. _semantic versioning: https://semver.org/
.. _PyPA version specifiers: https://packaging.python.org/en/latest/specifications/version-specifiers/


Before the release
^^^^^^^^^^^^^^^^^^

- CI is green on ``main``.
- The supported versions are up to date. They follow `SPEC 0`_: check the
  lower bounds of ``requires-python``, ``numpy`` and ``scipy`` and the Python
  classifiers in ``pyproject.toml`` against the SPEC 0 support table, and drop
  versions whose support window has ended.

.. _SPEC 0: https://scientific-python.org/specs/spec-0000/


Steps
^^^^^

1. **Set the version.** In ``src/openquad/__init__.py``, change ``__version__``
   from ``X.Y.Z+dev`` to the release version ``X.Y.Z``. Commit with the message
   ``Set version to X.Y.Z`` and merge it into ``main`` through a pull request.

2. **Tag the release** on that commit and push the tag:

   .. code-block:: bash

      git tag vX.Y.Z <commit>
      git push origin vX.Y.Z

3. **Build and upload** the package from the tagged commit, with a clean
   working tree. Uploading needs a PyPI API token for the project.

   .. code-block:: bash

      python -m pip install build twine
      python -m build
      python -m twine check dist/*
      python -m twine upload dist/*

   ``python -m build`` writes the source distribution and the wheel to
   ``dist/``. Remove older files from ``dist/`` first, so that only the new
   version is uploaded.

4. **Create a GitHub release** for the tag ``vX.Y.Z``, with a short list of
   the changes that refers to the issues and pull requests involved.

5. **Publish the documentation** on the ``gh-pages`` branch. Each release has
   its own folder, ``latest`` is a link to the newest one, and ``index.html``
   forwards to ``latest``:

   .. code-block:: bash

      python -m sphinx -b html doc build/docs
      git worktree add ../openquad-gh-pages gh-pages
      cp -r build/docs ../openquad-gh-pages/vX.Y.Z
      cd ../openquad-gh-pages
      ln -sfn vX.Y.Z latest

   On systems without symbolic links, git checks ``latest`` out as a small
   text file holding the folder name: replace that name instead. Update the
   version in the fallback link in ``index.html``, then commit with the message
   ``Add documentation for vX.Y.Z`` and push.

6. **Start the next development cycle.** Change ``__version__`` to
   ``X.Y.Z+dev``, commit with the message ``Bump version to vX.Y.Z+dev``, and
   merge it through a pull request.


Planned changes
^^^^^^^^^^^^^^^

Uploading is planned to move to a GitHub Actions workflow with `PyPI Trusted
Publishing`_, which needs no API token. It will replace step 3.

.. _PyPI Trusted Publishing: https://docs.pypi.org/trusted-publishers/
