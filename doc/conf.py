import re

from intersphinx_registry import get_intersphinx_mapping

import openquad

# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'OpenQuad'
copyright = '2024, OpenQuad developers'
author = 'Alexander Blech et al.'
release = openquad.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Report broken links as warnings
nitpicky = True
# nitpick_ignore = [("py:class", "callable")]

extensions = [
    'autoapi.extension',
    'sphinx.ext.intersphinx',
    'sphinx.ext.doctest',
    'sphinx_design',
    'sphinx_copybutton',
    'numpydoc',
]

intersphinx_mapping = get_intersphinx_mapping(
packages={"python", "numpy", "scipy"}
)
#intersphinx_disabled_reftypes = ["*"]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Will be included in the beginning of every RST file
rst_prolog = """
.. |openquad| replace:: OpenQuad
"""

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'pydata_sphinx_theme'

html_theme_options = {
    "logo": {
        "alt_text": "OpenQuad - home",
        "image_light": "_static/images/logo.svg",
        "image_dark": "_static/images/logo.svg",
    },
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/qucontrol/openquad",
            "icon": "fa-brands fa-github",
        },
        {
            "name": "PyPi",
            "url": "https://pypi.org/project/openquad",
            "icon": "fa-custom fa-pypi",
        },
    ],
    "external_links": [
        {"name": "Changelog", "url": "https://github.com/qucontrol/openquad/releases"},
    ],
    "footer_start": ["copyright"],
    "footer_center": ["sphinx-version"],
    "navbar_end": [
        "theme-switcher",
        #"version-switcher",
        "navbar-icon-links",
    ],
    "use_edit_page_button": False,
    "navbar_align": "content",
    #"show_version_warning_banner": True,
    "show_toc_level": 2,
    "show_nav_level": 2, #TODO
    "navigation_depth": 4, #TODO
    #"collapse_navigation": True,
    "announcement": "https://raw.githubusercontent.com/qucontrol/openquad/main/doc/_templates/announcement_banner.html"
}


html_static_path = ['_static']
html_css_files = ['openquad.css']
html_copy_source = False
html_show_sourcelink = False

# -- Other options       ------------------------------------------------------

# sphinx-copybutton settings
copybutton_prompt_text = ">>> "

# sphinx-autoapi settings
autoapi_type = "python"
autoapi_template_dir = "_templates/autoapi"
autoapi_dirs = ["../src/openquad"]
autoapi_root = "api/autoapi"
autoapi_add_toctree_entry = False
autoapi_keep_files = False
autoapi_member_order = "groupwise"
autoapi_own_page_level = "class"
autoapi_python_class_content = "class"
autoapi_options = [
    "undoc-members",
    "show-module-summary",
    "imported-members",
    "inherited-members",
]

def skip_submodules(app, what, name, obj, skip, options):
    if what == "module":
        skip = True
    elif what == "package" and name != "openquad":
        skip = True
    return skip

# -----------------------------------------------------------------------------
def setup(sphinx):
    sphinx.connect("autoapi-skip-member", skip_submodules)
