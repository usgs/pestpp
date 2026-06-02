import os
import sys

# -- Project information -----------------------------------------------------
project = "PEST++"
copyright = "2024, PEST++ Development Team"
author = "PEST++ Development Team"
release = "5"

# -- General configuration ---------------------------------------------------
extensions = [
    "breathe",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "myst_parser",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "_doxygen", "Thumbs.db", ".DS_Store"]

# -- Breathe configuration ---------------------------------------------------
# Path is relative to conf.py (i.e. docs/_doxygen/xml)
breathe_projects = {
    "pestpp": os.path.join(os.path.dirname(__file__), "_doxygen", "xml"),
}
breathe_default_project = "pestpp"
breathe_default_members = ("members", "undoc-members")

# -- HTML output -------------------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_logo = None

# Copy the users manual media folder into the HTML output so that
# image references in pestpp_users_manual.md resolve correctly.
# The pre-build step (or local setup) copies documentation/media/ -> docs/media/
# and Sphinx places it at _build/html/media/.
# _manual_assets/media/ is copied to _build/html/media/ so that
# the ./media/imageN.png references in the users manual resolve correctly.
html_extra_path = ["_manual_assets"]

# -- MyST (Markdown support) -------------------------------------------------
myst_enable_extensions = ["colon_fence", "html_image"]
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
