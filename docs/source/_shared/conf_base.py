# Shared Sphinx configuration for the sqc documentation site.
#
# This module holds every setting common to BOTH language trees (en/ and zh/).
# Each language's conf.py does `from conf_base import *` and then overrides only
# `language` and the display title. Keeping one source of truth here prevents the
# two trees from drifting apart.
from __future__ import annotations

import sys
from pathlib import Path

# --- Make `sqc` importable for autodoc ------------------------------------
# conf_base.py lives at docs/source/_shared/; the repo root is three levels up.
_REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_REPO_ROOT))

from sqc import __version__ as _sqc_version  # noqa: E402

# --- Project metadata ------------------------------------------------------
project = "sqc"
author = "Zip"
copyright = "2026, Zip"
version = _sqc_version
release = _sqc_version

# --- Extensions ------------------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",        # NumPy-style docstrings (used across sqc/)
    "sphinx.ext.intersphinx",
    # NOTE: sphinx.ext.viewcode is intentionally NOT enabled. It renders every
    # module's full source -- including not-yet-public stubs and internal
    # imports (e.g. ZCrosstalkWorkflow) that D3 keeps off the v1 public face.
    # The source remains available in the GitHub repository.
    "myst_parser",
    "nbsphinx",
]

# --- Autodoc / autosummary -------------------------------------------------
autosummary_generate = True
# Respect each subpackage's __all__ when expanding modules: names deliberately
# excluded from a package's public API (e.g. ZCrosstalkWorkflow, D3) never enter
# the generated API pages.
autosummary_ignore_module_all = False
autodoc_default_options = {
    "members": True,
    "inherited-members": False,
    "show-inheritance": True,
}
autodoc_typehints = "description"
autodoc_member_order = "bysource"
napoleon_numpy_docstring = True
napoleon_google_docstring = False

# --- nbsphinx --------------------------------------------------------------
# Never re-execute notebooks during the build: the tutorial runs heavy QuTiP
# simulations and CI must stay fast and deterministic.
nbsphinx_execute = "never"

# --- MyST ------------------------------------------------------------------
myst_enable_extensions = ["colon_fence", "deflist", "dollarmath", "amsmath"]
# NOTE: do not map ".ipynb" here -- nbsphinx registers its own source parser for
# notebooks, and forcing an entry interferes with its title extraction.
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

# --- intersphinx -----------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "qutip": ("https://qutip.readthedocs.io/en/stable/", None),
}

# --- HTML output (furo) ----------------------------------------------------
html_theme = "furo"
html_title = f"sqc {release}"
templates_path = ["../_shared/_templates"]
html_static_path = ["../_shared/_static"]
html_css_files = ["sqc.css"]

# Inject the language switcher above furo's default sidebar components.
html_sidebars = {
    "**": [
        "sidebar/brand.html",
        "sidebar/search.html",
        "sidebar/language-switch.html",
        "sidebar/scroll-start.html",
        "sidebar/navigation.html",
        "sidebar/scroll-end.html",
    ]
}

# --- Build strictness ------------------------------------------------------
nitpicky = False

# Docstrings in sqc/ use plain-text math notation (absolute-value bars |x|,
# df/dPhi, etc.) that docutils tries to parse as RST substitution references.
# Suppressing the "docutils" category keeps `sphinx-build -W` meaningful for
# structural problems we own (broken toctrees, dangling cross-references,
# missing files) without failing on these source-docstring quirks. Cleaning
# the docstrings themselves is tracked as a separate source-quality follow-up.
suppress_warnings = ["docutils"]

# Each language tree is a self-contained source directory (en/ or zh/); the
# other language's pages are simply not present, so no cross-language orphan
# warnings under `sphinx-build -W`.


# --- Hide not-yet-implemented stubs from the API reference -----------------
# Some public classes (e.g. SensingWorkflow) carry methods that only raise
# NotImplementedError -- capabilities planned for a future release (see the
# project roadmap). These are NOT part of the v1 public surface and must not
# appear in the API docs. `__all__` filtering cannot reach *methods*, so we
# skip any member whose source body raises NotImplementedError. This is
# name-agnostic: it catches current and future stubs automatically.
def _skip_notimplemented_members(app, what, name, obj, skip, options):
    if skip:
        return None
    import inspect

    try:
        source = inspect.getsource(obj)
    except (TypeError, OSError):
        return None
    if "raise NotImplementedError" in source:
        return True
    return None


def setup(app):
    app.connect("autodoc-skip-member", _skip_notimplemented_members)
