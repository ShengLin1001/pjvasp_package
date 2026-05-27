# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

sys.path.insert(0, os.path.abspath('../..'))

project = 'mymetal'
copyright = '2024-2026, J. Pei'
author = 'J. Pei'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'myst_parser',
    'autoapi.extension',
]

templates_path = ['_templates']
exclude_patterns = [
    'modules.rst',
    'mymetal.rst',
    'mymetal.calculate*.rst',
    'mymetal.io*.rst',
    'mymetal.ml.rst',
    'mymetal.post.rst',
    'mymetal.universial*.rst',
]

# Keep the documentation build independent of optional simulation runtimes.
autodoc_mock_imports = [
    'brokenaxes',
    'hetbuilder',
    'monty',
    'myvasp',
    'ovito',
    'palettable',
    'prettytable',
    'pymatgen',
    'torch',
    'torchvision',
]
autodoc_member_order = 'bysource'
autodoc_typehints = 'description'

autoapi_type = 'python'
autoapi_dirs = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), '../../mymetal'))
]
autoapi_root = 'autoapi'
autoapi_options = [
    'members',
    'undoc-members',
    'show-inheritance',
    'show-module-summary',
]



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_title = 'mymetal Documentation'
