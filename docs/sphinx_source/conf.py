# Configuration file for the Sphinx documentation builder.

# -- Path setup --------------------------------------------------------------
import os
import sys
sys.path.insert(0, os.path.abspath('../../post-processing/pylbo'))


# -- Project information -----------------------------------------------------
project = 'Pylbo'
copyright = 'Niels Claes, Jordi De Jonghe and Rony Keppens (2020)'
author = 'Niels Claes, Jordi De Jonghe and Rony Keppens'
release = '1.0'


# -- General configuration ---------------------------------------------------
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.autosummary',
              'sphinx.ext.mathjax',
              'sphinx.ext.viewcode',
              'sphinx.ext.napoleon',
              'numpydoc']
napoleon_include_private_with_doc = True
autosummary_generate = True

#templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
html_theme = 'classic'
#html_static_path = ['_static']
