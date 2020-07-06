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
              'sphinx.ext.intersphinx',
              'numpydoc',
              'sphinx_rtd_theme',]
autodoc_default_options = {
     'members': True,
     'undoc-members': True,
 }
numpydoc_show_inherited_class_members = True
napoleon_include_private_with_doc = True
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('http://docs.scipy.org/doc/numpy/', None),
    'matplotlib': ('http://matplotlib.sourceforge.net/', None)
}
autosummary_generate = False
autosummary_imported_members = False

#templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
github_url = 'https://github.com/n-claes/legolas/tree/master/post_processing/pylbo'
#html_static_path = ['_static']
