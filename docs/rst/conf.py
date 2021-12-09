# Configuration file for the Sphinx documentation builder.
# For a full documentation: https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
project = 'QSF'
copyright = '2021, Michal Mandrysz'
author = 'Michal Mandrysz'
# The full version, including alpha/beta/rc tags
release = '0.5'

import subprocess, os
import sphinx_rtd_theme

# Check if we're running on Read the Docs' servers
on_rtd = os.environ.get('READTHEDOCS', None) == 'True'
if on_rtd:
    subprocess.call('cd ..; doxygen', shell=True)
    
# -- General configuration ---------------------------------------------------
extensions = [ "breathe" ]
breathe_default_project = project
breathe_projects={'QSF': '../xml'}
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
htmlhelp_basename = 'qsfdoc'
html_static_path = ['_static']
html_theme = 'sphinx_rtd_theme'
# html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
