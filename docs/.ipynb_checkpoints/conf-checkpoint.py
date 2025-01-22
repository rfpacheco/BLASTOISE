# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

project = 'BLASTOISE'
copyright = '2023, Ronny Pacheco'
author = 'Ronny Pacheco'
release = '0.1'

# Here it's the location path of my program
sys.path.insert(0, os.path.abspath("../blastoise"))

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# These extensions are used for the ..audomodule:: and such
# may be 'sphinx-rtd-theme' as well as 'cloud_sptheme' and 'raku-rtd_theme'
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.napoleon',
]

autodoc_mock_imports = ["biopython", "Bio"]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', ".env/*",
                    ".ipynb_checkpoints/*"]



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Changed theme to "sphinx_rtd_theme"; install with "pip install sphinx_rtd_theme"

# html_theme = 'sphinx_rtd_theme' #  For the moment, because sphinx-autobuild doesn't support it
# html_theme = 'cloud' #  For the moment, because sphinx-autobuild doesn't support it
# html_theme = 'raku' #  For the moment, because sphinx-autobuild doesn't support it
html_theme = 'sphinx_rtd_theme'
html_theme_path = ['_templates'] # This is why it didn't work before
html_static_path = ['_static']

