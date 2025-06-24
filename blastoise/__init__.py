"""
BLASTOISE: A Tool for Repetitive Sequence Discovery in Genomic Data
====================================================================

Author: R. Pacheco
Version: 0.4.2
License: MIT
"""

__version__ = "0.4.2"
__author__ = "R. Pacheco"

# Make the main function available at the package level
from .main import main

__all__ = ["main", "modules", "extra"]