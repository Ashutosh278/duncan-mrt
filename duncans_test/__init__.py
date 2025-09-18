

"""
duncans_test

A Python package for performing Duncan's Multiple Range Test,
including statistical analysis and publication-ready plots.
"""

# Import the core Duncan's test logic and result class.
from .duncans import duncan_test, DuncanResults

# Import the plotting functions, now from the dedicated 'plots' module.
from .plotting import plot_bar, plot_cld, plot_heatmap

# This controls what is imported when a user runs 'from duncans_test import *'.
__all__ = [
    'duncan_test',
    'DuncanResults',
    'plot_bar',
    'plot_cld',
    'plot_heatmap',
]

