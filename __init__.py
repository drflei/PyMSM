"""PyMSM - Python implementation for calculating geomagnetic rigidity cutoff.

This package provides tools for calculating the geomagnetic rigidity cutoff
and transmission functions for particles in the Earth's magnetosphere.
"""

__version__ = '0.1.1'
__author__ = 'drflei'
__license__ = 'LGPL'

from .pymsm import PyMSM, MapDB

__all__ = ['PyMSM', 'MapDB']