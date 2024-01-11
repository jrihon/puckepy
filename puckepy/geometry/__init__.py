"""
puckepy.geometry
=================

Simple functions to calculate geometric data on coordinate systems

There is a specific type throughout the module : `Coordinates3D`
This is just an array of `length == 3`, of floating point numbers, 
representing a 3D coordinate.
>>> from typing import TypeAlias
>>> Coordinates3D: TypeAlias = tuple[float, float, float]
"""

# To get to sub modules
from . import geometry
from .geometry import *

# this .__all__ method is created in the geometry __init__.pyi and geometry.pyi
__all__ = geometry.__all__.copy() 
#__all__ = ["dihedral", "bondangle", "bondlength"] 
