"""
puckepy.formalism
=================

Read in pdb or xyz format files into the buffer.
Calculate the desired puckering coordinate for the prompted molecular system

There is a specific type throughout the module : `Coordinates3D`
This is just an array of `length == 3`, of floating point numbers, 
representing a 3D coordinate.
>>> from typing import TypeAlias
>>> Coordinates3D: TypeAlias = tuple[float, float, float]
"""

# To get sub modules
from . import formalism
from .formalism import *

__all__ = formalism.__all__.copy()

