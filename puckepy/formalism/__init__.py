"""
puckepy.formalism
=================

Read in pdb or xyz format files into the buffer.
Calculate the desired puckering coordinate for the prompted molecular system
"""

# To get sub modules
#from . import formalism
from .formalism import *

#__all__ = formalism.__all__.copy()
__all__ = ["Pdb", "Xyz"]

