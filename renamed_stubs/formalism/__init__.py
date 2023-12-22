"""
puckepy.formalism
=================

Read in pdb or xyz format files into the buffer.
Calculate the desired puckering coordinate for the prompted molecular system
"""

from typing import List, TypeAlias 
Coordinate3D: TypeAlias = list[type[float]]

class Pdb:
    """ Pdb Class """
    atom_names: List[str]
    coordinates: List[Coordinate3D] 
    def __init__(self, filename):
        """ Pdb Class constructor.
        ==

        Reads from a pdb-formatted file. Suited for single molecule files.
        ATOM      2  C6'  A     10      24.803  51.735  41.199  1.00  0.00           C  
        ATOM      5  C5'  A     10      25.097  52.567  42.397  1.00  0.00           C  
        ---------------
        where Coordinates3D : [float, float, float]
        
        >>> Pdb(filename)
        """


class Xyz:
    """ Xyz Class """
    coordinates: List[Coordinate3D] 
    def __init__(self, filename):
        """ Xyz Class constructor.
        ==

        Reads from a xyz-formatted file. Suited for single molecule files.
        O   3.76770440038636      1.71999235396699      1.14581624607411
        C   2.53548022010070      2.32709191442346      0.78140278302649
        ---------------
        where 3DCoordinates : [float, float, float]

        >>> Xyz(filename)
        """

#class CP:
#    """ Cremer-Pople Class for Fiverings """
#    amplitude: float
#    phase_angle: float
#    def __init__(self, amplitude, phase_angle):
#        """ Cremer-Pople Class constructor for Fivering systems.
#        ==
#        Create a class to manipulate Cremer-Pople coordinates.
#        ---------------
#        @attribute amplitude: float [ 0. <= amplitude <= 1. ]
#        @attribute phase_angle: float [ 0. <= phase_angle <= 360. ]
#
#        >>> CP(0.35, 288.) # => C3' Endo
#        """
#
#    def to_as(self) : ...

class AS:
    """ Altona-Sundaralingam Class for Fiverings """
    amplitude: float
    phase_angle: float
    def __init__(self, amplitude, phase_angle):
        """ Altona-Sundaralingam Class constructor for Fivering systems.
        ==
        Create a class to manipulate Altona-Sundaralingam coordinates.
        ---------------
        @attribute amplitude: float [ 0. <= amplitude <= 1. ]
        @attribute phase_angle: float [ 0. <= phase_angle <= 360. ]

        >>> AS(0.35, 18.) # => C3' Endo
        """


    def to_cp(self) : ...


