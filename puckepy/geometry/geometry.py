from puckepy import puckepy # this imports the puckepy.abi3.so binary

from typing import TypeAlias
Coordinates3D: TypeAlias = tuple[float, float, float]

__all__ = ["dihedral", "bondangle", "bondlength"] 


def dihedral(p0: Coordinates3D, p1: Coordinates3D, p2: Coordinates3D, p3: Coordinates3D) -> float :
    """ Calculate the dihedral between four coordinate points
        --------------------
        A simple query can be made to calculate explicit coordinates :
        >>> dihedral = dihedral([ [6.105, 8.289, 4.633],
        >>>                       [7.360, 7.768, 4.827],
        >>>                       [7.390, 6.603, 5.551],
        >>>                       [6.301, 5.942, 6.079] ])

        One can also index the Pdb() class' coordinate attribute : 
        >>> pdb = Pdb("nucleoside.pdb")
        >>> dihedral = dihedral([ pdb.coordinate[1]],
        >>>                       pdb.coordinate[4]],
        >>>                       pdb.coordinate[2]],
        >>>                       pdb.coordinate[7]] ])
    """
    return puckepy.geometry.dihedral(p0, p1, p2, p3)

def bondangle(p0: Coordinates3D, p1: Coordinates3D, p2: Coordinates3D)  :
    """ Calculate the bondangle between three coordinate points
        --------------------
        A simple query can be made to calculate explicit coordinates :
        >>> bondangle = bondangle([ [6.105, 8.289, 4.633],
        >>>                         [7.360, 7.768, 4.827],
        >>>                         [6.301, 5.942, 6.079] ])

        One can also index the Pdb() class' coordinate attribute : 
        >>> pdb = Pdb("nucleoside.pdb")
        >>> bondangle = bondangle([ pdb.coordinate[1]],
        >>>                         pdb.coordinate[4]],
        >>>                         pdb.coordinate[7]] ])

    """
    return puckepy.geometry.bondangle(p0, p1, p2)

def bondlength(p0: Coordinates3D, p1: Coordinates3D) -> float :
    """ Calculate the bondlength between two coordinate points
        --------------------
        A simple query can be made to calculate explicit coordinates :
        >>> bondlength = bondlength([ [6.105, 8.289, 4.633],
        >>>                           [6.301, 5.942, 6.079] ])

        One can also index the Pdb() class' coordinate attribute : 
        >>> pdb = Pdb("nucleoside.pdb")
        >>> bondlength = bondlength([ pdb.coordinate[1]],
        >>>                           pdb.coordinate[7]] ])

    """
    return puckepy.geometry.bondlength(p0, p1)
