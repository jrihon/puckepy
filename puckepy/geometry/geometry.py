from puckepy import puckepy # this imports the puckepy.abi3.so binary




__all__ = ["dihedral", "bondangle", "bondlength"] 
#
#
def dihedral(p0: list[float], p1: list[float], p2: list[float], p3: list[float]) -> float :
    """ Calculate the dihedral between four coordinate points

    """
    return puckepy.geometry.dihedral(p0, p1, p2, p3)

def bondangle(p0: list[float], p1: list[float], p2: list[float])  :
    """ Calculate the bondangle between three coordinate points

    """
    return puckepy.geometry.bondangle(p0, p1, p2)

def bondlength(p0: list[float], p1: list[float]) -> float :
    """ Calculate the bondlength between two coordinate points

    """
    return puckepy.geometry.bondlength(p0, p1)
