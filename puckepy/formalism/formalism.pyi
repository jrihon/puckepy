#from typing import TypeAlias 
#Coordinate3D: TypeAlias = list[float]

__all__: list[str]

class Pdb:
    coordinates: list[list[float]]
    atomnames: list[str]

    def __new__(cls, filename: str) -> Pdb : ...
    def parse(self) : ...
    """ Is parse on Pdb being recognised?
    """

class Xyz:

    def __new__(cls, filename: str) -> Xyz : ...
    def parse(self) -> list[list[float]] : ...
    """ Is parse on Xyz being recognised?
    """
