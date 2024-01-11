from puckepy import puckepy

__all__ = ["Pdb", "Xyz", "CP5", "CP6","AS", "SP"]   # Classes
__all__.extend(["write_to_pdb", "write_to_xyz"])    # Function

class Pdb:

    atom_names: list[str]
    coordinates: list[list[float]] 

    def __new__(cls, filename: str) :
        """ Pdb Class constructor.
        ==

        Reads from a pdb-formatted file. Suited for single molecule files.
        ATOM      2  C6'  A     10      24.803  51.735  41.199  1.00  0.00           C  
        ATOM      5  C5'  A     10      25.097  52.567  42.397  1.00  0.00           C  
        ---------------
        where Coordinates3D : [float, float, float]
        
        >>> Pdb(filename)
        """
        return puckepy.formalism.Pdb(filename)

    def parse(self) -> None : 
        """ Parses the queried Pdb file for the 
            - atomnames 
            - coordinates(x, y, z)
            This populates the attributes.

            The function mutates the object in place
            >>> pdb = Pdb("molecule.pdb")
            >>> pdb.parse()
        """
        self.parse()

    def parse_by_monomers(self) -> list[puckepy.formalism.Pdb] : 
        """ Parses the queried Pdb file for the 
            - atomnames 
            - coordinates(x, y, z)
            This populates the attributes.

            The function mutates the object in place
            >>> pdb = Pdb("molecule.pdb")
            >>> listPdbs = pdb.parse_by_monomers()
        """
        return self.parse_by_monomers()

class Xyz:

    def __new__(cls, filename: str):
        """ Xyz Class constructor.
        ==

        Reads from a xyz-formatted file. Suited for single molecule files.
        O   3.76770440038636      1.71999235396699      1.14581624607411
        C   2.53548022010070      2.32709191442346      0.78140278302649
        ---------------
        where 3DCoordinates : [float, float, float]

        >>> Xyz(filename)
        """
        return puckepy.formalism.Xyz(filename)

    def parse(self) -> list[list[float]]:
        """ Parses the queried Xyz file for the 
            - coordinates 
            @returns list[list[float]]
        """ 
        return self.parse()


class CP5:

    amplitude: float
    phase_angle: float

    def __new__(cls, amplitude: float = 0. , phase_angle: float = 0. ):
        """ Cremer-Pople Class constructor for Fivering systems.
        ==
        Create a class to manipulate Cremer-Pople coordinates.
        ---------------
        @attribute amplitude: float [ 0. <= amplitude <= 1. ]
        @attribute phase_angle: float [ 0. <= phase_angle <= 360. ]

        >>> CP(0.35, 288.) # => C3' Endo
        """
        return puckepy.formalism.CP5(amplitude, phase_angle)

    def from_atomnames(self, pdb: Pdb, query_names: list[str]) -> tuple[float, float] : 
        """ Get Cremer-Pople coordinates by querying from the atom names of the prompted Pdb(). 
        """

        return self.from_atomnames(pdb, query_names)

    def from_indices(self, coordinates: list[list[float]], indices: list[int])  -> tuple[float, float] : 
        """ Get Cremer-Pople coordinates by querying from the indices of the prompted coordinates. 
        """

        return self.from_indices(coordinates, indices)

    def invert(self) -> list[list[float]]:
        """ Perform an inversion of the Cremer-Pople coordinates and get returned the 
            molecular conformation the five-membered ring results in.
        """
        return self.invert()


class CP6:

    amplitude: float
    phase_angle: float
    theta: float

    def __new__(cls, amplitude: float = 0. , phase_angle: float = 0. , theta: float = 0.):
        """ Cremer-Pople Class constructor for Fivering systems.
        ==
        Create a class to manipulate Cremer-Pople coordinates.
        ---------------
        @attribute amplitude: float [ 0. <= amplitude <= 1. ]
        @attribute phase_angle: float [ 0. <= phase_angle <= 360. ]
        @attribute theta: float [ 0. <= theta <= 180. ]

        >>> CP6(0.35, 90., 90.) # => (O5', C3')^Boat
        """
        return puckepy.formalism.CP6(amplitude, phase_angle, theta)

    def from_atomnames(self, pdb: Pdb, query_names: list[str]) -> tuple[float, float, float] : 
        """ Get Cremer-Pople coordinates by querying from the atom names of the prompted Pdb(). 
        """

        return self.from_atomnames(pdb, query_names)

    def from_indices(self, coordinates: list[list[float]], indices: list[int])  -> tuple[float, float, float] : 
        """ Get Cremer-Pople coordinates by querying from the indices of the prompted coordinates. 
        """

        return self.from_indices(coordinates, indices)

    def invert(self) -> list[list[float]]:
        """ Perform an inversion of the Cremer-Pople coordinates and get returned the 
            molecular conformation the six-membered ring results in.
        """
        return self.invert()

class AS:

    amplitude: float
    phase_angle: float

    def __new__(cls, amplitude: float = 0. , phase_angle: float = 0. ):
        """ Altona-Sundaralingam Class constructor for Fivering systems.
        ==
        Create a class to manipulate Altona-Sundaralingam coordinates.
        ---------------
        @attribute amplitude: float [ 0. <= amplitude <= 1. ]
        @attribute phase_angle: float [ 0. <= phase_angle <= 360. ]

        >>> CP(0.35, 288.) # => C3' Endo
        """
        return puckepy.formalism.AS(amplitude, phase_angle)

    def from_atomnames(self, pdb: Pdb, query_names: list[str]) -> tuple[float, float] : 
        """ Get Altona-Sundaralingam coordinates by querying from the atom names of the prompted Pdb(). 
        """

        return self.from_atomnames(pdb, query_names)

    def from_indices(self, coordinates: list[list[float]], indices: list[int])  -> tuple[float, float] : 
        """ Get Altona-Sundaralingam coordinates by querying from the indices of the prompted coordinates. 
        """

        return self.from_indices(coordinates, indices)

#    def invert(self) -> list[list[float]]:
#        """ Converts to AS coordinates into CP5 and then performs
#            an inversion of the Cremer-Pople coordinates and get returned
#            the molecular conformation the five-membered ring results in.
#        """
#        return self.invert()


class SP :

    def __new__(cls) : 
        """ Strauss-Pickett Class constructor for Sixring systems.
        ==
        Create a class to calculate Strauss-Pickett coordinates.
        ---------------
        """ 
        return puckepy.formalism.SP()

    def from_atomnames(self, pdb: Pdb, query_names: list[str]) -> tuple[list[float], list[float]]  : 
        """ Get Strauss-Pickett coordinates by querying from the atom names of the prompted Pdb(). 
        """

        return self.from_atomnames(pdb, query_names)

    def from_indices(self, coordinates: list[list[float]], indices: list[int])  -> tuple[list[float], list[float]]  : 
        """ Get Strauss-Pickett coordinates by querying from the indices of the prompted coordinates. 
        """

        return self.from_indices(coordinates, indices)



def write_to_pdb(coordinates: list[list[float]], atomnames: list[str]) -> None :
    """ Write a set of coordinates with their respective atomnames to a `.pdb` 
        formatted file. 
    """

    puckepy.formalism.write_to_pdb(coordinates, atomnames)


def write_to_xyz(coordinates: list[list[float]]) -> None :
    """ Write a set of coordinates to an `.xyz` formatted file. 
    """

    puckepy.formalism.write_to_xyz(coordinates)
