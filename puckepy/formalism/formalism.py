from puckepy import puckepy

from typing import TypeAlias
Coordinates3D: TypeAlias = tuple[float, float, float]

__all__ = ["Pdb", "Xyz", "CP5", "CP6","AS", "SP"]   # Classes
__all__.extend(["write_to_pdb", "write_to_xyz"])    # Function

class Pdb:

    atom_names: list[str]
    coordinates: list[list[float]] 

    def __new__(cls, filename: str) :
        """ Pdb Class constructor.

            Reads from a pdb-formatted file. Suited for single molecule files.
            ATOM      2  C6'  A     10      24.803  51.735  41.199  1.00  0.00           C  
            ATOM      5  C5'  A     10      25.097  52.567  42.397  1.00  0.00           C  
            ---------------
            self.atom_names : list[float]
            self.coordinates : list[list[float]]
            
            >>> pdb = Pdb(filename)
        """
        return puckepy.formalism.Pdb(filename)

    def parse(self) -> puckepy.formalism.Pdb : 
        """ Parses the queried `monomer` Pdb file to populates the attributes.

            The method mutates the object in place

            >>> pdb = Pdb("molecule.pdb")
            >>> pdbContent = pdb.parse()
        """
        return self.parse()

    def parse_by_monomers(self) -> list[puckepy.formalism.Pdb] : 
        """ Parses the queried `duplex` Pdb file to populates the attributes.
            The method mutates the object in place

            Returns a list[] of Pdb() objects

            >>> pdbs = Pdb("duplex.pdb")
            >>> listPdbs = pdbs.parse_by_monomers()
            >>> [pdb.parse() for pdb in pdbs] # populate the respective fields of the Pdb() objects
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
            This class does not have public attributes

            >>> xyz = Xyz(filename)
        """
        return puckepy.formalism.Xyz(filename)

    def parse(self) -> list[Coordinates3D]:
        """ Parses the queried Xyz file to return the coordinates 
            of the `xyz` file

            Returns the coordinates of the queried `xyz` file

            >>> xyz = Xyz(filename)
            >>> xyzCoordinates = xyz.parse()
        """ 
        return self.parse()


class CP5:

    def __new__(cls, amplitude: float = 0. , phase_angle: float = 0. ):
        """ Cremer-Pople Class constructor for five-membered ring systems.
            ==
            Create a class to manipulate Cremer-Pople coordinates.
            ---------------
            amplitude: float [ 0. <= amplitude <= 1. ] `radians`
            phase_angle: float [ 0. <= phase_angle <= 360. ] `degrees`
            ---------------
            This class does not have public attributes

            >>> cp5 = CP5(0.35, 288.) # => C3' Endo

            >>> cp5 = CP5() # => Defaults to CP5(0., 0.)
        """
        return puckepy.formalism.CP5(amplitude, phase_angle)

    def from_atomnames(self, pdb: Pdb, query_names: list[str]) -> tuple[float, float] : 
        """ Get Cremer-Pople coordinates by querying from the atom names of the prompted Pdb(). 

            Returns the `amplitude` and `phase angle`

            >>> pdb = Pdb("adenosine.pdb")
            >>> amplitude, phaseAngle = CP5().from_atomnames(pdb=pdb, query_names=["O4'", "C1'", "C2'", "C3'", "C4'"])
        """
        return self.from_atomnames(pdb, query_names)

    def from_indices(self, coordinates: list[Coordinates3D], indices: list[int])  -> tuple[float, float] : 
        """ Get Cremer-Pople coordinates by querying from the indices of the prompted coordinates. 
            NOTE: Indexing is 0-based

            Returns the `amplitude` and `phase angle`

            >>> pdb = Pdb("adenosine.pdb")
            >>> pdbContent = pdb.parse()
            >>> amplitude, phaseAngle = CP5().from_indices(coordinates=pdbContent.coordinates, indices=[7, 8, 26, 24, 5])

            >>> xyz = Xyz("adenosine.xyz")
            >>> coordinates = xyz.parse()
            >>> amplitude, phaseAngle = CP5().from_indices(coordinates=coordinates, indices=[7, 8, 26, 24, 5])
        """
        return self.from_indices(coordinates, indices)

    def invert(self) -> list[Coordinates3D]:
        """ Perform an inversion of the Cremer-Pople coordinates and get returned the 
            molecular conformation the five-membered ring results in.

            Returns the coordinates of the queried fivering conformation

            >>> endoC3 = CP5(0.35, 288.).invert()
        """
        return self.invert()


class CP6:

    def __new__(cls, amplitude: float = 0. , phase_angle: float = 0. , theta: float = 0.):
        """ Cremer-Pople Class constructor for six-membered ring systems.
            ==
            Create a class to manipulate Cremer-Pople coordinates.
            ---------------
            amplitude: float [ 0. <= amplitude <= 1. ] `radians`
            phase_angle: float [ 0. <= phase_angle <= 360. ] `degrees`
            theta: float [ 0. <= theta <= 180. ] `degrees`
            ---------------
            This class does not have public attributes

            >>> cp6 = CP6(0.35, 90., 90.) # => (O5', C3')^Boat

            >>> cp6 = CP6() # => Defaults to CP6(0., 0., 0.)
        """
        return puckepy.formalism.CP6(amplitude, phase_angle, theta)

    def from_atomnames(self, pdb: Pdb, query_names: list[str]) -> tuple[float, float, float] : 
        """ Get Cremer-Pople coordinates by querying from the atom names of the prompted Pdb(). 

            Returns the `amplitude`, `phase angle` and `theta`

            >>> pdb = Pdb("homodna_adenosine.pdb")
            >>> amplitude, phaseAngle, theta = CP6().from_atomnames(pdb=pdb, query_names=["O5'", "C1'", "C2'", "C3'", "C4'", "C5'"])
        """
        return self.from_atomnames(pdb, query_names)

    def from_indices(self, coordinates: list[Coordinates3D], indices: list[int])  -> tuple[float, float, float] : 
        """ Get Cremer-Pople coordinates by querying from the indices of the prompted coordinates. 
            NOTE: Indexing is 0-based

            Returns the `amplitude`, `phase angle` and `theta`

            >>> pdb = Pdb("homodna_adenosine.pdb")
            >>> pdbContent = pdb.parse()
            >>> amplitude, phaseAngle, theta = CP6().from_indices(coordinates=pdbContent.coordinates, indices=[7, 8, 26, 24, 5, 6])

            >>> xyz = Xyz("homodna_adenosine.xyz")
            >>> coordinates = xyz.parse()
            >>> amplitude, phaseAngle, theta = CP6().from_indices(coordinates=coordinates, indices=[7, 8, 26, 24, 5, 6])
        """
        return self.from_indices(coordinates, indices)

    def invert(self) -> list[list[float]]:
        """ Perform an inversion of the Cremer-Pople coordinates and get returned the 
            molecular conformation the six-membered ring results in.

            Returns the coordinates of the queried sixring conformation

            >>> boat3O = CP6(0.35, 90., 90.).invert()
        """
        return self.invert()

class AS:

    def __new__(cls, amplitude: float = 0. , phase_angle: float = 0. ):
        """ Altona-Sundaralingam Class constructor for five-membered ring systems.
            ==
            Create a class to manipulate Altona-Sundaralingam coordinates.
            ---------------
            amplitude: float [ 0. <= amplitude <= 1. ] `radians`
            phase_angle: float [ 0. <= phase_angle <= 360. ] `degrees`

            >>> altsun = AS(0.35, 18.) # => C3' Endo
        """
        return puckepy.formalism.AS(amplitude, phase_angle)

    def from_atomnames(self, pdb: Pdb, query_names: list[str]) -> tuple[float, float] : 
        """ Get Altona-Sundaralingam coordinates by querying from the atom names of the prompted Pdb(). 
            NOTE: Indexing is 0-based

            Returns the `amplitude` and `phase angle`

            >>> pdb = Pdb("adenosine.pdb")
            >>> amplitude, phaseAngle = AS().from_atomnames(pdb=pdb, query_names=["O4'", "C1'", "C2'", "C3'", "C4'"])
        """

        return self.from_atomnames(pdb, query_names)

    def from_indices(self, coordinates: list[Coordinates3D], indices: list[int])  -> tuple[float, float] : 
        """ Get Altona-Sundaralingam coordinates by querying from the indices of the prompted coordinates. 

            Returns the `amplitude` and `phase angle`

            >>> pdb = Pdb("adenosine.pdb")
            >>> pdbContent = pdb.parse()
            >>> amplitude, phaseAngle = AS().from_indices(coordinates=pdbContent.coordinates, indices=[7, 8, 26, 24, 5])

            >>> xyz = Xyz("adenosine.xyz")
            >>> coordinates = xyz.parse()
            >>> amplitude, phaseAngle = AS().from_indices(coordinates=coordinates, indices=[7, 8, 26, 24, 5])
        """

        return self.from_indices(coordinates, indices)



class SP :

    def __new__(cls) : 
        """ Strauss-Pickett Class constructor for Sixring systems.
            ==
            Create a class to calculate Strauss-Pickett coordinates.
            ---------------
            This class does not have public attributes

            >>> sp = SP()
        """ 
        return puckepy.formalism.SP()

    def from_atomnames(self, pdb: Pdb, query_names: list[str]) -> tuple[tuple[float, float, float], tuple[float, float, float]]  : 
        """ Get Strauss-Pickett coordinates by querying from the atom names of the prompted Pdb(). 
            
            To keep true to the formalism, please adhere to the specified sequence unless confident.

            >>> pdb = Pdb("homodna_adenosine.pdb")
            >>> alphas, betas = SP().from_atomnames(pdb=pdb, query_names=["O5'", "C1'", "C2'", "C3'", "C4'", "C5'"])

            >>> xyz = Xyz("homodna_adenosine.xyz")
            >>> coordinates = xyz.parse()
            >>> alphas, betas = SP().from_indices(coordinates=coordinates, query_names=[7, 8, 26, 24, 5, 6])
        """
        return self.from_atomnames(pdb, query_names)

    def from_indices(self, coordinates: list[list[float]], indices: list[int])  -> tuple[tuple[float, float, float], tuple[float, float, float]]  : 
        """ Get Strauss-Pickett coordinates by querying from the indices of the prompted coordinates. 

            >>> xyz = Xyz("homodna_adenosine.xyz")
            >>> coordinates = xyz.parse()
            >>> alphas, betas = SP().from_indices(coordinates=coordinates, indices=[7, 8, 26, 24, 5, 6])
        """
        return self.from_indices(coordinates, indices)



def write_to_pdb(filename: str, coordinates: list[Coordinates3D], residuename: str) -> None :
    """ Write a set of coordinates with their respective atomnames to a `.pdb` 
        formatted file. 

        >>> conf_coordinates = CP5(0.35, 90.).invert() # 2' endo conformation
        >>> write_to_pdb(filename: "2endo.pdb", 
        >>>              coordinates=conf_coordinates,
        >>>              residuename="ABC"
        >>>              )
    """

    puckepy.formalism.write_to_pdb(filename, coordinates, residuename)


def write_to_xyz(filename: str, coordinates: list[Coordinates3D]) -> None :
    """ Write a set of coordinates to an `.xyz` formatted file. 

        >>> conf_coordinates = CP5(0.35, 90.).invert() # 2' endo conformation
        >>> write_to_pdb(filename: "2endo.xyz", 
        >>>              coordinates=conf_coordinates,
        >>>              )
    """
    puckepy.formalism.write_to_xyz(filename, coordinates)
