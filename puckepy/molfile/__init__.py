
"""
Read in pdb or xyz format files into the buffer.
"""

class FromPdb:
    """ Pdb Class """
    def __init__(self, filename):
        """ Pdb Class constructor.
        ==
        @attribute atom_names: List[str],
        @attribute coordinates: List[3DCoordinates] 
            where 3DCoordinates : [float, float, float]
        ---------------
        Read in a molecule file with the '.pdb' file extension.
            Should be pdb-formatted 
        >>> Pdb(filename)
        """


class FromXyz:
    """ Xyz Class """
    def __init__(self, filename):
        """ Xyz Class constructor.
        ==
        @attribute coordinates: List[3DCoordinates] 
            where 3DCoordinates : [float, float, float]
        ---------------
        Read in a molecule file with the '.xyz' file extension.
            Should be xyz-formatted 

        >>> Xyz(filename)
        """




