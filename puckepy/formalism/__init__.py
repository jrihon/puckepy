
"""
puckepy.formalism
=================

Read in pdb or xyz format files into the buffer.
Calculate the desired puckering coordinate for the prompted molecular system
"""

class Pdb:
    """ Pdb Class """
    def __init__(self, filename):
        """ Pdb Class constructor.
        ==

        Reads from a pdb-formatted file. Suited for single molecule files.
        ATOM      2  C6'  A     10      24.803  51.735  41.199  1.00  0.00           C  
        ATOM      5  C5'  A     10      25.097  52.567  42.397  1.00  0.00           C  
        ---------------
        @attribute atom_names: List[str],
        @attribute coordinates: List[3DCoordinates] 
            where 3DCoordinates : [float, float, float]

        >>> Pdb(filename)
        """


class Xyz:
    """ Xyz Class """
    def __init__(self, filename):
        """ Xyz Class constructor.
        ==

        Reads from a xyz-formatted file. Suited for single molecule files.
        O   3.76770440038636      1.71999235396699      1.14581624607411
        C   2.53548022010070      2.32709191442346      0.78140278302649
        ---------------
        @attribute coordinates: List[3DCoordinates] 
            where 3DCoordinates : [float, float, float]

        >>> Xyz(filename)
        """




