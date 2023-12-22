""" 
puckepy.confsampling
====================

Import the confsampling module in order to generate the landscape
of sample the molecular system of choice
"""
#from . import confsampling
#from confsampling import Fivering, Sixring, Peptide
from typing import List




class Fivering :
    """ Use the Fivering class. Takes in (interval) as parameter """
    nu1 : List[float]
    nu3 : List[float]
    def __init__(self, interval: int): 
        """ Construct a set of axis and torsions for sampling 
            fivering space 
            ---------------
            @params: interval -> int

            @output : Fivering
            @attribute : nu1 : List[float]
            @attribute : nu3 : List[float]
        """
    
class Peptide :
    """ Use the Peptide class. Takes in (interval) as parameter """
    phi : List[float]
    psi : List[float]
    def __init__(self, interval: int): 
        """ Construct a set of axis and torsions for sampling 
            peptide space 
            ---------------
            @params: interval -> int

            @output : Peptide
            @attribute : phi : List[float]
            @attribute : psi : List[float]
        """

class Sixring :
    """ Use the Sixring class. Takes in (amount) as parameter """
    alpha1 : List[float]
    alpha2 : List[float]
    alpha3 : List[float]
    def __init__(self, amount: int): 
        """ Construct a set of axis and torsions for sampling 
            sixring space 
            ---------------
            @params: amount -> int

            @output : Sixring
            @attribute : alpha1 : List[float]
            @attribute : alpha2 : List[float]
            @attribute : alpha3 : List[float]
        """
#class Fivering(interval: int) -> Axis, Torsions : ...
#""" Use the Fivering class. Takes in (interval) as parameter """
#    
#class Sixring(amount: int) -> Axis, Torsions : ...
#""" Use the Sixring class. Takes in (amount) as parameter """
#
#class Peptide(interval: int) -> Axis, Torsions : ...
#""" Use the Sixring class. Takes in (interval) as parameter """
