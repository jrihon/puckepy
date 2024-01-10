from typing import List

class Fivering :
    """ Construct a set of axis and torsions for sampling 
        fivering space 
        ---------------
        @params: interval -> int

        @output : Fivering
        @attribute : nu1 : List[float]
        @attribute : nu3 : List[float]
    """
    nu1 : List[float]
    nu3 : List[float]
    def __new__(cls, interval: int) -> None :  ...
#        nu3 : List[float]
    
class Peptide :
    """ Construct a set of axis and torsions for sampling 
        peptide space 
        ---------------
        @params: interval -> int

        @output : Peptide
        @attribute : phi : List[float]
        @attribute : psi : List[float]
    """
    phi : List[float]
    psi : List[float]
    def __new__(cls, interval: int) -> None :  ...

class Sixring :
    """ Construct a set of axis and torsions for sampling 
        sixring space 
        ---------------
        @params: amount -> int

        @output : Sixring
        @attribute : alpha1 : List[float]
        @attribute : alpha2 : List[float]
        @attribute : alpha3 : List[float]
    """
    alpha1 : List[float]
    alpha2 : List[float]
    alpha3 : List[float]
    def __new__(cls, amount: int) -> None : ...
