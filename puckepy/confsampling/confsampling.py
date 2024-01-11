from puckepy import puckepy


__all__ = ["Fivering","Sixring", "Peptide"]


class Fivering :
    """ Construct a set of axis and torsions for sampling 
        fivering space 
        ---------------
        @params: interval -> int

        @output : Fivering
        @attribute : nu1 : List[float]
        @attribute : nu3 : List[float]
    """
    nu1 : list[float]
    nu3 : list[float]
    def __new__(cls, interval: int) :

        return puckepy.confsampling.Fivering(interval)

    
class Peptide :
    """ Construct a set of axis and torsions for sampling 
        peptide space 
        ---------------
        @params: interval -> int

        @output : Peptide
        @attribute : phi : List[float]
        @attribute : psi : List[float]
    """
    phi : list[float]
    psi : list[float]
    def __new__(cls, interval: int) : 
        return puckepy.confsampling.Peptide(interval)

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
    alpha1 : list[float]
    alpha2 : list[float]
    alpha3 : list[float]
    def __new__(cls, amount: int) :
        return puckepy.confsampling.Sixring(amount)
