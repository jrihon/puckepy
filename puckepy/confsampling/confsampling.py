from puckepy import puckepy


__all__ = ["Fivering","Sixring", "Peptide"]


class Fivering :

    """ Construct a set of axis and torsions for sampling 
        fivering space. 

        The `interval` parameter uses linear_space() function to calculate 
        the returned parameters. This the landscape itself is 2D, this would amount in 
        `interval * interval` of pairs of restraints.

        The extent of the range is : [-60, 60, `interval`]
        ---------------
        self.nu1 : List[float]
        self.nu3 : List[float]

        >>> fivering = Fivering(21) # Every 6 degrees
        >>> for nu1, nu3 in zip(fivering.nu1, fivering.nu3)
        >>>     print(nu1, nu3)
    """
    nu1 : list[float]
    nu3 : list[float]

    def __new__(cls, interval: int) :

        return puckepy.confsampling.Fivering(interval)

    
class Peptide :
    """ Construct a set of axis and torsions for sampling 
        peptide space 

        The `interval` parameter uses linear_space() function to calculate 
        the returned parameters. This the landscape itself is 2D, this would amount in 
        `interval * interval` of pairs of restraints.

        The extent of the range is : [0, 360, `interval`]
        ---------------
        self.phi : List[float]
        self.psi : List[float]

        >>> peptide = peptide(37) # Every 10 degrees
        >>> for phi, psi in zip(peptide.phi, peptide.psi)
        >>>     print(nu1, nu3)
    """
    phi : list[float]
    psi : list[float]

    def __new__(cls, interval: int) : 
        return puckepy.confsampling.Peptide(interval)

class Sixring :
    """ Construct a set of axis and torsions for sampling 
        sixring space 

        The `amount` parameter will be used to cover the surface of the 
        Cremer-Pople globe with points and approximate the `amount` to an evenly
        distributed set of points.
        ---------------
        self.alpha1 : List[float]
        self.alpha2 : List[float]
        self.alpha3 : List[float]

        >>> sixring = sixring(631) # Generate 630 points
        >>> for a1, a2, a3 in zip(sixring.alpha1,sixring.alpha2, sixring.alpha3)
        >>>     print(a1, a2, a3)
    """
    alpha1 : list[float]
    alpha2 : list[float]
    alpha3 : list[float]

    def __new__(cls, amount: int) :
        return puckepy.confsampling.Sixring(amount)
