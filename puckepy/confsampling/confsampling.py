from puckepy import puckepy


__all__ = ["Fivering","Sixring", "Peptide", "FiveringAxes","SixringAxes", "PeptideAxes"]


class Fivering :

    """ Construct a set of torsions for sampling fivering space. 

        The `interval` parameter uses linear_space() function to calculate 
        the returned parameters. This the landscape itself is 2D, this would amount in 
        `interval * interval` of pairs of restraints.

        The extent of the range is : [-60, 60, `interval`]
        ---------------
        self.nu1 : list[float]
        self.nu3 : list[float]

        >>> fivering = Fivering(21) # Every 6 degrees
        >>> for nu1, nu3 in zip(fivering.nu1, fivering.nu3)
        >>>     print(nu1, nu3)
    """
    nu1 : list[float]
    nu3 : list[float]

    def __new__(cls, interval: int) :

        return puckepy.confsampling.Fivering(interval)

    
class Peptide :
    """ Construct a set torsions for sampling peptide space 

        The `interval` parameter uses linear_space() function to calculate 
        the returned parameters. This the landscape itself is 2D, this would amount in 
        `interval * interval` of pairs of restraints.

        The extent of the range is : [0, 360, `interval`]
        ---------------
        self.phi : list[float]
        self.psi : list[float]

        >>> peptide = Peptide(37) # Every 10 degrees
        >>> for phi, psi in zip(peptide.phi, peptide.psi)
        >>>     print(phi, psi)
    """
    phi : list[float]
    psi : list[float]

    def __new__(cls, interval: int) : 
        return puckepy.confsampling.Peptide(interval)

class Sixring :
    """ Construct a set of torsions for sampling sixring space 

        The `amount` parameter will be used to cover the surface of the 
        Cremer-Pople globe with points and approximate the `amount` to an evenly
        distributed set of points.
        ---------------
        self.alpha1 : list[float]
        self.alpha2 : list[float]
        self.alpha3 : list[float]

        >>> sixring = Sixring(631) # Generate 630 points
        >>> for a1, a2, a3 in zip(sixring.alpha1,sixring.alpha2, sixring.alpha3)
        >>>     print(a1, a2, a3)
    """
    alpha1 : list[float]
    alpha2 : list[float]
    alpha3 : list[float]

    def __new__(cls, amount: int) :
        return puckepy.confsampling.Sixring(amount)

class FiveringAxes :

    """ Construct a set of axes for mapping fivering space. 

        The `interval` parameter uses linear_space() function to calculate 
        the returned parameters. This the landscape itself is 2D, this would amount in 
        `interval * interval` of pairs of restraints.

        The extent of the range is : [-60, 60, `interval`]
        ---------------
        self.zx : list[float]
        self.zy : list[float]

        >>> fivering_axes = FiveringAxes(21) # Every 6 degrees
        >>> for zx, zy in zip(fivering_axes.zx, fivering_axes.zy)
        >>>     print(zx, zy)
    """
    zx : list[float]
    zy : list[float]

    def __new__(cls, interval: int) :

        return puckepy.confsampling.FiveringAxes(interval)

    
class PeptideAxes :
    """ Construct a set of axes for mapping peptide space 

        The `interval` parameter uses linear_space() function to calculate 
        the returned parameters. This the landscape itself is 2D, this would amount in 
        `interval * interval` of pairs of restraints.

        The extent of the range is : [0, 360, `interval`]
        ---------------
        self.x : list[float]
        self.y : list[float]

        >>> peptide_axes = PeptideAxes(37) # Every 10 degrees
        >>> for x, y in zip(peptide_axes.x, peptide_axes.y)
        >>>     print(x, y)
    """
    x : list[float]
    y : list[float]

    def __new__(cls, interval: int) : 
        return puckepy.confsampling.PeptideAxes(interval)

class SixringAxes :
    """ Construct a set of axes for mapping sixring space 

        The `amount` parameter will be used to cover the surface of the 
        Cremer-Pople globe with points and approximate the `amount` to an evenly
        distributed set of points.
        ---------------
        self.rho : int
        self.theta : list[float]
        self.phi : list[float]

        >>> sixring_axes = SixringAxes(631) # Generate 630 points
        >>> for theta, phi in zip(sixring_axes.theta, sixring_axes.phi)
        >>>     print(sixring_axes.rho, theta, phi)
    """
    rho : int
    theta : list[float]
    phi : list[float]

    def __new__(cls, amount: int) :
        return puckepy.confsampling.SixringAxes(amount)

