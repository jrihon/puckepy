r"""   
 This messages is the declarative statement 
 of the puckepy module when being imported
"""

def helloworld(): ...



def add_two(a: int) -> int:
    """ Add 2 to the queried value `a` """

class Number:
    """  Testing some stuff, how to return stuff """
    nummie: int

    def __init__(self, num: int) -> None: ...
    """ Insert a valid `i32` value as an integer """

    def get_number(self) -> int: ...
    """ Get the value of nummie """


#class Peptide:
#    """
#    The restraints for peptide-like structures.
#    """
#    phi: list[int]
#    psi: list[int]
#    def __init__(self, num: int) -> None: ...
#    """
#    Given a specific value, this will return a set of [phi, psi]
#    """
#
