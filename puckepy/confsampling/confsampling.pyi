
__all__: list[str]


class Fivering :
    nu1: list[float]
    nu3: list[float]

    def __new__(cls, interval: int) -> Fivering : ...

class Peptide :
    phi: list[float]
    psi: list[float]

    def __new__(cls, interval: int) -> Peptide : ...

class Sixring :
    alpha1: list[float]
    alpha2: list[float]
    alpha3: list[float]

    def __new__(cls, amount: int) -> Sixring : ...
