from typing import TypeAlias
Coordinates3D: TypeAlias = tuple[float, float, float]

__all__: list[str]

def dihedral(p0: Coordinates3D, p1: Coordinates3D, p2: Coordinates3D, p3: Coordinates3D) -> float : ...

def bondangle(p0: Coordinates3D, p1: Coordinates3D, p2: Coordinates3D) -> float : ...

def bondlength(p0: Coordinates3D, p1: Coordinates3D) -> float : ...
