# Documentation


## Conformational sampling

### Sampling of the *Peptide* space
```python
from puckepy.confsampling import Peptide, PeptideAxes
pepspace = Peptide(37)
for phi, psi in zip(pepspace.phi, pepspace.psi): 
    print(phi, psi)

pepaxes = PeptideAxes(37)
for x, y in zip(pepaxes.x, pepaxes.y): 
    print(x, y)
```

###  Sampling of the *Fivering* space
```python
from puckepy.confsampling import Fivering, FiveringAxes
fivespace = Fivering(21)
for nu1, nu3 in zip(fivespace.nu1, fivespace.nu3): 
    print(nu1, nu3)

fiveaxes = FiveringAxes(21)
for zx, zy in zip(fiveaxes.zx, fiveaxes.zy): 
    print(zx, zy)
```

### Sampling of the *Sixring* space
```python
from puckepy.confsampling import Sixring, SixringAxes
sixspace = Sixring(21)
for a1, a2, a3 in zip(sixspace.alpha1, sixspace.alpha2, sixspace.alpha3): 
    print(a1, a2, a3)

sixaxes = Sixring(21)
for theta, phi in zip(sixaxes.theta, sixaxes.phi): 
    print(sixaxes.rho, theta, phi)
```
<br>

## Puckering formalism

### Calculating the pucker coordinates of a *five-membered ring* by the *Cremer-Pople* formalism, from a *pdb* formatted file
```python
from puckepy.formalism import Pdb, CP5
pdb = Pdb("./fivering_3endo.pdb").parse()
amplitude, phaseangle = CP5().from_atomnames(pdb=pdb,
                                             query_names=["O4'","C1'", "C2'", "C3'", "C4'"]
                                             )
print(amplitude, phaseangle)
>>>  0.352266 85.6
```

### Calculating the pucker coordinates of a *five-membered ring* by the *Cremer-Pople* formalism, from a *xyz* formatted file
```python
from puckepy.formalism import Xyz, CP5
xyz = Xyz("./fivering_3endo.xyz").parse()
amplitude, phaseangle = CP5().from_indices(coordinates=xyz,
                                             indices=[8, 9, 28, 26, 6]
                                             )
print(amplitude, phaseangle)
>>>  0.352266 85.6
```

### Calculating the pucker coordinates of a *five-membered ring* by the *Altona-Sundaralingam* formalism, from a *xyz* formatted file
```python
from puckepy.formalism import Xyz, AS
xyz = Xyz("./fivering_3endo.xyz").parse()
amplitude, phaseangle = AS().from_indices(coordinates=xyz,
                                             indices=[8, 9, 28, 26, 6]
                                             )
print(amplitude, phaseangle)
>>>  0.352266 355.6
```

### Calculating the pucker coordinates of a *six-membered ring* by the *Altona-Sundaralingam* formalism, from a *pdb* formatted file
```python
from puckepy.formalism import Pdb, CP6
pdb = Pdb("./sixring_chair.pdb").parse()
amplitude, phaseangle, theta = CP6().from_atomnames(pdb=pdb,
                                             query_names=["O5'","C1'", "C2'", "C3'", "C4'", "C5'"]
                                             )
print(amplitude, phaseangle, theta)
>>>  0.6587 120.7 1.2
```
### Calculating the pucker coordinates of a *six-membered ring* by the *Strauss-Pickett* formalism, from a *pdb* formatted file
```python
from puckepy.formalism import Pdb, CP6
pdb = Pdb("./sixring_chair.pdb").parse()
alphas, betas = SP().from_atomnames(pdb=pdb,
                                    query_names=["O5'","C1'", "C2'", "C3'", "C4'", "C5'"]
                                    )
print(alphas[0], alphas[1], alphas[2])
>>>  139.8816 146.5372 139.8816
```

### Calculating the pucker coordinates of a *duplex of five-membered ring residues* by the *Cremer-Pople* formalism, from a *pdb* formatted file
```python
from puckepy.formalism import Pdb, CP5
listOfPdb = Pdb("./fivering_3endo.pdb").parse_by_monomers()
for pdbMonomer in listOfPdb :
    amplitude, phaseangle = CP5().from_atomnames(pdb=pdbMonomer,
                                                 query_names=["O4'","C1'", "C2'", "C3'", "C4'"]
    print(amplitude, phaseangle)
>>>  0.352266 85.6
>>>  0.368456 87.1
>>>  ...      ...
```
<br>


## Invert puckering coordinates

### Inverting *Cremer-Pople* coordinates of a *five-membered ring* to its molecular structure
```python
  from puckepy import formalism
  cp5 = formalism.CP5(0.35, 90.)
  inversion = cp5.invert()
  
  formalism.write_to_pdb("inverted_fivering", inversion, "FIV")
```
### Inverting *Cremer-Pople* coordinates of a *six-membered ring* to its molecular structure
```python
  from puckepy import formalism
  cp6 = formalism.CP6(0.67, 120., 1.5)
  inversion = cp6.invert()
  
  formalism.write_to_pdb("inverted_sixring", inversion, "SIX")
```
<br>

## Basic geometry operations

### Calculate *geometry* properties of the *desired molecule*
```python
from puckepy.formalism import Xyz
from puckepy import geometry
coordinates = Xyz("fivering_3endo.xyz").parse()

# Calculate dihedral
dihedral = geometry.dihedral(coordinates[10],
                             coordinates[11],
                             coordinates[12],
                             coordinates[13]
                             )

# Calculate bondangle
angle = geometry.angle(coordinates[10],
                       coordinates[11],
                       coordinates[12]
                       )

# Calculate bondlength
bondlength = geometry.bondlength(coordinates[10],
                                 coordinates[11]
                                 )
```

### Calculate *geometry* properties from *arbitrary data*
```python
from puckepy import geometry

# Calculate dihedral
dihedral = geometry.dihedral([2.23, 0.23, 1.],
                             [3.23, -0.23, 1.],
                             [5.76, 0.01, 2.69],
                             [0.25, 3.44 , 1.93]
                             )

# Calculate bondangle
bondangle = geometry.bondangle([2.23, 0.23, 1.],
                               [3.23, -0.23, 1.],
                               [5.76, 0.01, 2.69]
                               )

# Calculate bondlength
bondlength = geometry.bondlength([2.23, 0.23, 1.],
                                 [3.23, -0.23, 1.]
                                 )
```
