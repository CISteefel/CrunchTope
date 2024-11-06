# Porosity Block


This keyword block is used to set parameters and options related to the treatment of the porosity.
The main contrast is between simulations in which the porosity is fixed and those in which it evolves as the mineral volume fractions evolve.


[Basic keywords](#basic-keywords)
- [fix_porosity](#fix_porosity)
- [minimum_porosity](#minimum_porosity)
- [porosity_update](#porosity_update)
- [read_PorosityFile](#read_porosityfile)

## Basic keywords

### fix_porosity

#### Syntax
```
fix_porosity  [value]
```
[value] is a real number (Default:  If not provided, the code calculates porosity based on the mineral volume fractions).


#### Explanation
The keyword parameter fix_porosity is used to set a global fixed porosity for a simulation.
When provided, it disables the calculation of porosity based on the sum of the volume fractions.
Accordingly, no porosity update based on these quantities is possible.
When this parameter is set in the POROSITY keyword block, set_porosity keywords within individual geochemical conditions can be used once the geochemical condition is distributed in space within the INITIAL_CONDITION keyword block.

#### Example

```
fix_porosity   0.45
```

### minimum_porosity

#### Syntax
```
minimum_porosity  [value]
```
[value] is a real number (Default: $10^{-14}$).

#### Explanation
This keyword is only used to set a minimum to which the porosity can evolve as the mineral volume fractions change.
This option is only used where porosity_update is true and neither fix_porosity nor read_porosity are set.

#### Example

```
minimum_porosity   1e-3
```

### porosity_update

#### Syntax
```
porosity_update [logical]
```
[logical] = true or yes or false or no (Default=false) 

#### Explanation
This keyword followed by a logical (true or yes or false or no) and is used to determine whether the porosity is updated as the mineral volume fractions evolve.
This keyword is only used where neither the fix_porosity nor read_porosity keywords have been set.
When true, it instructs the code to update the porosity as the mineral volume fractions evolve according to
$$\phi = \sum_{n=1}^{N_m} 1 - \phi_m$$
where $N_m$ is the number of minerals in the system and $\phi_m$ is the volume fraction of the individual mineral.

#### Example

```
porosity_update true
```

### read_PorosityFile

#### Syntax
```
read_PorosityFile [filename] [format]
```
[filename] gives the name of the file (up to 132 characters) containing temperatures distributed in space to be read (Default: none). 
[format]= SingleColumn, ContinuousRead, FullFormat, or Unformatted (Default: SingleColumn) 

#### Explanation
This keyword provides the name of file containing porosity values defined at the center of the grid cell for the entire spatial domain.
This option supersedes all other temperature specifications (???).
The format of the file depends on the optional format specified (see Format of Additionnal Input Files).
If no format is specified, the code assumes a single column (SingleColumn) of values in the order 1-NX, 1-NY, 1-NZ:
```
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      READ(52,*) p(jx,jy,jz)
    END DO
  END DO
END DO
```

#### Example

```
read_PorosityFile MyPorosityFile.dat ContinuousRead
```
