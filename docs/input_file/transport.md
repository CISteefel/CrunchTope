# Transport Block

The TRANSPORT keyword block reads in the inputs controlling transport of species.
Flow is treated in the FLOW keyword block.
The most important inputs concern molecular diffusion and hydrodynamic dispersion, along with transport properties of the medium (e.g., the “cementation factor” or tortuosity) that affect these processes.


[Basic keywords](#basic-keywords)
- [Fix_diffusion](#fix_diffusion)
- [Calculate_diffusion](#calculate_diffusion)
- [Cementation_exponent](#cementation_exponent)
- [Formation_factor](#formation_factor)
- [Diffusion_activation](#diffusion_activation)
- [D_25](#D_25)
- [Gas_diffusion](#gas_diffusion)
- [Dispersivity](#dispersivity)
- [Constant_tortuosity](#constant_tortuosity)
- [Read_TortuosityFile](#read_tortuosityfile)
- [Tortuosity](#tortuosity)
- [Anisotropy_ratioY](#anisotropy_ratioy)
- [Anisotropy_ratioZ](#anisotropy_ratioz)
- [Threshold_porosity](#threshold_porosity)
- [Tortuosity_below](#tortuosity_below)
- [Tortuosity_above](#tortuosity_above) 

## Basic keywords

### Fix_diffusion

#### Syntax
```
fix_diffusion  [value]
```
Keyword followed by a value for the molecular diffusion coefficient in the space and time units set by the space_units and time_units keywords within the TRANSPORT block.
[value] (Default=0.0 (default units of m2/yr))

#### Explanation
If this parameter is set, then additional instructions to calculate diffusion (calculate_diffusion) based on a 25°C value and an activation energy are ignored.
If no species-specific diffusion coefficients are specified (see D_25 keyword) With setting of this parameter, all aqueous species are given the same diffusion coefficient and only Fick’s law is implemented in the calculation (the electrochemical migration term is zero when all species have the same diffusion coefficient in any case).
If any species-specific diffusion coefficients are set, then the full form of the Nernst-Planck equation is solved (see section XXX) and the value given by fix_diffusion is used for any species not given a specific diffusion coefficient value with the D_25 keyword.

#### Example

```
fix_diffusion 1.0e-09
```

### Calculate_diffusion

#### Syntax
```
calculate_diffusion  [value]
```
Keyword followed by a value for the molecular diffusion coefficient in the space and time units set by the space_units and time_units keywords within the TRANSPORT block.
[value] (Default=0.0 (default units of m2/yr))

#### Explanation
This parameter is used to calculate a diffusion coefficient at 25°C.
Diffusion coefficients at other temperatures are calculated with the value set in the diffusion_activation keyword.
The fix_diffusion keyword over-rides this keyword (for clarity, both should not normally be set).

#### Example

```
calculate_diffusion
```

### Cementation_exponent

#### Syntax
```
cementation_exponent  [value]
```
Keyword followed by a value for the cementation exponent in Archie’s Law.
[value] (Default=1.0)

#### Explanation
This parameter provides an exponent to be used in Archie’s Law to specify a porosity-dependent tortuosity according to:
$$D_{eff} = \phi^m D_w$$
where $m$ is the cementation exponent, $D_w$ is the diffusion coefficient in water, and $D_{eff}$ is the effective diffusion coefficient in porous medium.
As defined here, $m = 1$ would be the formulation in porous media in the case where there is no tortuosity effect and the diffusion coefficient is corrected only for the cross-sectional area actually made up of pores.
According to this formulation, the formation factor would be defined as:
$$F = \phi^{-m}$$


#### Example

```
cementation_exponent 1.0
```

### Formation_factor

#### Syntax
```
formation_factor [value]
```
Keyword followed by a value for the formation factor.
[value] (Default=1.0)

#### Explanation
This parameter provides an a formation factor, which is defined as:
$$D_{eff} = \frac{D_w}{F}$$
where $F$ is the formation factor, $D_w$ is the diffusion coefficient in water, and $D_{eff}$ is the effective diffusion coefficient in porous medium.
According to the definition of the formation factor in use in the marine geochemistry literature, for example, the cementation exponent should be set to 0 (i.e, the porosity dependence of the effective diffusion coefficient is incorporated into the formation factor.  

#### Example

```
formation_factor 1.0
```

### Diffusion_activation

#### Syntax
```
diffusion_activation [value]
```
Keyword followed by a value for the molecular diffusion activation energy in kcal/mole.
[value] (Default=5.0 (default units of kcal/mole))

#### Explanation
This parameter sets the activation energy to be used in calculating diffusion coefficients at temperatures other than 25°C.

#### Example

```
diffusion_activation 5.0
```

### D_25

#### Syntax
```
D_25 [species_name] [value]
```
Keyword specifying a diffusion coefficient at 25°C for a specific aqueous species.
[value] (Default=value set in fix_diffusion keyword)

#### Explanation
This keyword provides a diffusion coefficient in the units based on the space_units and time_units keywords set in the TRANSPORT block.
When even one species_specific diffusion coefficient is given, the code (GIMRT only) solves the Nernst-Planck equation, which includes an explicit calculation of electrochemical migration.
In this case, any species not listed with this keyword are given the diffusion coefficient specified in the fix_diffusion keyword.

#### Example

```
D_25
```

### Gas_diffusion

#### Syntax
```
gas_diffusion [value]
```
Keyword followed by a value for the gas diffusion coefficient.
[value] (Default=0.0 (default units of m2/yr))

#### Explanation
This keyword sets a value for the gas diffusion coefficient in the units given in the time_units and space_units keywords set within the TRANSPORT block.

#### Example

```
Gas_diffusion 0.0
```

### Dispersivity

#### Syntax
```
dispersivity [longitudinal_value] [transverse_value]
```
[value] (0.0 (default units of m))

#### Explanation
This keyword sets a value for the longitudinal (first entry) and transverse (second value) dispersivities.
Space units may be changed with the space_units keyword set within the TRANSPORT block.
If the problem is one-dimensional, the transverse dispersivity may be omitted.

#### Example

```
Dispersivity 0.01 0.001
```

### Constant_tortuosity

#### Syntax
```
constant_tortuosity [value] or [constant_tourtuosity option]
```
Keyword followed by a value for the tortuosity, or by a string giving the name of the tortuosity option.
[value] (1.0)

#### Explanation
This keyword sets a constant value for the tortuosity, which is then applied throughout the spatial domain.
Definitions of the tortuosity vary—in most cases, it does not include the porosity that multiplies the diffusion coefficient and this is the convention adopted here.
The tortuosity is also considered here to be a number less than 1 (some formulations divide by the tortuosity, so that it is always greater than 1).
The tortuosity ($\tau$) as used in CrunchFlow, therefore, can be viewed as a parameter which allows an effective diffusion coefficient in porous media (DPM) to be calculated from the diffusion coefficient in pure water ($D_w$) and the porosity according to:
$$D_{PM} = \phi \tau D_w$$.
If an option rather than a numerical value is given, the code checks the option specified against the list currently implemented.
At the time of this manual, the only option that is accepted is CostaRica, which is a porosity-tortuosity option to published soon.
With an option set, the tortuosity is calculated dynamically at run time (so the keyword constant_tortuosity is a bit of a misnomer here), but it avoids a read of tortuosity from an input file (read_tortuosity) or by field (tortuosity).


#### Example

```
constant_tortuosity 1.0
```

### Read_TortuosityFile

#### Syntax
```
read_tortuosityFile [filename] [format]
```
Keyword followed by a file name containing tortuosity values over the entire spatial domain.
Default=None

#### Explanation
This keyword provides a filename for a file containing tortuosity values over the entire spatial domain.
The tortuosity is considered here to be a number less than 1 (some formulations divide by the tortuosity, so that it is always greater than 1).
The tortuosity ($\tau$) as used in CrunchFlow, therefore, can be viewed as a parameter which allows an effective diffusion coefficient in porous media (DPM) to be calculated from the diffusion coefficient in pure water (Dw) and the porosity according to:
$$D_{PM} = \phi \tau D_w$$.
The format of the file depends on the format specified after the file name, but the values of the tortuosity are assumed to be listed with NX varying first, then NY, and then NZ.
For those familiar with FORTRAN and using the SingleColumn format (the default), the actual source code is:
```
Do jz = 1,nz
  Do jy = 1,ny
    Do jx = 1,nx
      READ(52,*) tortuosity(jx,jy,jz)
    End Do
  End Do
End Do
```

#### Example

```
read_tortuosityFile
```

### Tortuosity

#### Syntax
```
tortuosity [value] zone  [jxbegin-jxend]   [jybegin-jyend]   [jzbegin-jzend]
```
or
```
tortuosity [value] default
```
Keyword used to specify tortuosity by zone in the input file.
[value] (Default=1.0)

#### Explanation
The tortuosity keyword uses an input format similar to other parameters specified by zone (e.g., pressure, permeability).
The keyword is followed by the value of the tortuosity, which is then followed by either the label zone or default.
If zone is specified, then the code expects the X, Y, and Z coordinates in the form of a beginning JX, and ending JX, a beginning JY and an ending JY, and a beginning JZ and an ending JZ, in each case separated by a hyphen.
If default is specified instead of zone, the code will initialize the entire spatial domain (i.e., all grid cells) to the value.
More than one specification of tortuosity can be (and normally is) provided, for example
```
tortuosity  0.1   default
tortuosity  0.01  zone  1-20 1-1 1-1
tortuosity  0.02  zone  1-1 1-20 1-1
```
Zones specified later in the sequence overwrite zones specified earlier—in the example above, the default specification sets the tortuosity value to 0.1 initially in all of the grid cells, while the next tortuosity keyword changes the value to 0.01 in grid cells JX=1 to JX=20.
The convention adopted here is that all three dimensions must be specified, even in a one-dimensional problem.
The tortuosity is considered here to be a number less than 1 (some formulations divide by the tortuosity, so that it is always greater than 1).
The tortuosity ($\tau$) as used in CrunchFlow, therefore, can be viewed as a parameter which allows an effective diffusion coefficient in porous media (DPM) to be calculated from the diffusion coefficient in pure water (Dw) and the porosity according to:
$$D_{PM} = \phi(jx, jy, jaz) \tau(jx, jy, jz) D_w$$.


#### Example

```
tortuosity
```

### Anisotropy_ratioY

#### Syntax
```
Anisotropy_ratioY [value]
```
[value] (Default=1.0)

#### Explanation
This keyword defines the ratio of the diffusion coefficient in the Y direction relative to the diffusion coefficient in the X direction.
The anisotropy term, therefore, is dimensionless.
If the anisotropy factor in the Y direction is defined as $\eta_Y$, then the effective diffusion coefficient in the Y direction (JY) is given by:
$$D_{PM} = \eta_Y \phi \tau D_w$$.

#### Example

```
Anisotropy_ratioY 1.0
```

### Anisotropy_ratioZ

#### Syntax
```
Anisotropy_ratioZ [value]
```
[value] (Default=1.0)

#### Explanation
This keyword defines the ratio of the diffusion coefficient in the Z direction relative to the diffusion coefficient in the X direction.
The anisotropy term, therefore, is dimensionless.
If the anisotropy factor in the Z direction is defined as $\eta_Z$, then the effective diffusion coefficient in the Z direction (JZ) is given by:
$$D_{PM} = \eta_Z \phi \tau D_w$$.

#### Example

```
Anisotropy_ratioZ 1.0
```

### Threshold_porosity

#### Syntax
```
threshold_porosity [value]
```
Keyword followed by a value defining a porosity threshold separating two different tortuosity values.
[value] (Default=0.0)

#### Explanation
This keyword sets a threshold porosity, either side of which the tortuosity takes on differing values given by the keywords tortuosity_below and tortuosity_above.

#### Example

```
threshold_porosity
```

### Tortuosity_below

#### Syntax
```
tortuosity_below [value]
```
Keyword followed by a value giving the tortuosity below the porosity threshold set in the keyword threshold_porosity.
[value] (Default=0.0 (only used where threshold_porosity is set))

#### Explanation
This keyword sets the value of the tortuosity to be used below the threshold porosity.
This keyword is only used if threshold_porosity is set.

#### Example

```
tortuosity_below
```

### Tortuosity_above

#### Syntax
```
tortuosity_above [value]
```
Keyword followed by a value giving the tortuosity above the porosity threshold set in the keyword threshold_porosity.
[value] (Default=0.0 (only used where threshold_porosity is set))

#### Explanation
This keyword sets the value of the tortuosity to be used above the threshold porosity.
This keyword is only used if threshold_porosity is set.

#### Example

```
tortuosity_avobe
```
