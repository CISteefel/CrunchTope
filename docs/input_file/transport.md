## TRANSPORT

The TRANSPORT keyword block reads in the inputs controlling transport of
species. Flow is treated in the FLOW keyword block. The most important
inputs concern molecular diffusion and hydrodynamic dispersion, along
with transport properties of the medium (e.g., the "cementation factor"
or tortuosity) that affect these processes.

#### Fix_diffusion

Keyword followed by a value for the molecular diffusion coefficient in
the space and time units set by the *space_units* and *time_units*
keywords within the TRANSPORT block.

    Syntax:  fix_diffusion   value

<u>Default:</u>  &nbsp; 0.0 (default units of $m^2$/yr)

<u>Explanation:</u> &nbsp;  If this parameter is set, then additional
instructions to calculate diffusion (*calculate_diffusion*) based on a
25°C value and an activation energy are ignored. If no species-specific diffusion coefficients 
are specified (see *D_25*
keyword) With setting of this parameter, all aqueous species are given
the same diffusion coefficient and only Fick's law is implemented in the
calculation (the electrochemical migration term is zero when all species
have the same diffusion coefficient in any case). If any
species-specific diffusion coefficients are set, then the full form of
the Nernst-Planck equation is solved and the value
given by *fix_diffusion* is used for any species not given a specific
diffusion coefficient value with the *D_25* keyword.

#### Calculate_diffusion

Keyword followed by a value for the molecular diffusion coefficient in
the space and time units set by the *space_units* and *time_units*
keywords within the TRANSPORT block.

    Syntax:  calculate_diffusion  value

<u>Default:</u> &nbsp; 0.0 (default units of $m^2$/yr)

<u>Explanation:</u> &nbsp;  This parameter is used to calculate a
diffusion coefficient at 25°C. Diffusion coefficients at other
temperatures are calculated with the value set in the
*diffusion_activation* keyword The fix_diffusion keyword over-rides this
keyword (for clarity, both should not normally be set).

#### Cementation_exponent

Keyword followed by a value for the cementation exponent in Archie's
Law.

    Syntax:  cementation_exponent  value

<u>Default:</u> &nbsp; 1.0

<u>Explanation:</u> &nbsp;  This parameter provides an exponent to be
used in Archie's Law to specify a porosity-dependent tortuosity
according to:

$$D_{eff} = \phi^{m}D_{w}$$

where *m* is the cementation exponent, $D_{w}$ is the diffusion
coefficient in water, and *$D_{effective}$is the effective diffusion
coefficient in porous medium. As defined here, *m* = 1 would be the
formulation in porous media in the case where there is no tortuosity
effect and the diffusion coefficient is corrected only for the
cross-sectional area actually made up of pores. According to this
formulation, the formation factor would be defined as:

$$F = \phi^{- m}$$

#### Formation_factor

Keyword followed by a value for the formation factor.

    Syntax:  formation_factor  value

<u>Default:</u> &nbsp; 1.0

<u>Explanation:</u> &nbsp;  This parameter provides an a formation
factor, which is defined as:

$$D_{eff} = \frac{D_{w}}{F}$$

where *F* is the formation factor, $D_{w}$ is the diffusion coefficient
in water, and $D_{effective}$ is the effective diffusion coefficient in porous
medium. According to the definition of the formation factor in use in
the marine geochemistry literature, for example, the cementation
exponent should be set to 0 (i.e, the porosity dependence of the
effective diffusion coefficient is incorporated into the formation
factor.

#### Diffusion_activation

Keyword followed by a value for the molecular diffusion activation
energy in kcal/mole.

    Syntax:  diffusion_threshold  value

<u>Default:</u> &nbsp; $D_{w}$ = 5.0 (default units of kcal/mole)

<u>Explanation:</u> &nbsp;  This parameter sets the activation energy to
be used in calculating diffusion coefficients at temperatures other than
25°C.

#### D_25

Keyword specifying a diffusion coefficient at 25°C for a specific
aqueous species.

    Syntax:  D_25 species_name  value

<u>default:</u> &nbsp;  value set in fix_diffusion keyword

<u>Explanation:</u> &nbsp;  This keyword provides a diffusion coefficient
in the units based on the space_units and time_units keywords set in the
TRANSPORT block. When even one species_specific diffusion coefficient is
given, the code (GIMRT only) solves the Nernst-Planck equation, which
includes an explicit calculation of electrochemical migration. In this
case, any species not listed with this keyword are given the diffusion
coefficient specified in the *fix_diffusion* keyword.

#### Gas_diffusion

Keyword followed by a value for the gas diffusion coefficient.

    Syntax:  gas_diffusion  value

<u>Default:</u> &nbsp;  0.0 (default units of $m^2$/yr)

<u>Explanation:</u> &nbsp;  This keyword sets a value for the gas
diffusion coefficient in the units given in the *time_units* and
*space_units* keywords set within the TRANSPORT block.

#### Dispersivity

Keyword followed by a value for the gas diffusion coefficient.

    Syntax:  dispersivity  <longitudinal value>  <transverse value>

<u>Default:</u> &nbsp;  0.0 (default units of m)

<u>Explanation:</u> &nbsp;  This keyword sets a value for the
longitudinal (first entry) and transverse (second value) dispersivities.
Space units may be changed with the *space_units* keyword set within the
TRANSPORT block. If the problem is one-dimensional, the transverse
dispersivity may be omitted.

#### Constant_tortuosity

Keyword followed by a value for the tortuosity, or by a string giving
the name of the tortuosity option.

    Syntax:  constant_tortuosity <value> or constant_tortuosity <option>

<u>Default:</u> &nbsp; 1.0

<u>Explanation:</u> &nbsp;  This keyword sets a constant value for the
tortuosity, which is then applied throughout the spatial domain.
Definitions of the tortuosity vary---in most cases, it does not include
the porosity that multiplies the diffusion coefficient and this is the
convention adopted here. The tortuosity is also considered here to be a
number less than 1 (some formulations divide by the tortuosity, so that
it is always greater than 1). The tortuosity (τ) as used in CrunchFlow,
therefore, can be viewed as a parameter which allows an effective
diffusion coefficient in porous media (*D~PM~*) to be calculated from
the diffusion coefficient in pure water (*D~w~*) and the porosity
according to:

$$D_{PM} = \phi\tau D_{w}$$

If an option rather than a numerical value is given, the code checks the
option specified against the list currently implemented. At the time of
this manual, the only option that is accepted is *CostaRica*, which is a
porosity-tortuosity option to published soon. With an option set, the
tortuosity is calculated dynamically at run time (so the keyword
*constant_tortuosity* is a bit of a misnomer here), but it avoids a read
of tortuosity from an input file (*read_tortuosity*) or by field
(*tortuosity*).

#### Read_TortuosityFile

Keyword followed by a file name containing tortuosity values over the
entire spatial domain.

    Syntax:  read_TortuosityFile  filename  format

<u>Backwards compatibility:</u> &nbsp; read_tortuosity

<u>Default:</u> &nbsp; None

<u>Explanation:</u> &nbsp;  This keyword provides a filename for a file
containing tortuosity values over the entire spatial domain. The
tortuosity is considered here to be a number less than 1 (some
formulations divide by the tortuosity, so that it is always greater than
1). The tortuosity (τ) as used in CrunchFlow, therefore, can be viewed
as a parameter which allows an effective diffusion coefficient in porous
media ($D_{PM}$) to be calculated from the diffusion coefficient in pure
water ($D_{w}$) and the porosity according to:

$$D_{PM} = \phi\tau D_{w}$$

The format of the file depends on the format specified after the file
name, but the values of the tortuosity are assumed to be listed with NX
varying first, then NY, and then NZ. For those familiar with FORTRAN and
using the *SingleColumn* format (the default), the actual source code
is:

    Do jz = 1,nz
      Do jy = 1,ny
        Do jx = 1,nx
          READ(52,*) tortuosity(jx,jy,jz)
        End Do
      End Do
    End Do

#### Tortuosity

Keyword used to specify tortuosity by zone in the input file.

    Syntax:  tortuosity  value  zone  jxbegin-jxend  jybegin-jyend  jzbegin-jzend

or

    Syntax:  tortuosity value

<u>Default:</u> &nbsp; 1.0

<u>Explanation:</u> &nbsp;  The *tortuosity* keyword uses an input format
similar to other parameters specified by zone (e.g., pressure,
permeability). The keyword is followed by the value of the tortuosity,
which is then followed by either the label *zone* or *default*. If
*zone* is specified, then the code expects the X, Y, and Z coordinates
in the form of a beginning jx, and ending jx, a beginning jy and an
ending jy, and a beginning jz and an ending jz, in each case separated
by a hyphen. If *default* is specified instead of *zone*, the code will
initialize the entire spatial domain (i.e., all grid cells) to the
value. More than one specification of *tortuosity* can be (and normally
is) provided, for example

    tortuosity   0.1  default
    tortuosity  0.01  zone  1-20  1-1  1-1
    tortuosity  0.02  zone  1-1  1-20  1-1

Zones specified later in the sequence overwrite zones specified
earlier---in the example above, the default specification sets the
tortuosity value to 0.1 initially in all of the grid cells, while the
next tortuosity keyword changes the value to 0.01 in grid cells jx=1 to
jx=20. *The convention adopted here is that all three dimensions must be
specified, even in a one-dimensional problem.

The tortuosity is considered here to be a number less than 1 (some
formulations divide by the tortuosity, so that it is always greater than
1). The tortuosity (τ) as used in CrunchFlow, therefore, can be viewed
as a parameter which allows an effective diffusion coefficient in porous
media ($D_{PM}$) to be calculated from the diffusion coefficient in pure
water ($D_{w}$) and the porosity according to:

$$D_{PM} = \phi(jx,jy,jz)\tau(jx,jy,jz)D_{w}$$

#### Anisotropy_ratioY

Keyword followed by a value defining the anisotropy ratio in the Y
direction relative to the value in the X direction.

    Syntax:  Anisotropy_ratioY  value

<u>Default:</u> &nbsp; 1.0

<u>Explanation:</u> &nbsp;  This keyword defines the ratio of the
diffusion coefficient in the Y direction relative to the diffusion
coefficient in the X direction. The anisotropy term, therefore, is
dimensionless. If the anisotropy factor in the Y direction is defined as
ε~Y~, then the effective diffusion coefficient in the Y direction (JY)
is given by:

$$D_{PM} = \varepsilon_{Y}\phi\tau D_{w}$$

#### Anisotropy_ratioZ

Keyword followed by a value defining the anisotropy ratio in the Z
direction relative to the value in the X direction.

    Syntax:  anisotropy_ratioZ  value

<u>Default:</u> &nbsp; 1.0

<u>Explanation:</u> &nbsp;  This keyword defines the ratio of the
diffusion coefficient in the Z direction relative to the diffusion
coefficient in the X direction. The anisotropy term, therefore, is
dimensionless. If the anisotropy factor in the Z direction is defined as
ε~Z~, then the effective diffusion coefficient in the Z direction (JZ)
is given by:

$$D_{PM} = \varepsilon_{Z}\phi\tau D_{w}$$

#### Threshold_porosity

Keyword followed by a value defining a porosity threshold separating two
different tortuosity values.

    Syntax:  threshold_porosity  value

<u>Default:</u> &nbsp; 0.0

<u>Explanation:</u> &nbsp;  This keyword sets a threshold porosity,
either side of which the tortuosity takes on differing values given by
the keywords *tortuosity_below* and *tortuosity_above*.

#### Tortuosity_below

Keyword followed by a value giving the tortuosity below the porosity
threshold set in the keyword *threshold_porosity*..

    Syntax:  tortuosity_below  value

<u>Default:</u> &nbsp;  *.0 (only used where threshold_porosity is set)

<u>Explanation:</u> &nbsp;  This keyword sets the value of the tortuosity
to be used below the threshold porosity. This keyword is only used if
*threshold_porosity* is set.

#### Tortuosity_above

Keyword followed by a value giving the tortuosity above the porosity
threshold set in the keyword *threshold_porosity*..

    Syntax:  tortuosity_above  value

<u>Default:</u> &nbsp;  0.0 (only used where threshold_porosity is set)

<u>Explanation:</u> &nbsp;  This keyword sets the value of the tortuosity
to be used above the threshold porosity. This keyword is only used if
*threshold_porosity* is set.