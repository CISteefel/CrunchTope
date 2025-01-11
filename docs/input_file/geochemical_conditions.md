### GEOCHEMICAL CONDITIONS

Geochemical conditions are keyword blocks which may occur more than once
and are used to set up boundary and initial conditions and source terms.
In each case, they provide information on the actual concentrations of
primary species and minerals needed to carry out a distribution of
species calculation.

Only after a particular geochemical condition is read and "processed" is
it allocated as a particular boundary or initial condition or as a
source term.

Geochemical conditions are input via a CONDITION keyword block. The
CONDITION string is followed in each case by a label which identifies
the condition elsewhere in the input file. Most of the information in a
condition block pertains to aqueous and solid phase
concentrations---however, certain parameters may be input which apply to
the condition as a whole. These are discussed in the following section.

#### Units

Concentration units may be set using the units keyword within a
particular condition. The change in concentration units only applies to
the condition within which it is set. The format for input of different
units is:

    units <concentration units>

The default concentration units for aqueous species is molality (moles
solute per kg water). The additional options include:

Table 3: Options for specifying concentrations in geochemical conditions
  | Concentration Units | Keyword Designator   |
  | ------------------- | ------------------   |
  | moles/kg water      |  mol/kg              |
  | millimoles/kg water |  mmol/kg             |
  | micromoles/kg water | umol/kg, micromol/kg |
  | parts per million   | ppm                  |

In each case, the use of the "ppm" units refers to the parts per million
of the solute using the molecular weight.of the primary species to which
it refers. Therefore, if the primary species is NO$_3^-$, then the
concentration in ppm is calculated using the molecular weight of NO$_3^-$
(not, for example, using N).

#### Equilibrate_surface

When specified on a separate line, this keyword instructs the code to
partition all of the total concentrations between the aqueous and solid
phases (i.e., surface hydroxyl sites and exchangers). It does not
partition any mass through precipitation-dissolution reaction. This
instruction will only apply to species for which a "total concentration"
constraint is used (described below). This instruction can also be
specified on a species by species basis as discussed below. If the
equilibrate_surface instruction is absent, then the code assumes that
the "total concentration" refers only to the aqueous phase and the
exchanger and surface hydroxyl concentrations are calculated
accordingly.

#### Temperature

Keyword used to set a temperature within an individual geochemical
condition.

    Syntax:  temperature  <value>

<u> Default </u>:  Provided by value in TEMPERATURE keyword block

Explanation:  This keyword in the CONDITION overrides the value set in the
**TEMPERATURE** keyword block. If absent in the CONDITION, the temperature specified in the **TEMPERATURE** keyword block is used for the
geochemical condition.

#### Set_porosity

Keyword used to set a porosity within an individual geochemical
condition.

    Syntax:  set_porosity <value>

<u> Default </u>:  Provided by value in POROSITY keyword block

Explanation:  This keyword overrides the value set in the
**POROSITY** keyword block. If absent in the CONDITION, that value is used for the
geochemical condition.

<u>NOTE:</u> This option only works if **fix_porosity** is set
in the **POROSITY** keyword block. Otherwise, porosity is calculated
from the 1 minus the sum of the volume fractions of the minerals.*

#### Set_saturation

Keyword used to set a liquid saturation within an individual geochemical
condition.

    Syntax:  set_saturation <value>

<u> Default </u>:  Provided by value in RUNTIME keyword block,
otherwise 1.0

Explanation:  This keyword overrides the value set in the
**RUNTIME** keyword block. If absent, that value is used for the
geochemical condition.

#### Types of Constraints: Aqueous Species

Various types of constraints may be specified for aqueous species in the
distribution of species calculations carried out on each of the
geochemical conditions. These constraints are summarized in Table 4.

Table 4: Constraints available for specifying aqueous species concentrations
  | Constraint Type | Use | Example |
  | --------------- | --- | ------- |
  | Total | Mole balance on total aqueous concentration or total aqueous+sorbed concentrations | Na+ 0.001 |
  | Species Concentration | Specify an individual primary species concentration | Na+ 0.001 species |
  | Species Activity      | Specify an individual primary species activity | Na+ 0.001 activity
  | pH      | Specify activity of hydrogen ion  | pH 7.8 |
  | Gas     | Equilibrate primary species with a gas | O2(aq) O2(g) 0.20 |
  | Mineral | Equilibrate primary species with a mineral | O2(aq) Pyrite |
  | Charge  | Charge balance the solution using the primary species | Na+ charge |
  ------------------ ----------------------------- ----------------------
The most commonly used constraint is the "total concentration"
constraint. Analyses of natural waters generally report total
concentrations rather than individual species concentrations. Or someone
may know how much total mass of a chemical element they added to an
experimental system. While specification of individual species
concentrations is relatively uncommon, specification of species activity
is common, especially for hydrogen ion (pH) and for molecular oxygen
(O$_{2(aq)}$).

#### Total Concentration Constraint

The total concentration constraint signals to CrunchFlow to carry out a
total mole balance calculation. If the balance is over the aqueous phase
only (the default), the equation takes the form:

$$T_{j} = C_{j} + \sum_{i = 1}^{Ns}{\nu_{ij}C_{i}}$$

where $T_j$ is the total concentration of the primary species, $C_j$ is
the individual primary species concentration, $\nu_{ij}$ is the
stoichiometric coefficient, and $C_i$ is the concentration of the *Ns*
secondary species. If the keyword *equilibrate_surface* is specified
either as a global parameter within a condition or for an individual
primary species, then the mole balance equation also includes the mass
of the primary species incorporated in ion exchange sites and as surface
complexes:

$$T_{j} = C_{j} + \sum_{i = 1}^{Ns}{\nu_{ij}C_{i}} + \sum_{k = 1}^{Nexs}{\nu_{kj}C_{k}} + \sum_{l = 1}^{Nsurf}{\nu_{lj}C_{l}}$$

where $C_k$ refers to the *Nex* exchange species and $C_l$ refers to the
*Nsurf* surface complexes. To specify a total concentration within a
condition, one uses the form:

    <PrimarySpeciesName> <Concentration>

or equivalently,

    <PrimarySpeciesName> <Concentration> total

where the "total" is assumed to be the default. To distribute the
specified mass over the exchangers and surface hydroxyl sites as well,
the user follows the concentration value with the keyword
*equilibrate_surface* as in the first three entries of the following
example:

>     Am+++ 1.00e-10 equilibrate_surface
>     Np++++ 1.00e-10 equilibrate_surface
>     Pu++++241 1.00e-10 equilibrate_surface
>     K+ 8.70e-5
>     Mg++ 2.06e-5

In this example, the Am, Np, and Pu will be partitioned between the
aqueous phase and the exchange and surface hydroxyl sites (if present),
while the total concentration in the case of K and Mg will be
distributed completely within the aqueous phase.

#### Species Concentration Constraint

The species concentration constraint is straightforward and need not be
discussed here other than to show how it is specified in the input file.
Since the total concentration constraint is considered the default, it
is necessary to specify the keyword species after the concentration
value:

     <PrimarySpeciesName> <Concentration> species

#### Species Activity Constraint

Specifying a species activity is also straightforward and takes the
form:

    <PrimarySpeciesName> <Concentration> activity

A specification of the pH is a special case which does not require the
trailing *activity* keyword

>     pH 7.8

#### Gas Constraint

The gas constraint is used to calculate a primary species activity so as
to obtain equilibrium with respect to a gas at a specified partial
pressure/fugacity. At this point, no fugacity corrections are available
in the code, so effectively partial pressure = fugacity. Partial
pressure is input in bars. The format for using the gas constraint is to
follow the primary species name by the name of the gas with which the
species is to be equilibrated and then the partial pressure in bars:

    <primary species name> <gas name> <partial pressure (bars)>

To adjust the O2(aq) activity to be at equilibrium with respect to
atmospheric oxygen, one would use the following form:

    O2(aq) O2(g) 0.20

#### Mineral Constraint

The mineral constraint operates in a similar fashion to the gas
constraint except that no trailing value is needed. This is because
minerals (or solid phases) are assumed at this stage to be pure, that
is, their activities = 1. Unlike some other codes like PHREEQC,
CrunchFlow does not carry out a formal reaction path calculation to
attain equilibrium with respect to the mineral. The mineral constraint
is used together with the other constraints in the geochemical system to
provide the same number of constraint equations as there are unknown and
this system of equations is solved without any path dependence. The
format for a mineral constraint is simply:

    <primary species name> <mineral name>

#### Charge Constraint

The charge constraint may be somewhat more difficult to implement
because the user may choose to balance on an anion when the remainder of
the solution is already negatively charged or vice versa. In this case
(without the proper safeguards), the concentration of the
charge-balancing species would go to 0. If this situation arises, an
error message will result and the run will terminate. The format for
specifying a charge balance is:

    <primary species name> charge

#### Solid Phase Inputs

In most cases, information on mineral concentrations (or volume
fractions) and mineral surface areas are needed to completely describe a
geochemical condition. Those cases which don't require such information
are typically boundary conditions or source terms where only the aqueous
phase is entering the domain. In this case it is possible to leave the
mineral input out of the condition.

The code does not check to see whether mineral input is or is not
needed---if it doesn't find any information on a mineral or solid phase
which has been loaded in the **MINERALS** keyword block, it will assume
its volume fraction = 0 and use a default bulk surface area of 100 m^2^
mineral/m^3^ porous medium. In this way, a mineral which is neglected in
a particular condition is treated as a "potential secondary phase", that
is, it is not initially present but it could form if supersaturated.
Since boundary conditions and source terms represent "ghost cells or
zones" outside of the domain which are therefore not computed as a
function of time, supersaturation of a mineral phase will have no
effect.

#### Mineral Volume Fractions

The first values to be input following a mineral/solid phase name is the
volume fraction. The ability to add mineral concentrations in other
units (mass fraction, etc.) will be added in future releases. The format
is:

    <name of solid phase> <volume fraction>

Surface area information is appended to the same line and is described
below.

#### Mineral Surface Area Options

Mineral surface area is important in two ways: 1) it gives the "reactive
surface area" to be used in mineral dissolution and precipitation
reactions, and 2) it provides the surface area upon which adsorption
(surface complexation) may occur (see section on Surface Complexation
Input). There are currently two options for specifying mineral surface
area: 1) *bulk surface area* in units of m^2^ solid phase/m^3^ porous
medium, or 2) *specific surface area* in units of m^2^/g. The format is:

    <name of solid phase> <volume fraction> <surface area option> <value>

where the surface area options are *bulk_surface_area* or *bsa* and
*specific_surface_area* or *ssa* respectively for the bulk and specific
surface area options. If the surface area option is left off, the code
assumes by default that the trailing numerical value refers to the bulk
surface area. Specification of a bulk surface area will lead to a
calculation of specific surface area according to:

$$A_{specific} = \frac{A_{bulk}V_{m}}{\phi_{m}MW_{m}}$$

where *V~m~* is the molar volume of the solid phase, $\phi_{m}$is the
initial volume fraction, and *MW~m~* is the molecular weight of the
phase. The specific surface area is then used together with other
surface complexation parameters to calculate concentrations of surface
hydroxyl sites per volume porous medium.

The two surface area options lead to different methods for computing
mineral surface area as a function of time.

[Bulk surface area:]{.underline} In the case where the bulk surface area
is specified, the reactive surface area of a solid phase as a function
of time is calculated according to:

$$A_{bulk} = A_{bulk}^{initial}\left( \frac{\phi_{m}}{\phi_{m}^{initial}} \right)^{2/3}\left( \frac{\phi}{\phi_{}^{initial}} \right)^{2/3}$
(dissolution)$$

$$A_{bulk} = A_{bulk}^{initial}\left( \frac{\phi}{\phi_{}^{initial}} \right)^{2/3}$
(precipitation)$$

where $\phi$ refers to the porosity and $\phi_{m}$ refers to the
individual mineral volume fraction. The inclusion of a 2/3 dependence on
porosity is chiefly to ensure that as the porosity goes to 0, so too
does the mineral surface area available for reaction. This formulation
is used primarily for primary minerals (that is, minerals with initial
volume fractions \> 0). For secondary minerals which precipitate, the
value of the *initial* bulk surface area specified is used as long as
precipitation occurs---if this phase later dissolves, the above
formulation is used, but with an arbritrary "initial volume fraction" of
0.01. Clearly the expressions above are only a very approximate means of
calculating the reduction in surface area of a secondary phase if later
dissolution occurs. More accurate for secondary minerals (and therefore
recommended) is the use of specific surface area.

[Specific surface area:]{.underline} When specific surface area is
selected in the input file, the reactive surface area used for mineral
precipitation and dissolution is calculated from:

$$A_{bulk} = \frac{\phi_{m}A_{specific}MW_{m}}{V_{m}}$$

This is the preferred approach for precipitation of secondary phases in
particular, since the bulk surface area (BSA) approach currently uses
the initial surface area only (i.e., the reactive surface area of a
precipitating phase does not change as the secondary mineral
concentration evolves). In contrast, with the specific surface area
option, evolution of the mineral volume fraction causes the bulk surface
area to evolve according to the equation above.

If the initial volume fraction of the solid phase is \> 0, then the
$A_{bulk}^{initial}$ is calculated directly from the specific surface
area. If the initial volume fraction = 0 (that is, we are dealing with a
secondary mineral initially not present at all), the specific surface
area cannot be evaluated without some other procedure. To handle this
case of a secondary phase not initially present, one can specify a
"threshold mineral volume fraction", which is used to calculate the bulk
surface area until the computed time-evolving volume fraction exceeds
the threshold value. A threshold value for the case of an initial volume
fraction = 0 is set by

    <name of solid phase> 0.0 specific_surface_area <value> <threshold volume fraction>

So, for example, we could specify that the bulk surface area of calcite
(which is initially not present) be calculated from a threshold volume
fraction of 0.0001 and a specific surface area of 2 m^2^/g with the
following form:

    Calcite 0.0 specific_surface_area 2.0 0.0001

With this formulation, the bulk surface area of calcite would be
calculated directly from the calcite volume fraction and specific
surface area once it exceeded the threshold value of 0.0001 set her.
This is a simple if not particularly mechanistic way to capture a
nucleation event.

#### Surface Hydroxyl Concentrations

If surface complexation is included in a particular problem, then
surface hydroxyl concentrations must be specified in each of the
geochemical conditions---there is currently no default value available.
The information which needs to be provided is the surface hydroxyl site
density, $\rho_{sites}$, in units of moles sites per m^2^ mineral. This
information is combined with the specific surface area, $A_{specific}$,
to calculate a concentration of surface hydroxyl sites per m^3^ porous
medium:

$$C_{sites}^{\ }\  = \frac{\rho_{sites}A_{specific}MW_{m}\phi_{m}}{V_{m}}$$

The format for inputting the molar density of sites is (mol/m^2^):

    <Surface complex> Value

#### Exchanger Concentrations

The concentration of ions on exchangers can be set in several ways, with
the most important distinction being between 1) exchange on the bulk
material, and 2) exchange on specific minerals. In the case of exchange
on a specific mineral, the code will link the concentration of exchange
sites (i.e., the cation exchange capacity, or CEC) to the concentration
of that mineral. Thus, an increase in the volume fraction of a secondary
phase (e.g., clay as a result of chemical weathering or Fe-oxyhydroxide
as a result of corrosion) will increase the cation exchange capacity of
the bulk material.

*Exchange on bulk material*: This option is used where no specific
mineral has been identified as the exchanger phase using the *on
mineral* syntax in the **ION_EXCHANGE** keyword block. To calculate
total concentration of exchange sites per bulk porous medium, the cation
exchange capacity represented by that site needs to be specified along
with the density and the total volume fraction of the bulk solid phase.
This information must be provided in each geochemical

To specify the cation exchange capacity associated with each exchanger
in a geochemical condition, the following syntax is used:

    Exchanger -cec Value

where the CEC is given in units of equivalents per gram solid. The name
of the exchanger must correspond to one of those identified in the
**ION_EXCHANGE** keyword block.

Solid density may be specified one of three ways:

1.  Direct specification of the bulk solid density (kg/m^3^)

>     SolidDensity Value

2.  Specification of the solid:solution ratio (g/kg water \~ g/L):

>     SolidDensity -ss_ratio Value

3.  Calculation based on the sum of the mineral volume fractions:

>     SolidDensity CalculateFromMinerals

*Exchange on a specfic mineral:* Here, the total number of equivalents
of exchange sites is calculated from the combination of the mineral
volume fraction and the cation exchange capacity per gram of that
mineral. The solid density, therefore, is not relevant (the keyword is
therefore ignored), although the values of the porosity and liquid
saturation affect the result (see below). Mineral densities are
calculated from the combination of molar volumes and molecular weights
in the database. As in the case of exchange on the bulk material, the
CEC is specified with a similar format, but here the units are
equivalents per gram mineral:

    Exchanger -cec Value