## CrunchFlow

**Software for Modeling Multicomponent Reactive Flow and Transport**

*USER'S MANUAL*

**Carl I. Steefel**\
**Toshiyuki Bandai**\
**Sergi Molins**

Energy Geosciences Division\
Lawrence Berkeley National Laboratory\
Berkeley, CA 94720 USA

CISteefel@lbl.gov

***January 13, 2025***

### Features of CrunchFlow

Before running the code, the user of **GIMRT** needs to decide what
physical process(es) he or she wants to model and design the simulation
accordingly. This involves choosing the boundary conditions for the
problem (e.g., no flux boundaries or fixed concentration boundaries) and
setting up the initial grid and/or initial mineral zones. The user must
decide how many grid points are needed to solve a particular problem.
Typically one proceeds by carrying out the simulations on a relatively
coarse grid (e.g., 30 to 50 grid points) and only refining the grid
later when one has a general idea of the behavior of the system. The
user also needs to decide what level of chemical complexity they wish to
consider.

This can involve choices in the number of independent components
(primary species), the number of secondary species, and the number of
reacting minerals. The discussion in the sections which follow is
designed to make the basis for these decisions clearer. These sections
are then followed by a step by step discussion of the input file which
is needed to run **GIMRT**.

**NOTE:** &nbsp; OS3D temporarily disabled until water activity and ionic strength are
folded in as master variables.

### Units

Although there are some advantages to allowing the input of data in any
self-consistent set of units, many of the parameters used in **GIMRT**
are given commonly in certain units, so we have decided to require these
units. All units for input parameters (e.g., velocities, reaction rates
etc.) are given in the description of the input file. All distances are
in meters. The code assumes that the reaction rates for minerals are in
units of mol $m^{-2}$ s$^{-1}$ (these are converted in the code to
mol $m^{-2}$ yr$^{-1}$. Diffusion coefficients are in units of
$m^{-2}$ s$^{-1}$. All other references to time (time for plot
files, timesteps etc.) are in years. Velocities, for example, are in
units of m yr$^{-1}$. Concentrations for aqueous species are in units
of mol kg$^{-1}$ water (i.e., molality). The partial pressure of
gases are specified in bars.

### Boundary Conditions

At the present time there are three possible boundary conditions which
can be used in **GIMRT** and **OS3D**. The most straightforward way to
describe the possible boundary conditions in the codes is to actually
begin with the *third* type of boundary (also sometimes called a
Danckwerts boundary conditon). This is a boundary condition which
requires continuity of flux across the boundary. Both codes assume that
the flux at these boundaries is purely advective (i.e., no diffusion or
dispersion occurs upstream). Where the flow is into the system, the flux
at the first interior nodes downstream is assumed to be just $uC_{ext}$,
the same as the advective flux at the upstream, exterior node. In this
special case where the dispersive and diffusive flux are neglected, the
*third* type of boundary condition is the same as a Dirichlet or fixed
concentration boundary condition. At the boundary nodes where the flow
passes out of the domain the flux is also assumed to be purely
advective, so that it passes out of the system with the downstream,
exterior node having no effect on the system. Note that when the flow
across a boundary is non-zero, we always assume that the dispersive and
diffusive fluxes are equal to zero.

The second type of boundary condition is the no-flux condition which is
applied when there is no flow across the boundary. In this case, both
the advective flux and the dispersive/diffusive fluxes are equal to
zero.

The last type (often called the *first* type or a Dirichlet boundary
condition) is normally applied when the advective flux is equal to zero,
but there is a dispersive or diffusive flux across the boundary (if the
advective flux is non-zero, then the code automatically assumes a zero
dispersive/diffusive flux). This condition is many cases is not
physically reasonable, since it assumes that the diffusive and/or
dispersive flux through the boundary does not affect the concentration
of the exterior node (i.e., it is fixed). This boundary condition would
normally be applied where the physical processes exterior to the domain
are such that the diffusive/dispersive flux through the boundary is not
enough to alter the concentrations there. One example might be a
boundary along a fracture within which the flow is so rapid that the
diffusive flux through the fracture wall has a negligible effect (e.g.,
Steefel and Lichtner, 1994). Another example might be a boundary
corresponding to an air-water interface where a gas like O$_{2}$ or
CO$_{2}$ is present in such large abundance that it effectively
buffers the concentrations in the aqueous phase. Another example which
is often mentioned in this regard is the case where mineral equilibria
fixes the concentrations. This is dangerous, however, since the system
is only completely determined in the special case where the phase
equilibria results in an invariant point, i.e., no thermodynamic degrees
of freedom exist. This is extremely rare in open systems, despite the
frequent invocation of this condition in metamorphic petrology.
Normally, since the concentrations will not be uniquely determined, a
rigorous treatment would require that a speciation calculation with
equilibrium constraints be carried out to actually determine what the
concentrations at the exterior node should be. This, however, is not
implemented at the present time.

**OS3D** and **GIMRT** differ slightly in their application of a
Dirichlet boundary condition where the advective flux across the
boundary is zero. In the case of **GIMRT**, the user may specify that an
exterior node be a fixed concentration node and that it affect the
interior node (i.e., there will be a non-zero diffusive and/or
dispersive flux through the boundary). In **OS3D**, however, the code
*always* assumes the dispersive/diffusive flux across the boundary is
zero (due to the way the transport algorithm is formulated). If one
wants to use a Dirichlet condition in the case of **OS3D**, one does so
by fixing the concentrations at the first node or nodes interior to the
boundary. This will have the same effect as allowing an exterior node
influence the interior domain, except that one can think of the boundary
as having been moved inward by one grid cell. The procedure for doing
this is described in the section below explaining the input file
threedin.

To summarize the boundary condition options in **OS3D** and **GIMRT**, a
non-zero flow across the boundary means that the code will use the
upstream, exterior boundary concentrations to give a purely advective
flux at the boundary and the code will allow a fluid packet to pass out
of the system with no influence from a downstream, exterior nodes.
Exterior nodes across a boundary through which there is no flow will not
influence the interior domain, unless the Dirichlet boundary condition
option is used (usage described below). The default is therefore no
diffusive/dispersive flux boundaries which is only overridden where the
Dirichlet option is invoked.

### Choice of Primary and Secondary Species

After setting up the problem physically, one needs to decide on the
chemistry to be included. **GIMRT** requires that one specify an initial
choice of primary species which then determine the number of independent
chemical components in the system. This choice of primary species is not
unique. The user must choose,
however, which secondary species are to be included in the simulation.
This can be an advantage for multicomponent reaction--transport modeling
since one would like the option of carrying out chemically simplified
simulations. However, the user must take care that important species are
not neglected. A sweep of the database can be carried out with the 
*database_sweep* option.

The only cautionary note here is that one must include either among the
primary or secondary species a species which is used in writing the
reaction in the **EQ3/EQ6** database. Note that this database is
provided with the code in a reformatted form. For example, if one wants
to include aluminum in the simulation, then the species $Al^{+ 3}$ must
appear among *either* the primary or secondary species since that
species is used as the primary or basis species in the **EQ3/EQ6** data
base. In the same way, one must include $O_{2(g)}$ among the list of
gases if $O_{2(aq)}$ is to be present in the system since much of the
**EQ3** database is written in terms of $O_{2(g)}$. If this is not done,
the code will not be able to find the reaction when searching the data
base. One should not, however, include any secondary species which
cannot be written completely in terms of the basis species specified in
the input file. The simplest way to grasp these features is to check the
database to see which species are used to write the reactions involving
a species of interest.

With the recent release of the software, it is necessary to include $H_2O$ 
as a primary species (i.e., it becomes one of the primary unknowns). The 
code will calculate the activity coefficient of $H_2O$ according to

$$ \gamma_{H_2O} = 1.0 - 0.17 \sum m_j $$

where the $m_j$ are the molalities of the species primary and secondary 
species in the simulation.

### Redox Reactions

There are many possibilities for choices of primary and secondary
species to represent redox reactions. **GIMRT** writes all redox
reactions as whole cell reactions, i.e., electrons are *not* used as
either primary or secondary species. It is easy to show that the
electrons are completely unnecessary if the solutions are electrically
neutral, and the inclusion of electrons requires that one include a
*species* the mass of which should everywhere be zero. There are cases
where electrons need to be included, as in corrosion problems where an
actual electron flux occurs. When the solution is electrically neutral,
however, no special provision for redox species and reactions is
necessary when whole cell reactions are used (in contrast to statements
by Yeh and Tripathi (1989)), although the inclusion of redox reactions
implies that we have included an additional balance equation
representing the conservation of exchangeable charge.

Since $O_{2(g)}$ is used as the basis species in redox reactions in the
**EQ3/EQ6** data base, it must be included in the list of gases in the
input file if redox reactions are to be considered. A gas, however,
cannot currently be used as a primary species in **GIMRT**. One could 
use $O_{2(aq)}$, however, or any number of other species 
that form a redox couple, for
example, $Fe^{+ 2}$ and $Fe^{+ 3}$, or $HS^{-}$ and $SO_{4}^{- 2}$.

### Acid- Base Reactions

Acid--base reactions require no special treatment other than in the
initialization process where one normally fixes or calculates pH rather
than an $H^{+}$ concentration. Once a pH has been either fixed or
calculated in the initialization procedure, however, $H^{+}$ can be
treated like any other component. Again, no special considerations are
needed in contrast to the statements by Yeh and Tripathi (1989).

### Activity Coefficient Model

At this time, the user may choose to run the code either with unit
activity coefficients or with an extended Debye--Hückel formulation to
calculate activity coefficients for the ionic species. The Debye--Hückel
formulation in its simplest form is given by

$$\log\gamma_{i} = - \frac{AZ_{i}^{2}(I)^{1/2}}{1 + B{\overset{\circ}{a}\ }_{i}(I)^{1/2}} + \dot{b}I$$

The coefficients for the temperature dependence of the parameters are
taken from the **EQ3** data base.

### Initialization Procedure

**GIMRT** initializes the solute concentrations at individual grid
points over the domain by carrying out a distribution of species using a
number of options. The initial concentrations of the primary species may
be determined by

-  carrying out a mass balance on a total concentration (e.g.,
  aluminum),

-  carrying out a charge balance on a species (e.g., Cl$\ ^{-}$),

-  a mineral equilibrium constraint (e.g., $H^{+}$ required to be in
  equilibrium with calcite),

-  fixing the partial pressure of a gas (e.g., CO$\ _{2(g)}$),

-  fixing the activity (e.g., pH) of a species,

-  fixing the concentration of a species, and

-  fixing the number of moles surface hydroxyl sites per mole mineral
  and the intrinsic surface area of the mineral (i.e., m$\ ^{2}$/g).

This last option is available in **GIMRT** and it allows one to specify
the concentration of surface complexes. These are specified by giving a
site density (mole per $m^{2}$ mineral) and a surface area of mineral
per gram mineral. The total concentration of surface hydroxyls developed
on the mineral (converted to mol/kg solvent) is then computed from this
information and the abundance of the mineral specified in the startup
file. If a mineral is initially is not present but one would like to
consider adsorption on this mineral, then the total concentration of
surface hydroxyls on this mineral is set to a small number
(10$^{-10}$ mol/kg) and its concentration only builds up once the
secondary mineral precipitates.

How each of the above initialization constraints is applied is described
in the section on the input file below.

### Temperature Gradients

Since **GIMRT** is based on finite difference methods, it has no problem
handling temperatures which are non--constant either in time or in
space. If either a non--default thermodynamic data base or the
master25.data data base is not used at startup, then the code defaults
to the mastertemp.data file where equilibrium constants are given at
temperatures of 0$^{\circ}$C, 25$^{\circ}$C, 60$^{\circ}$C,
100$^{\circ}$C, 150$^{\circ}$C, 200$^{\circ}$C, 250$^{\circ}$C,
and 300$^{\circ}$C at pressures corresponding to the water saturation
curve.  The code then fits a polynomial to the thermodynamic data so that
equilibrium constants can be generated at any temperature between
0$^{\circ}$C and 300$^{\circ}$C along the water saturation curve.
In addition, the algorithm for fluid
density must be modified as well. Presently the fluid density is
calculated as a function of temperature from a polynomial fit to the
data in Helgeson and Kirkham (1974). Fluid viscosity is calculated from
the data presented by Bruges and others (1966).

Custom databases can be developed at higher pressures using the 
software *pygcc* developed by Awolayo and Tutolo (2022).

A temperature dependence is also included in the calculation of the
reaction rate constants, the diffusion coefficients, and
the fluid densities. Reaction rate constants at temperature are obtained
normally by extrapolating the rate constants for 25$^{\circ}$C given
in the input file. If one would rather input directly the desired rate
constant at temperature, than the simplest procedure (short of changing
the source code) is to input the rate *as if it were the
25*$^{\circ}$*C data* while changing the activation energy to 0
kcal/mole. There is also an option for directly inputting the diffusion
coefficient at temperature.

There are two options for temperature presently implemented in
**GIMRT**. The first is simply to run the simulation isothermally,
either at 25$^{\circ}$C or at temperature. The second is to specify a
linear temperature gradient which operates only in the X direction. Many
other ways of handling temperature are possible, including allowing for
a transient temperature field, but this must be provided by the
user/programmer at the present time. This can be done with a simple read
of temperatures, making sure that the dimensions of the temperature
field match those specified in the code. The positions in the code where
a user should place a routine for calculating or inputting a temperature
are marked in the routine gimrt.f.

### Running CrunchFlow

CrunchFlow reads a user-provided input file on startup which provides
the necessary physical and chemical parameters needed for a run. The
input file, the name of which is specified by the user, is keyword-based
so that the order of appearance does not matter. In many cases, it is
possible to leave various optional parameters out altogether and thus
allowing the code to use default parameters. Keywords are grouped
broadly into keyword blocks**,** which in turn include a variety of
keyword parameters. Keyword blocks may appear anywhere in the input file
and keyword parameters may appear anywhere within a particular block,
but certain keyword parameters must appear in the appropriate blocks.
The possible keyword blocks include:

**TITLE
RUNTIME
OUTPUT
PRIMARY_SPECIES
SECONDARY_SPECIES
GASES
MINERALS
AQUEOUS_KINETICS
ION_EXCHANGE
SURFACE_COMPLEXATION
PEST
DISCRETIZATION
INITIAL_CONDITIONS
BOUNDARY_CONDITIONS
TRANSPORT
FLOW
POROSITY
TEMPERATURE
CONDITION**

With the exception of the keyword block **CONDITION**, the blocks should
appear only once in the input file. Each keyword block is terminated by
an **END**. The keyword block **CONDITION** is a special case in that it
can occur multiple times. Each occurrence of **CONDITION** specifies a
separate geochemical condition **(**these may be boundary or initial
conditions or source terms**)** containing the geochemical input needed
to describe a particular problem. The various keyword blocks will be
discussed individually in more detail below.

### Reading Input File Name from a File

Normally, CrunchFlow will prompt the user to enter the name of an input
file interactively. The user can provide the name of the input file,
however, by including it in a file names *PestControl.ant*. The code
will check to see if this file exists in the directly from which
CrunchFlow is invoked, and if it does, will attempt to read the file
name from it. Upon successfully reading the input file name, the code
will then check to see that this file exists before reading from it. The
use of the *PestControl.ant* to provide an input file name will then
skip the requirement of interactive input from the user, an option that
is particularly useful when running PEST.

### Keyword Blocks

Each keyword block is initiated with the appropriate keyword **(**case
insensitive**)** on a line of its own and terminated by **END**, also on
a line of its own. The contents of any keyword block, then, is anything
between the keyword block name and the **END** specifier. Blank lines
within the keyword block will be ignored. Lines beginning with **!** are
treated as comments and ignored.

### Space and Time Units

Keyword parameters involving space and time units are used in several of
the keyword blocks. These take the form:

    time_units    unit
    space_units   unit

where *unit* in the case of time may be years (the default), days,
hours, minutes, or seconds, and in the case of space may be meters (the
default), kilometers, centimeters, millimeters, and micrometers.
Non-metric units (e.g., feet) are not recognized. As with any other
keyword parameter, they may occur anywhere within a particular block to
which they apply. A specification of a space or time unit, however, is
not global and applies only to the keyword block within which it occurs.

Table 1: List of possible space units.
  | Space Units  | Alternate Acceptable Specifications |
  | -----------------| ------------------|
  | meters (default) | meter, m |
  | kilometers       | kilometer, km |
  | centimeters      | centimeter, cm |
  | millimeters      | millimeter, mm |
  | micrometers      | micrometer, micron, microns, um |

Table 2: List of time units
  | Space Units  | Alternate Acceptable Specifications |
  | -----------------|---------------------------------| 
 | years (default)   | year                            |
 | days              | day                             |
 | hours             | hour                            |
 | minutes           | minute                          |
 | seconds           | second                          |
