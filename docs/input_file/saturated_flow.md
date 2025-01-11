### FLOW {#flow .SteefelHeading3}

#### **Constant_flow** {#constant_flow .SteefelHeading4}

Keyword followed by a value(s) giving the X, Y, and Z components of the
Darcy flux.

Syntax:  **constant_flow** *qx \<qy\> \<qz\>*

<u> Default </u>:  *0.0 (m^3^/m^2^/yr)*

Explanation:  This keyword is used to set a constant Darcy
flux. As a vector, it includes X, Y, and Z components. In a
one-dimensional problem in X, the Y and Z velocities are optional. Units
are changed inside the FLOW block with the *space_units* and
*time_units* keywords.

#### **Constant_gasflow** {#constant_gasflow .SteefelHeading4}

Keyword followed by a value(s) giving the X, Y, and Z components of the
gas flux.

Syntax:  **constant_gasflow** *qxgas \<qygas\> \<qzgas\>*

<u> Default </u>:  *0.0 (m^3^/m^2^/yr)*

Explanation:  This keyword is used to set a constant gas
flux. As a vector, it includes X, Y, and Z components. In a
one-dimensional problem in X, the Y and Z velocities are optional. Units
are changed inside the FLOW block with the *space_units* and
*time_units* keywords.

#### **Pump** {#pump .SteefelHeading4}

Keyword followed by a value(s) giving pumping rate, the name of the
geochemical condition to be used in the case of an injection well, and
its grid coordinates.

Syntax:  **pump rate label** *jx \<jy\> \<jz\>*

<u> Default </u>: *None*

Explanation:  This keyword beginning with *pump* sets a
pumping rate (positive for injection, negative for withdrawal) in units
of l/sec and assigns a geochemical condition to it. This is followed by
the coordinates of the well, with JY and JZ being optional in the case
of a one-dimensional problem. The geochemical condition is identified by
its label, although the condition is only used in the case of an
injection well.

#### **Gaspump** {#gaspump .SteefelHeading4}

Keyword followed by a value(s) giving gas pumping rate, the name of the
geochemical condition to be used in the case of an injection well, and
its grid coordinates.

Syntax:  **gaspump rate label** *jx \<jy\> \<jz\>*

<u> Default </u>: *None*

Explanation:  This keyword beginning with *gaspump* sets a
gas pumping rate (positive for injection, negative for withdrawal) in
units of l/sec and assigns a geochemical condition to it. This is
followed by the coordinates of the well, with JY and JZ being optional
in the case of a one-dimensional problem. The geochemical condition is
identified by its label, although the condition is only used in the case
of an injection well.

#### **Infiltration** {#infiltration .SteefelHeading4}

Keyword followed by a value(s) giving the infiltration or recharge rate,
the name of the geochemical condition to be used, and the coordinate
direction (X, Y, or Z).

Syntax:  **infiltration rate label** *coordinate*

<u> Default </u>: *None*

Explanation:  This keyword specifies an infiltration rate
which is distributed over the boundary of the system (at JX, JY, or JZ =
0). The geochemical condition to be used is given by its **label**.

#### **Read_VelocityFile** {#read_velocityfile .SteefelHeading4}

Keyword followed by a name for a file containing velocity vectors over
the entire spatial domain.

Syntax:  **read_VelocityFile** *filename format*

[Backwards compatible]{.underline}*:* **read_velocity**

<u> Default </u>: *None*

Explanation:  This keyword provides a filename prefix for
ASCII files containing velocity vector values over the entire spatial
domain. The extensions *vx*, *vy*, and *vz* are assumed to indicate the
files containing the X, Y, and Z velocities, respectively. For example,
the syntax

> [read_VelocityFile copperleach SingleColumn]{.underline}

would open and attempt to read the files:

> copperleach.vx\
> copperleach.vy\
> copperleach.vz

using the SingleColumn format (see Section 7.1.2). Units are set with
the *space_units* and *time_units* keywords---otherwise, default units
are m^3^/m^2^/yr.

The format of the file depends on the format chosen. In the absence of a
specification of the file format, a single column (*SingleColumn*)
format is assumed. Velocities are defined at grid interfaces (not at the
center of the nodes), so NX+1 X velocities are needed in the X
coordinate direction, NY+1 Y velocities in the Y coordinate direction,
and NZ+1 Z velocities in the Z coordinate direction. CrunchFlow uses the
convention that velocities for a particular grid cell are defined at the
right hand interface (the so-called "right hand rule"). So, for example,
the velocity, *qx(2,1,1)* is assumed to be the X velocity at the
interface between grid cell 2 and 3 (Figure XX). Accordingly, the first
entry in the *velocityfile.vx* file will be the X velocity defined at
the interface between grid cell (0,1,1) and (1,1,1). So, in the case of
where *ContinuousRead* is specified as the file format, the code would
read:

READ(52,\*) (qx(jx,jy,jz),jx=0,nx),jy=1,ny),jz=1,nz)

READ(53,\*) (qy(jx,jy,jz),jx=1,nx),jy=0,ny),jz=1,nz)

READ(54,\*) (qz(jx,jy,jz),jx=1,nx),jy=1,ny),jz=0,nz)

#### **Read_GasVelocityFile** {#read_gasvelocityfile .SteefelHeading4}

Keyword followed by a name for a file containing gas velocities over the
entire spatial domain.

Syntax:  **read_GasVelocityFile** *filename format*

[Backwards compatible]{.underline}*:* **read_gasvelocity**

<u> Default </u>: *None*

Explanation:  This keyword provides a filename prefix for
ASCII files containing gas velocities over the entire spatial domain.
The extensions *vx*, *vy*, and *vz* are assumed to indicate the files
containing the X, Y, and Z velocities, respectively. For example, the
syntax

> [read_GasVelocityFile co2stream ContinuousRead]{.underline}

would open and attempt to read the files:

> [co2stream]{.underline}.vx\
> [co2stream]{.underline}.vy\
> [co2stream]{.underline}.vz

using the ContinuousRead format (see Section 7.1.2) Units are set with
the space_units and time_units keywords---otherwise, default units are
m^3^/m^2^/yr.

The format of the file depends on the format chosen. In the absence of a
specification of the file format, a single column (*SingleColumn*)
format is assumed. Velocities are defined at grid interfaces (not at the
center of the nodes), so NX+1 X velocities are needed in the X
coordinate direction, NY+1 Y velocities in the Y coordinate direction,
and NZ+1 Z velocities in the Z coordinate direction. CrunchFlow uses the
convention that velocities for a particular grid cell are defined at the
right hand interface (the so-called "right hand rule"). So, for example,
the gas velocity, *qxgas(2,1,1)* is assumed to be the X velocity at the
interface between grid cell 2 and 3 (Figure XX). Accordingly, the first
entry in the gas*velocityfile.vx* file will be the X velocity defined at
the interface between grid cell (0,1,1) and (1,1,1). So, in the case of
where *ContinuousRead* is specified as the file format, the code would
read:

READ(52,\*) (qxgas(jx,jy,jz),jx=0,nx),jy=1,ny),jz=1,nz)

READ(52,\*) (qygas(jx,jy,jz),jx=1,nx),jy=0,ny),jz=1,nz)

READ(52,\*) (qzgas(jx,jy,jz),jx=1,nx),jy=1,ny),jz=0,nz)

#### **Calculate_flow** {#calculate_flow .SteefelHeading4}

Keyword followed by *true* or *false,* which indicates whether flow will
be calculated within the code or not.

Syntax:  **calculate_flow** *logical*

<u> Default </u>: *False*

Explanation:  This keyword indicates whether flow will be
calculated according to Darcy's Law or not. At this time, only fully
saturated flow calculations are possible. Once this logical is set to
*true*, the code will look for inputs on the pressure and permeability
fields.

#### **Read_PermeabilityFile** {#read_permeabilityfile .SteefelHeading4}

Keyword followed by a file name containing permeability values over the
entire spatial domain.

Syntax:  **read_PermeabilityFile** *filename format*

[Backwards compatible:]{.underline} **read_permeability**

<u> Default </u>: *None*

Explanation:  This keyword provides a filename prefix for
ASCII files containing permeability values over the entire spatial
domain. The extensions *x*, *y*, and *z* are assumed to indicate the
files containing the X, Y, and Z permeabilities, respectively. For
example, the syntax

> [read_PermeabilityFile copperleach SingleColumn]{.underline}

would open and attempt to read the files:

> copperleach.x\
> copperleach.y\
> copperleach.z

using the *SingleColumn* (see Section 7.1.2) format. Units of the
permeability are m^2^.

The format of the file depends on the file format specified. In none is
provided, a single column format consisting of a single permeability per
line, with NX varying first, then NY, and then NZ is assumed. Unlike the
fluid velocity, permeabilities are defined at the center of the nodes,
since at least in some cases they are linked to properties like the
porosity, which in turn may be affected by mineral dissolution and
precipitation. Interface permeabilities are calculated as harmonic means
between adjacent cells. For example, between nodes *i* and *i+1*, the
harmonic mean of the permeabililty is given by:

$k_{i + \frac{1}{2}} = \frac{2\left( k_{i}k_{i + 1} \right)}{\left( k_{i} + k_{i + 1} \right)}$.

Since permeabilities need to be defined at the boundaries to determine
the amount of flow across them, "ghost cell" X permeabilities are
required at the 0 and nx+1 nodes, and so on. For 1D problems, *permy*
and *permz* files are not needed.

For those familiar with FORTRAN, the actual source code in the case
where a *SingleColumn* format is specified is given by (compare to the
*ContinuousRead* format in the case of the *read_VelocityFile* keyword):

> DO jz = 1,nz\
> DO jy = 1,ny\
> DO jx = 0,nx+1\
> READ(52,\*) permx(jx,jy,jz)\
> END DO\
> END DO\
> END DO
>
> DO jz = 1,nz\
> DO jy = 0,ny+1\
> DO jx = 1,nx\
> READ(53,\*) permy(jx,jy,jz)\
> END DO\
> END DO\
> END DO
>
> DO jz = 0,nz+1\
> DO jy = 1,ny\
> DO jx = 1,nx\
> READ(54,\*) permz(jx,jy,jz)\
> END DO\
> END DO\
> END DO

#### **Permeability_x** {#permeability_x .SteefelHeading4}

Keyword used to specify the permeability in the X coordinate direction
by zone in the input file.

Syntax:  **permeability_x** *value zone jxbegin-jxend
jybegin-jyend jzbegin-jzend*

or **permeability_x** *value default*

<u> Default </u>: *None*

Explanation:  The *permeability_x* keyword uses an input
format similar to other parameters specified by zone (e.g., pressure,
tortuosity). The keyword is followed by the value of the X permeability,
which is then followed by either the label *zone* or *default*. If
*zone* is specified, then the code expects the X, Y, and Z coordinates
in the form of a beginning JX, and ending JX, a beginning JY and an
ending JY, and a beginning JZ and an ending JZ, in each case separated
by a hyphen. If *default* is specified instead of *zone*, the code will
initialize the entire spatial domain (i.e., all grid cells) to the
value. More than one specification of *permeability_x* can be (and
normally is) provided. For example, in a 2D system with 20 grid cells in
the X coordinate direction and 20 in the Y coordinate direction, a
higher X permeability zone running down the center of the domain could
be specified with

> [permeability_x 1.E-13]{.underline} default\
> [permeability_x]{.underline} [1.E-12]{.underline} zone 0-21 10-10 1-1

Note that values of the X permeability at grid cell 0 and NX+1 have been
specified, since the harmonic mean of the permeability is calculated at
the system boundaries using these "ghost cell" values. Zones specified
later in the sequence overwrite zones specified earlier---in the example
above, the default specification sets the *permeability_x* value to
1.E-13 m^2^ initially in all of the grid cells, while the next
*permeability_x* keyword creates a higher permeability zone in grid
cells JY=10 over the entire X length of the domain. To create a no-flow
boundary at JX=0 over the entire Y length, one could use:

> [permeability_x 1.E-13]{.underline} default\
> [permeability_x]{.underline} [1.E-12]{.underline} zone 0-0 1-20 1-1

or more selectively with something like:

> [permeability_x 1.E-13]{.underline} default\
> [permeability_x]{.underline} [1.E-12]{.underline} zone 0-0 1-8 1-1\
> [permeability_x]{.underline} [1.E-12]{.underline} zone 0-0 12-20 1-1

which would then create a no-flow boundary at JX=0 except for JY nodes
9, 10, and 11.

*The convention adopted here is that all three dimensions must be
specified as zones, even in a one-dimensional problem, although
permeability_y and permeability_z are not needed.*

#### **Permeability_y** {#permeability_y .SteefelHeading4}

Keyword used to specify the permeability in the Y coordinate direction
by zone in the input file.

Syntax:  **permeability_y** *value zone jxbegin-jxend
jybegin-jyend jzbegin-jzend*

or **permeability_y** *value default*

<u> Default </u>: *None*

Explanation:  The *permeability_y* keyword uses an input
format similar to other parameters specified by zone (e.g., pressure,
tortuosity). See the keyword *permeability_x* for details on how this
keyword is used to specify permeabilities by zone. *The convention
adopted here is that all three dimensions must be specified, even in a
one-dimensional problem.*

#### **Permeability_z** {#permeability_z .SteefelHeading4}

Keyword used to specify the permeability in the Z coordinate direction
by zone in the input file.

Syntax:  **permeability_z** *value zone jxbegin-jxend
jybegin-jyend jzbegin-jzend*

or **permeability_z** *value default*

<u> Default </u>: *None*

Explanation:  The *permeability_z* keyword uses an input
format similar to other parameters specified by zone (e.g., pressure,
tortuosity). See the keyword *permeability_x* for details on how this
keyword is used to specify permeabilities by zone. *The convention
adopted here is that all three dimensions must be specified, even in a
one-dimensional problem.*

#### **Gravity** {#gravity .SteefelHeading4}

Keyword used to indicate the coordinate direction of the gravity vector.

Syntax:  **gravity** *\<Angle to X \> \<Angle to Y \>
\<Angle to Z \>* *direction*

<u> Default </u>: *90°*

Explanation:  This keyword provides the angle relative to
the X, Y, and Z coordinate directions for the gravity vector. The
influence of gravity is then calculated from

$g = \cos\theta$[,]{.underline}

where θ is the angle of the gravity vector relative to the specific
coordinate direction (X, Y, or Z). Direction can be either *down* or
*up* and is used to indicate whether the direction the gravity vector is
assumed to operate (for example, toward X\[NX\] or toward X\[0\]). So,
for example, the following input

> gravity 90.0 0.0 90.0 down

would be used to set the gravity vector in the Y direction, with gravity
oriented down (toward JY = 0).

#### **Pressure** {#pressure .SteefelHeading4}

Keyword used to specify the fluid pressure by zone in the input file.

Syntax:  **pressure** *value zone jxbegin-jxend
jybegin-jyend jzbegin-jzend \[fix\]*

or **pressure** *value default*

<u> Default </u>: *None*

Explanation:  The *pressure* keyword uses an input format
similar to other parameters specified by zone (e.g., tortuosity,
permeability). The keyword is followed by the value of the fluid
pressure in units of Pascals, which is then followed by either the label
*zone* or *default*. If *zone* is specified, then the code expects the
X, Y, and Z coordinates in the form of a beginning JX, and ending JX, a
beginning JY and an ending JY, and a beginning JZ and an ending JZ, in
each case separated by a hyphen. If *default* is specified instead of
*zone*, the code will initialize the entire spatial domain (i.e., all
grid cells) to the value. More than one specification of *pressure* can
be (and normally is) provided, for example

> pressure 1000.0 default\
> pressure 0.0 zone 1-20 21-21 1-1 fix

Zones specified later in the sequence overwrite zones specified
earlier---in the example above, the default specification sets the fluid
pressure to 1000 Pa initially in all of the grid cells, while the next
keyword sets the pressure at JY = 21 to zero. In the case of a 20 by 20
X-Y simulation, the grid cell 21 corresponds to node NY+1 and is
therefore a ghost cell. This is the normal way in which fixed pressure
boundary conditions are set. Following the zone coordinates with *fix*
instructs the code to lock the value in for the course of the
simulation---not including *fix* would mean that the pressure is only
set as an initial condition, after which it is free to evolve according
to the physics of the problem. *The convention adopted here is that all
three dimensions must be specified, even in a one-dimensional problem.*

To create a hydrostatic pressure gradient in the Y direction, with depth
increasing with increasing JY, one could use:

> pressure 1000.0 default\
> pressure 0.0 zone 1-20 0-0 1-1 fix\
> gravity 90.0 0.0 90.0 up

With this input, the code would time step until the hydrostatic pressure
gradient was achieved in the absence of other forcings. To initialize
the system to hydrostatic pressure before time stepping begins, use:

> pressure 1000.0 default\
> pressure 0.0 zone 1-20 0-0 1-1 fix\
> gravity 90.0 0.0 90.0 up\
> initialize_hydrostatic true

This will create initially hydrostatic conditions even if pressure
boundary conditions or source/sink terms are such that hydrostatic
conditions will not prevail once time stepping begins.

#### **Initialize_hydrostatic** {#initialize_hydrostatic .SteefelHeading4}

Keyword followed by *true* or *false* which selects whether or not to
initialize the domain to hydrostatic conditions.

Syntax:  **initialize_hydrostatic** *logical*

*logical* is a standard Fortran logical (true or false).

<u> Default </u>: *False*

Explanation:  The keyword parameter
*initialize_hydrostatic* is followed by *true* or *false* (or *yes* or
*no*) and is used to select, when true, the option to initialize the
pressure field to hydrostatic conditions before time stepping begins.
When set to true, hydrostatic conditions will be initialized no matter
what other pressure boundary conditions or source/sink terms take effect
once time stepping begins. However, a coordinate direction for the
gravity vector must be selected with the *gravity* keyword and
appropriate boundary conditions (typically just a fixed pressure at one
end of the domain).