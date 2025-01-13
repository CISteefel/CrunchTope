
## RUNTIME

The **RUNTIME** keyword block contains a number of parameters for
running a particular simulation. All of them are optional, as is the
**RUNTIME** keyword block itself, since defaults values will be used if
none are provided by the user. All of the parameters, with the exception
of the species name given in the keyword parameter **master**, are case
insensitive. Within the **RUNTIME** keyword block, the following keyword
parameters are possible.

- **coordinate**

- **correction_max**

- **courant_number**

- **database**

- **database_sweep**

- **debye-huckel**

- **density_module**

- **dissolution_max**

- **fix_saturation**

- **generic_rates**

- **gimrt**

- **gimrt_pc**

- **gimrt_pclevel**

- **gimrt_solver**

- **graphics**

- **hindmarsh**

- **lag_activity**

- **later_inputfiles**

- **master**

- **pc**

- **pc_level**

- **precipitation_max**

- **reaction_path**

- **read_saturationfile**

- **restart**

- **save_restart**

- **screen_output**

- **solver**

- **speciate_only**

- **steady-state**

- **timestep_max**

- **timestep_init**

- **time_tolerance**

These keyword parameters are described in detail below.

    Example:
    RUNTIME
    time_units years
    timestep_init 1.E-09
    timestep_max 0.01
    time_tolerance 0.05
    gimrt true
    debye-huckel true
    database_sweep false
    speciate_only false
    master H+
    screen_output 10
    END

#### Coordinate

Keyword followed by an option, either *rectangular, cylindrical,* or
*spherical*

Syntax: **coordinate**  option

    <u> Default </u>:  Rectangular

Explanation: Rectangular coordinates are the default, but the user may
also choose cylindrical or spherical coordinates.

#### **Correction_max**

Absolute magnitude of the maximum correction to the natural logarithm of
a component species in any one Newton step.

Syntax:  **correction_max** *value*

*value* is a positive real number

<u> Default </u>: 2.0

Explanation:  This is an optional parameter which controls
the absolute magnitude of the maximum allowable correction to a
component species in any particular Newton step. The correction is in
units of the natural logarithm of the species molality. The limitation
of the correction is a form of damping of the Newton step which prevents
overflows in many cases. Correction factors over 10 are not recommended.
Excessively small correction factors (\< 1) will slow the rate of
convergence of the Newton method significantly.

#### Courant_number

Maximum Courant or Courant-Friedrichs-Lewy number (defined by
$CFL = \frac{v\Delta t}{\Delta x}$) for OS3D option (GIMRT *false)*

Syntax:  **Courant_number** *value*

*value* is a positive real number between 0.0 and 1.0

<u> Default </u>: 0.5

Explanation:  This is an optional parameter that controls
the maximum CFL number by limiting the timestep, Δt, when the OS3D
option is used. The keyword is ignored if the GIMRT option is selected.
The Courant number cannot be greater than 1.0 for the OS3D option, since
results in an unstable method.

#### Database

Option to specify a database file.

Syntax:  **database** *filename*

*filename* is an optional string (up to 132 characters) indicating the
name of the database file.

<u> Default </u>:  datacom.dbs

Explanation:  This specification of the name of the
database file may also go in the DATABASE keyword block.

#### database_sweep

Keyword followed by *true* or *false* which selects whether to carry out
a full sweep of the database, loading all possible species and minerals.

Syntax: **database_sweep** *logical*

*logical* is a standard FORTRAN logical (true or false).

<u> Default </u>: false

Explanation: The keyword parameter database_sweep is followed by true or
false (or yes or no) and is used to select, when true, the option to
sweep the entire thermodynamic database so as to load all possible
species and minerals. In this case, the user's choice of secondary
species and minerals listed in the input file is ignored. After the
sweep, each of the geochemical conditions is speciated and the code
stops. User's who want to use these species and minerals in an actual
reactive transport run need to copy the list of secondary species and/or
minerals over to the keyword blocks for these categories within the
input file and then deselect the database_sweep option. In this mode,
the code runs essentially like EQ3 does in terms of loading all possible
relevant species given a choice of primary component species. This
option supersedes whatever is set with the speciate_only keyword
described below.

#### **Debye-Huckel**

Keyword followed by *true* or *false* which selects model for
calculation of activity coefficients.

Syntax: debye-huckel logical

logical is a standard FORTRAN logical (true or false).

<u> Default </u>: true

Explanation: The keyword parameter debye-huckel is followed by true or
false (or yes or no) and is used to select the model to compute activity
coefficients. Currently, the only two options are the extended
Debye-Huckel model (Helgeson) or activity coefficients equal to 1 (no
activity corrections).

#### **Density_module** {#density_module .SteefelHeading4}

Keyword followed by an option, either temperature,
sodium[\_]{.underline}nitrate, sodium_chloride, or potassium_nitrate.

Syntax: density_module option

<u> Default </u>:  Temperature

Explanation: The default is to calculate the fluid density only as a
function of temperature. Limited options are available to calculate the
fluid density as a function of composition, with sodium nitrate and
sodium chloride solution being the two most common.

#### **Dissolution_max** {#dissolution_max .SteefelHeading4}

Maximum decrease in mineral volume fraction for any one time step.

Syntax: dissolution_max value

value is a real positive number.

<u> Default </u>:  0.001 (m^3^ mineral/m^3^ porous medium)

Explanation: This keyword parameter specifies the maximum decrease in
mineral volume fraction (dimensionless) for any one time step. This
option is applied in calculating the size of the time step for the nex
time step based on rates computed in the preceding one. It is not,
therefore, used currently as a firm control on time step since no repeat
of a time step is carried out if dissolution_max is exceeded in a step.

#### fix_saturation

A constant, uniform liquid saturation may be specified with this
keyword.

Syntax: fix_saturation *value*

*value* is a positive real number between 0.0 and 1.0

<u> Default </u>: 1.0

Explanation: This is an optional parameter that can be used to set a
constant, uniform liquid saturation for a problem. The value must be
greater than 0.0 and less than or equal to 1.0.

#### Generic_rates

Keyword followed by a value representing the logarithm of the reaction
rate to be used for all minerals in the MINERAL keyword block. This
option avoids the need for kinetic entries for the minerals in the
database

Syntax:  **generic_rates** *value*

*value* represents the logarithm of the mineral rate in units of
mol/m^2^/sec.

<u> Default </u>:  The option is not selected and individual kinetic
entries will be required in the database.

Explanation:  The keyword parameter **generic_rates** is
followed by a value for the logarithm of the mineral reaction rate that
will be used for all of the minerals given in the MINERAL keyword block.
In this mode, the code will not require that a kinetic rate law be
explicitly given in the database for each of the minerals considered.
This option may be useful where the interest is in simply tracking a
reaction path without detailed kinetic information.

#### Gimrt

Keyword followed by *true* or *false* which selects coupling method for
reaction and transport.

Syntax:  **gimrt** *logical*

*logical* is a standard FORTRAN logical (true or false).

<u> Default </u>: true

Explanation:  The keyword parameter **gimrt** is followed
by *true* or *false* (or *yes* or *no*) and is used to select a global
implicit coupling of reaction and transport (GIMRT) if true, or time
splitting of reaction and transport (OS3D) if false. Other restrictions
on the coupling method may apply---for example, a global implicit method
of coupling reaction and transport is currently allowed only up to 2
dimensions. The selection of a global implicit method (**gimrt** *true*)
eliminates any strict requirement that a CFL criteria be observed,
although other practical limitations on the time step (either in terms
of numerical stability or of time truncation error) may apply. See
Steefel and MacQuarrie (1996) for a fuller discussion of the pros and
cons of the two coupling methods.

#### Gimrt_pc

Keyword followed by an option, either *bjacobi,* or *ilu.*

Syntax:  **gimrt_pc** *option*

<u> Default </u>:  *bjacobi*

Explanation:  This sets the preconditioner method for the
GIMRT option, with the choices being a block Jacobi method (bjacobi) or
ILU. See the PETSc web pages at
[www.anl.gov/petsc](http://www.anl.gov/petsc) for further explanation.

#### Gimrt_pclevel

Keyword followed by an integer value giving the level of fill to be used
in preconditioning the global implicit reactive transport matrix*.*

Syntax:  **gimrt_pclevel** *integer value*

<u> Default </u>: *2*

Explanation:  This sets the preconditioner level of fill
for the GIMRT option as an integer value. Using higher levels of fill
than 2 may slow the speed of execution. See the PETSc web pages at
[www.anl.gov/petsc](http://www.anl.gov/petsc) for further explanation.

#### Gimrt_solver

Keyword followed by an option, either *gmres* or *bcgs.*

Syntax:  **gimrt_solver** *option*

<u> Default </u>: *gmres*

Explanation:  This sets the solver method for the GIMRT
option, with the choices being the GMRES method (*gmres*) or a stablized
biconjugate gradient (*bcgs*). See the PETSc web pages at
[www.anl.gov/petsc](http://www.anl.gov/petsc) for further explanation.

#### Graphics

Keyword followed by an option, either *kaleidagraph, tecplot,* or
*xmgr.*

Syntax:  **graphics** *option*

<u> Default </u>:  *kaleidagraph* for 1D, *tecplot* for 2-3D.

Explanation:  This sets the style of output for the various
graphics packages, with *Kaleidagraph* being the default for 1D,
*Tecplot* the default for 2D. The *xmgr* option can be used for either
the X windows based <u> Xmgr </u>, or for <u> Gnuplot </u>.

#### Hindmarsh

Keyword followed by *true* or *false* which selects the linear solver
for one-dimensional problems.

Syntax:  **hindmarsh** *logical*

*logical* is a standard FORTRAN logical (true or false).

<u> Default </u>: true

Explanation:  The keyword parameter **hindmarsh** is
followed by *true* or *false* (or *yes* or *no*) and is used to select
the linear solver for one-dimensional problems. Setting **hindmarsh** as
true turns on the direct block tridiagonal solver written by Alan
Hindmarsh (1977). Setting **hindmarsh** to false turns on the PETSc
solver.

#### Lag_activity

Keyword followed by *true* or *false* which selects when activity
coefficients are updated.

Syntax:  **lag_activity** *logical*

*logical* is a standard FORTRAN logical (true or false).

<u> Default </u>: true

Explanation:  The keyword parameter **lag_activity,** when
true, causes the code to update activity coefficients only at the
beginning of the time step, outside of the Newton iteration loop. In
this scheme, activity coefficients are based on the ionic strength
calculated from the previous time step. Since at this time ionic
strength is not included as an independent unknown in CrunchFlow, the
method yields accurate results in most cases and shows better
convergence behavior than does a scheme in which activity coefficients
are updated throughout a particular Newton iteration loop. The method
may not be completely accurate when large contrasts in ionic strength
occur in a problem, although even this depends on the time step used.

#### Later_inputfiles

Option to specify additional input files to be read on restart.

Syntax:  **later_inputfiles**
*\[filename1\],\[filename2\],\[filename3\]...*

*filename1, filename2...* are strings (up to 132 characters). indicating
the name of the input files to be read following restart.

<u> Default </u>:  No additional input files are read.

Explanation:  This option allows additional input files to
be read automatically restart. Normally this option is used in
conjunction with the keywords *save_restart* and *restart*, which
generate and read respectively restart of "pickup" files that can be
used to read in the spatial distribution of species and minerals from a
previous run, changing a boundary condition. There is no a a priori
restriction on the number of files. The *later_inputfiles* option should
be specified only in the first input file in the series to be read
(later specifications will be ignored).

#### Master

Master variable to be used for time step control.

Syntax:  **master** *name*

*name* is a string up to 132 characters representing the name of a
species.

<u> Default </u>: *O2(aq),* if present, otherwise *H^+^* if
present, otherwise species number 1.

Explanation:  The keyword parameter time_tolerance is
applied to this species. Normally, for redox problems *O2(aq)* is used
since it is the most sensitive "master variable". For non-redox
multicomponent problems, *H^+^* is normally used. The name of the
species is case-sensitive, unlike other parameters specified within the
**RUNTIME** keyword block.

#### OvershootTolerance

Keyword followed by a real number giving the tolerance for overshooting
a mineral volume fraction.

Syntax:  **OvershootTolerance** *value*

*value* is a positive real number

<u> Default </u>: *1.0E-05*

Explanation:  This keyword sets the tolerance for
overshooting a zero volume fraction. This is necessary because of the
way in which CrunchFlow handles mineral abundances. For any one time
step, mineral abundances and surface areas are assumed fixed---mineral
volume fractions, therefore, are updated with an explicit scheme at the
end of the time step. This approach, however, makes it possible to
overshoot the available mineral concentration, which can result in
overall mass balance errors. The tolerance set here specifies the
maximum overshoot in units of m^3^ mineral/m^3^ porous medium. If the
tolerance is exceeded, the time increment, *delt*, is cut and the time
step is repeated.

#### Pc

Keyword followed by an option, either *ilu, jacobi,* or *direct.*

Syntax:  **pc** *option*

<u> Default </u>: *ilu*

Explanation:  This sets the preconditioner method for PETSc
to solve the flow or diffusion equations. It does not apply to the
solution of the global implicit reactive transport equations in the
GIMRT option. The choices include *ilu*, with a partial fill set by the
keyword *pclevel*, *jacobi*, and *direct*. The *direct* option carries
out a direct solution of the matrix without iteration (i.e., it is not a
sparse matrix solver). See the PETSc web pages at
[www.anl.gov/petsc](http://www.anl.gov/petsc) for further explanation.

#### Pclevel

Keyword followed by an integer value giving the level of fill to be used
in preconditioning the flow or diffusion matrices*.*

Syntax:  **pclevel** *integer value*

<u> Default </u>: *5*

Explanation:  This sets the preconditioner level of fill
for solving the flow or diffusion equation. An integer value must be
provided. Normally, levels of fill of above 5 should only be used in the
case of an ill-conditioned matrix, as typically occurs with large
differences in permeability. Use of higher fill levels will slow
execution in the case of reasonably well-conditioned linear systems. See
the PETSc web pages at [www.anl.gov/petsc](http://www.anl.gov/petsc) for
further explanation.

#### Reaction_path

Keyword followed by *true* or *false* which selects whether or not to
run in "reaction path" mode as opposed to "batch" mode.

Syntax:  **reaction_path** *logical*

*logical* is a standard FORTRAN logical (true or false).

<u> Default </u>: false

Explanation:  The keyword parameter **reaction_path** is
followed by *true* or *false* (or *yes* or *no*) and is used to select,
when true, the option to run in "reaction path" mode. This mode is
contrasted with "batch" mode and applies only in 0 zero spatial
dimension problems (a single grid cell is used). In reaction path mode,
mineral volume fractions are not updated. This is meant to represent
(loosely) the case of pure advective transport in which secondary
mineral products are left behind in a flow path.

#### Read_SaturationFile

Option to specify a file with liquid saturations to be read.

Syntax:  **read_SaturationFile** *filename format*

*filename* gives the name of the file(up to 132 characters) contining
the liquid saturation distributed over the entire spatial domain.

<u> Default </u>: None.

Explanation:  This provides the name of a file containing
liquid saturations distributed in space to be read. This option
overrides the *fix_saturation* keyword. The format of the file depends
on the optional format specified (see Section 7.1.2). If no format is
specified, the code assumes a single column (SingleColumn) of values in
the order 1-NX, 1-NY, 1-NZ:

>     DO jz = 1,nz
>       DO jy = 1,ny
>         DO jx = 1,nx
>           READ(52,*) satliq(jx,jy,jz)
>         END DO
>       END DO
>     END D

#### Restart

Option to use a restart file.

Syntax:  **restart** *\[filename\] \[append\]*

*filename* is an optional string (up to 132 characters). indicating the
name of the restart file.

<u> Default </u>:  crunch.rst

Explanation:  This keyword parameter specifies that the
code should be restarted. If restart is not followed by anything, the
default restart file name of *crunch.rst* is used. Alternatively, the
user may specify an optional restart file name (this can be useful,
since the default file name will be overwritten with a new run). This
optional restart file, however, must have been created in an earlier run
with the *save_restart* keyword. Binary restart files will be created
each time a set of spatial profiles are written (see keyword
*spatial_profile* in the **OUTPUT** keyword block). With a normal
successful termination of the code, the restart file will be updated at
the end of the run. Subsequent output to the "breakthrough" files (see
keyword *time_series*) can be appended by adding the optional *append*
keyword after the filename. In the same way, the file counter "n" (for
example, pHn.dat) will be updated as well if this *append* option is
chose, so spatial profiles written after restart will reflect the files
written previously before restart.

#### Save_restart

Option to save a restart file.

[Syntax]{.underline}*: save\_***restart** *\[filename\]*

*filename* is an optional string (up to 132 characters). indicating the
name of the restart file to be written to.

<u> Default </u>:  *crunch.rst*

Explanation:  This keyword parameter specifies that the
code should create a restart file with the name provided. If the keyword
*save_restart* is not provided, restart information will be saved to the
default restart file name of *crunch.rst*.

#### Screen_output

Keyword followed by an integer value giving the interval at which run
time output is written to the screen.

Syntax:  **screen_output** *integer value*

<u> Default </u>: *1*

Explanation:  This gives the interval at which run time
output (time, time step size, number of Newton iterations required...)
is written to the screen. This is useful where the code is taking a
large number of time steps to complete a run and only periodic
information is needed. This number should be a minimum of 1 (output
every time step).

#### Solver

Keyword followed by an option, either *bicg, gmres, bcgs,* or *direct.*

Syntax:  **solver** *option*

<u> Default </u>: *bicg*

Explanation:  This sets the solver method for PETSc to
solve the flow or diffusion equations. It does not apply to the solution
of the global implicit reactive transport equations in the GIMRT option.
The choices include a bi-conjugate gradient solver, *bicg*, and the
additional options of *gmres, bcgs*, and *direct*. The *direct* option
carries out a direct solution of the matrix without iteration (i.e., it
is not a sparse matrix solver). See the PETSc web pages at
[www.anl.gov/petsc](http://www.anl.gov/petsc) for further explanation.

#### Speciate_only

Keyword followed by *true* or *false* which selects whether or not to
speciate the geochemical conditions and then stop.

Syntax:  **speciate_only** *logical*

*logical* is a standard Fortran logical (true or false).

<u> Default </u>: false

Explanation:  The keyword parameter **speciate_only** is
followed by *true* or *false* (or *yes* or *no*) and is used to select,
when true, the option to speciate the geochemical conditions and then to
stop without carrying out any reactive transport calculations. This
option is contrasted with database_sweep which, in addition to stopping
after the speciation calculations, also sweeps the database to load all
of the possible species and minerals. Using speciate_only, the user's
list of secondary species and minerals given in the input file is used.

#### Steady_state

Keyword followed by *true* or *false* which selects whether or not to
run to a quasi-steady state as a condition for stopping a particular
simulation.

Syntax:  **steady_state** *logical* \[*tolerance*\]

*logical* is a standard FORTRAN logical (true or false)

*tolerance* is an optional real number giving the criteria for
attainment of a steady state.

<u> Default </u>:  false, with a tolerance of 1.E-07 in units of
mol/kg/year normalized to the primary species total concentration (so
tolerance has units of 1/year).

Explanation:  The keyword parameter **steady_state** is
followed by *true* or *false* (or *yes* or *no*) and is used to select,
when true, the option to run to a quasi-steady state and then terminate
the program rather than to rely on the specific output times given in
the **OUTPUT** keyword block. If true, the logical specifier can be
followed by an optional convergence criteria in units of normalized
mol/kg/year for all of the primary species total concentrations (the
normalization then results in a convergence criteria with units of
1/year). Depending on the value of the convergence criteria (a value of
1.E-07 is the default), a quasi-steady state may be attained since in
many cases slow dissolution of primary mineral phases may prevent a true
steady-state from being achieved. In many cases, however, it is possible
to attain a "quasi-stationary state" (see Lichtner, 1988) in which
aqueous concentrations are nearly unchanging while mineral volume
fractions are slowly evolving.

#### Timestep_max

Indicates maximum time step allowed in problem.

Syntax:  **timestep_max** *value*

*value* is a positive real number

<u> Default </u>: 1.0 year

Explanation:  The default units are years, with the option
to reset the time units to one of those discussed under **Time Units**
(e.g., seconds, days, etc.) using the **time_units** parameter keyword.
This is the maximum allowed time step, but a smaller maximum may be used
if dictated by Courant number restrictions, as in the case of CrunchFlow
run in OS3D mode (time splitting of transport and reaction).

#### Timestep_init

Indicates initial time step to be used in simulation..

Syntax:  **timestep_init** *value*

*value* is a positive real number

<u> Default </u>: 1.E-09 years

Explanation:  The default units are years, with the option
to reset the time units to one of those discussed under **Time Units**
(e.g., seconds, days, etc.) using the **time_units** parameter keyword.
The initial time step also serves as the minimum time step which is only
superseded by a restriction on the Courant-Friedrich-Lewy (CFL) number.

#### Time_tolerance

Tolerance for the second derivative of the change in the natural
logarithm of the master species with respect to time.

Syntax:  **time_tolerance** *value*

*value* is a positive real number

<u> Default </u>: 0.001

Explanation:  Tolerance applies to the second derivative of
the natural logarithm of the master species (specified by the keyword
parameter master). This parameter can be used to either control the time
truncation error or to limit the increase in time steps for the sake of
numerical stability. With a time tolerance set very high, the code will
normally increase the time step by a factor of 2 every time step until
the maximum time step is reached.

