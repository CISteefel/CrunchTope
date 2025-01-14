
## INITIAL_CONDITIONS

In order to run a reactive transport or a reaction path calculations,
initial conditions must be set. Even in the case where the desired
result is a steady state simulation in which initial conditions are
irrelevant, the code must still achieve this steady state by stepping
through time. Only with the *speciate_only* or *database_sweep*
keywords set to true in the **RUNTIME** keyword block is it possible to
proceed without specifying initial conditions for the entire spatial
domain.

Required keyword block followed by a specification of the name of the
geochemical condition to be used, followed in turn by the grid cells
where the initial condition is to be set. The initial condition can be
fixed for the entire course of the simulation using the optional
appended keyword *fix*.

    Syntax:  condition_name  JX-JX [JY-JY] [JZ-JZ]  [fix]   

<u>Default:</u> &nbsp; None.  

<u>Explanation:</u>  &nbsp; Using this keyword block, geochemical
conditions are distributed over the spatial domain. Multiple
specifications can occur. Initial conditions listed later in the input
file will overwrite conditions specified above. The optional parameter
fix can be appended to the output so as to fix this geochemical
condition for the entire course of the simulation. The code now assumes
that if fix appears, it will occur immediately after the specification
of the X coordinates (*jx*) in a 1D problem, immediately after the Y
coordinates (*jy*) in a 2D problem, and after the Z coordinate (*jz*) in the
case of a 3D problem. In the example below, the code will verify that
the geochemical condition names (quartz_zone, portlandite_zone, and
nonreactive_zone) have been provided in the input file. This keyword
block, therefore, is read after the geochemical conditions specified in
the input file are read. The code will also verify that every grid cell
in the domain is specified, implying that this block is read after the
**DISCRETIZATION** keyword block.

*Example:*

    INITIAL_CONDITIONS
    quartz_zone 1-42 1-42
    portlandite_zone 22-33 1-42
    nonreactive_zone 42-42 42-42 fix
    END
