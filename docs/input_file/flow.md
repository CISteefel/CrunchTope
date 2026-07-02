# Flow Block

[Basic keywords](#basic-keywords)
- [calculate_flow](#calculate_flow)
- [constant_flow](#constant_flow)
- [gravity](#gravity)
- [infiltration](#infiltration)
- [permeability_x](#permeability_x)
- [permeability_y](#permeability_y)
- [permeability_z](#permeability_z)
- [read_PermeabilityFile](#read_PermeabilityFile)


[Keywords for saturated flow](#keywords-for-saturated-flow)
- [initialize_hydrostatic](#initialize_hydrostatic)
- [pressure](#pressure)

[Keywords for unsaturated flow](#keywords-for-unsaturated-flow)
- [Richards](#richards)
- [Richards_print](#richards_print)
- [hydraulic_function](#hydraulic_function)
- [richards_ic](#richards_ic)
- [read_richards_ic_file](#read_richards_ic_file)
- [psi_is_head](#psi_is_head)
- [theta_r_is_S_r](#theta_r_is_s_r)
- [theta_s_is_porosity](#theta_s_is_porosity)
- [vg_alpha](#vg_alpha)
- [vg_n](#vg_n)
- [vg_theta_r](#vg_theta_r)
- [vg_theta_s](#vg_theta_s)
- [read_vg_alpha](#read_vg_alpha)
- [read_vg_n](#read_vg_n)
- [read_vg_theta_r](#read_vg_theta_r)
- [boundary_condition](#boundary_condition)
- [boundary_condition_steady](#boundary_condition_steady)
- [x_begin_bc_type](#x_begin_bc_type)
- [x_end_bc_type ](#x_end_bc_type)
- [x_begin_bc_type_steady](#x_begin_bc_type_steady)
- [x_end_bc_type_steady](#x_end_bc_type_steady)
- [set_evaporation_boundary](#set_evaporation_boundary)
- [set_psi_0](#set_psi_0)
- [spatial_domain](#spatial_domain)
- [set_dpsi_max](#set_dpsi_max)
- [set_tol_a](#set_tol_a)
- [set_max_Newton](#set_max_Newton)
- [set_max_line_search](#set_max_line_search)

[Keywords for gas transport](#keywords-for-gas-transport)
- [constant_gasflow](#constant_gasflow)
- [gaspump](#gaspump)
- [read_GasVelocityFile](#read_gasvelocityfile)

## Basic keywords

### calculate_flow

#### Syntax
```
calculate_flow [logical]
```
[logical] = true or yes or false or no (Default=false)

#### Explanation
This keyword indicates whether flow will be calculated according to Darcy’s Law or not.
Once this logical is set to true, the code will look for inputs on the pressure and permeability fields.

#### Example

```
calculate_flow true
```

### constant_flow

#### Syntax
```
constant_flow [qx] [qy] [qz]
```
[qx] [qy] [qz] are real numbers giving the X, Y, and Z components of the Darcy flux (Default: 0.0 m3/m2/yr). 

#### Explanation
This keyword is used to set a constant Darcy flux.  As a vector, it includes X, Y, and Z components.
In a one-dimensional problem in X, the Y and Z velocities are optional.
Units are changed inside the FLOW block with the space_units and time_units keywords.

#### Example

```
constant_flow 0.5 0.2 0.0
```

### gravity

#### Syntax
```
gravity [Angle to X] [Angle to Y] [Angle to Z] [direction]
```
[Angle to X], [Angle to Y] and [Angle to Z] are real numbers  (Default: 90).  
[direction] takes the value down or up.

#### Explanation
This keyword provides the angle relative to the X, Y, and Z coordinate directions for the gravity vector.
The influence of gravity is then calculated from

$$g = \cos\theta$$

where $\theta$ is the angle of the gravity vector relative to the specific coordinate direction (X, Y, or Z).
Direction can be either down or up and is used to indicate whether the direction the gravity vector is assumed to operate (for example, toward X[NX] or toward X[0]).
So, for example, the following input

```
gravity    90.0  0.0  90.0  down
```

would be used to set the gravity vector in the Y direction, with gravity oriented down (toward JY = 0).

#### Example

```
gravity    90.0  0.0  90.0  down
```


### infiltration

#### Syntax
```
infiltration [rate] [label] X/Y/Z
```
[rate] is a real number giving pumping rate (Default: none).
[label] gives the name of the geochemical condition that applies to the boundary.
X/Y/Z is the coordinate direction to be used.

#### Explanation
This keyword specifies an infiltration or recharge rate which is distributed over the boundary of the system (at JX, JY, or JZ = 0).
The geochemical condition to be used is given by its label.
#### Example

```
infiltration 0.5 INF_TEST Y
```



### infiltration

#### Syntax
```
pump [rate] [label] [JX] [JY] [JZ]
```
[rate] is a real number giving pumping rate (Default: none).  
[label] gives the name of the geochemical condition that applies to the well.
[JX] [JY] [JZ] are three integers setting the coordinates of the well.

#### Explanation
This keyword beginning with pump sets a pumping rate (positive for injection, negative for withdrawal) in units of L/s and assigns a geochemical condition to it.  This is followed by the coordinates of the well, with JY and JZ being optional in the case of a one-dimensional problem.
The geochemical condition is identified by its label, although the condition is only used in the case of an injection well.
#### Example

```
pump -2 WELL2 2 4 1
```



### permeability_x

#### Syntax
```
permeability_x [value] zone  [jxbegin-jxend   jybegin-jyend   jzbegin-jzend]
```
or
```
permeability_x [value] default
```
[value] is a positive real number  (Default: 1.0). 
[jxbegin-jxend   jybegin-jyend   jzbegin-jzend] are X, Y, and Z coordinates in the form of a beginning JX, and ending JX, a beginning JY and an ending JY, and a beginning JZ and an ending JZ, in each case separated by a hyphen (-).

#### Explanation
The permeability_x keyword uses an input format similar to other parameters specified by zone (e.g., pressure, tortuosity).
The keyword is followed by the value of the X permeability, which is then followed by either the label zone or default.
If zone is specified, then the code expects the X, Y, and Z coordinates in the form of a beginning JX, and ending JX, a beginning JY and an ending JY, and a beginning JZ and an ending JZ, in each case separated by a hyphen.
If default is specified instead of zone, the code will initialize the entire spatial domain (i.e., all grid cells) to the value.
More than one specification of permeability_x can be (and normally is) provided.
For example, in a 2D system with 20 grid cells in the X coordinate direction and 20 in the Y coordinate direction, a higher X permeability zone running down the center of the domain could be specified with
```
permeability_x    1.E-13  default
permeability_x    1.E-12  zone  0-21 10-10 1-1
```
Note that values of the X permeability at grid cell 0 and NX+1 have been specified, since the harmonic mean of the permeability is calculated at the system boundaries using these “ghost cell” values.
Zones specified later in the sequence overwrite zones specified earlier—in the example above, the default specification sets the permeability_x value to 1.E-13 m2 initially in all of the grid cells, while the next permeability_x keyword creates a higher permeability zone in grid cells JY=10 over the entire X length of the domain.
To create a no-flow boundary at JX=0 over the entire Y length, one could use:
```
permeability_x    1.E-13  default
permeability_x    1.E-12  zone  0-0 1-20 1-1
```
or more selectively with something like:
```
permeability_x    1.E-13  default
permeability_x    1.E-12  zone  0-0 1-8 1-1
permeability_x    1.E-12  zone  0-0 12-20 1-1
```
which would then create a no-flow boundary at JX=0 except for JY nodes 9, 10, and 11.
The convention adopted here is that all three dimensions must be specified as zones, even in a one-dimensional problem, although permeability_y and permeability_z are not needed.

#### Example

```
permeability_x    1.E-12  zone  0-0 12-20 1-1
```


### permeability_y

#### Syntax
```
permeability_y [value] zone  [jxbegin-jxend   jybegin-jyend   jzbegin-jzend]
```
or
```
permeability_y [value] default
```
[value] is a positive real number  (Default: 1.0). 
[jxbegin-jxend   jybegin-jyend   jzbegin-jzend] are X, Y, and Z coordinates in the form of a beginning JX, and ending JX, a beginning JY and an ending JY, and a beginning JZ and an ending JZ, in each case separated by a hyphen (-).

#### Explanation
The permeability_y keyword uses an input format similar to other parameters specified by zone (e.g., pressure, tortuosity).
The keyword is followed by the value of the X permeability, which is then followed by either the label zone or default.
If zone is specified, then the code expects the X, Y, and Z coordinates in the form of a beginning JX, and ending JX, a beginning JY and an ending JY, and a beginning JZ and an ending JZ, in each case separated by a hyphen.
If default is specified instead of zone, the code will initialize the entire spatial domain (i.e., all grid cells) to the value.
More than one specification of permeability_y can be (and normally is) provided.  For example, in a 2D system with 20 grid cells in the X coordinate direction and 20 in the Y coordinate direction, a higher X permeability zone running down the center of the domain could be specified with
```
permeability_y    1.E-13  default
permeability_y    1.E-12  zone  0-21 10-10 1-1
```
Note that values of the X permeability at grid cell 0 and NX+1 have been specified, since the harmonic mean of the permeability is calculated at the system boundaries using these “ghost cell” values.
Zones specified later in the sequence overwrite zones specified earlier—in the example above, the default specification sets the permeability_y value to 1.E-13 m2 initially in all of the grid cells, while the next permeability_y keyword creates a higher permeability zone in grid cells JY=10 over the entire X length of the domain.
To create a no-flow boundary at JX=0 over the entire Y length, one could use:
```
permeability_y    1.E-13  default
permeability_y    1.E-12  zone  0-0 1-20 1-1
```
or more selectively with something like:
```
permeability_y    1.E-13  default
permeability_y    1.E-12  zone  0-0 1-8 1-1
permeability_y    1.E-12  zone  0-0 12-20 1-1
```
which would then create a no-flow boundary at JX=0 except for JY nodes 9, 10, and 11.
The convention adopted here is that all three dimensions must be specified as zones, even in a one-dimensional problem, although permeability_y and permeability_z are not needed.

#### Example

```
permeability_y    1.E-12  zone  0-0 12-20 1-1
```

### permeability_z

#### Syntax
```
permeability_z [value] zone  [jxbegin-jxend   jybegin-jyend   jzbegin-jzend]
```
or
```
permeability_z [value] default
```
[value] is a positive real number  (Default: 1.0). 
[jxbegin-jxend   jybegin-jyend   jzbegin-jzend] are X, Y, and Z coordinates in the form of a beginning JX, and ending JX, a beginning JY and an ending JY, and a beginning JZ and an ending JZ, in each case separated by a hyphen (-).

#### Explanation
The permeability_z keyword uses an input format similar to other parameters specified by zone (e.g., pressure, tortuosity).
The keyword is followed by the value of the X permeability, which is then followed by either the label zone or default.
If zone is specified, then the code expects the X, Y, and Z coordinates in the form of a beginning JX, and ending JX, a beginning JY and an ending JY, and a beginning JZ and an ending JZ, in each case separated by a hyphen.
If default is specified instead of zone, the code will initialize the entire spatial domain (i.e., all grid cells) to the value.  More than one specification of permeability_z can be (and normally is) provided.
For example, in a 2D system with 20 grid cells in the X coordinate direction and 20 in the Y coordinate direction, a higher X permeability zone running down the center of the domain could be specified with
```
permeability_z    1.E-13  default
permeability_z    1.E-12  zone  0-21 10-10 1-1
```
Note that values of the X permeability at grid cell 0 and NX+1 have been specified, since the harmonic mean of the permeability is calculated at the system boundaries using these “ghost cell” values.
Zones specified later in the sequence overwrite zones specified earlier—in the example above, the default specification sets the permeability_z value to 1.E-13 m2 initially in all of the grid cells, while the next permeability_z keyword creates a higher permeability zone in grid cells JY=10 over the entire X length of the domain.
To create a no-flow boundary at JX=0 over the entire Y length, one could use:
```
permeability_z    1.E-13  default
permeability_z    1.E-12  zone  0-0 1-20 1-1
```
or more selectively with something like:
```
permeability_z    1.E-13  default
permeability_z    1.E-12  zone  0-0 1-8 1-1
permeability_z    1.E-12  zone  0-0 12-20 1-1
```
which would then create a no-flow boundary at JX=0 except for JY nodes 9, 10, and 11.
The convention adopted here is that all three dimensions must be specified as zones, even in a one-dimensional problem, although permeability_y and permeability_z are not needed.

#### Example

```
permeability_z    1.E-12  zone  0-0 12-20 1-1
```


### read_PermeabilityFile

#### Syntax
```
read_PermeabilityFile [filename] [format]
```
[filename] gives the name of the file (up to 132 characters) containing permeability values over the entire spatial domain (Default: none).
The extensions vx, vy, and vz are assumed to indicate the files containing the X, Y, and Z permeabilities, respectively.
[format]= SingleColumn, ContinuousRead, FullFormat, or Unformatted (Default: SingleColumn) 

#### Explanation
This keyword provides a filename prefix for ASCII files containing permeability values over the entire spatial domain.
The extensions vx, vy, and vz are assumed to indicate the files containing the X, Y, and Z permeabilities, respectively.
For example, the syntax
```
read_PermeabilityFile  copperleach   SingleColumn
```
would open and attempt to read the files:
```
copperleach.vx
copperleach.vy
copperleach.vz
```
using the SingleColumn format.
The format of the file depends on the optional format specified (see Format of Additionnal Input Files).
Units of the permeability are m2.
The format of the file depends on the optional format specified (see Format of Additionnal Input Files). 
The format of the file depends on the file format specified.  In none is provided, a single column format consisting of a single permeability per line, with NX varying first, then NY, and then NZ is assumed.
Unlike the fluid velocity, permeabilities are defined at the center of the nodes, since at least in some cases they are linked to properties like the porosity, which in turn may be affected by mineral dissolution and precipitation.
Interface permeabilities are calculated as harmonic means between adjacent cells.
For example, between nodes i and i+1, the harmonic mean of the permeabililty is given by:
$$k_{i+1/2} = \frac{2k_i k_{i+1}}{k_i + k_{i+1}}$$

Since permeabilities need to be defined at the boundaries to determine the amount of flow across them, “ghost cell” X permeabilities are required at the 0 and nx+1 nodes, and so on.
For 1D problems, permy and permz files are not needed. 
For those familiar with FORTRAN, the actual source code in the case where a SingleColumn format is specified is given by (compare to the ContinuousRead format in the case of the read_VelocityFile keyword):
```
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 0,nx+1
      READ(52,*) permx(jx,jy,jz)
    END DO
  END DO
END DO
DO jz = 1,nz
  DO jy = 0,ny+1
    DO jx = 1,nx
      READ(53,*) permy(jx,jy,jz)
    END DO
  END DO
END DO
DO jz = 0,nz+1
  DO jy = 1,ny
    DO jx = 1,nx
      READ(54,*) permz(jx,jy,jz)
    END DO
  END DO
END DO

```

#### Example

```
read_PermeabilityFile  copperleach   SingleColumn
```

### read_VelocityFile

#### Syntax
```
read_VelocityFile [filename] [format]
```
[filename] gives the name of the file (up to 132 characters) containing velocity vector values over the entire spatial domain (Default: none).
The extensions vx, vy, and vz are assumed to indicate the files containing the X, Y, and Z velocities, respectively
[format]= SingleColumn, ContinuousRead, FullFormat, or Unformatted (Default: SingleColumn) 

#### Explanation
This keyword provides a filename prefix for ASCII files containing velocity vector values over the entire spatial domain.
The extensions vx, vy, and vz are assumed to indicate the files containing the X, Y, and Z velocities, respectively.
For example, the syntax
```
read_VelocityFile  copperleach   SingleColumn
```
would open and attempt to read the files:
```
copperleach.vx
copperleach.vy
copperleach.vz
```
using the SingleColumn format.
The format of the file depends on the optional format specified (see Format of Additionnal Input Files).
Units are set with the space_units and time_units keywords—otherwise, default units are m3/m2/yr.
The format of the file depends on the optional format specified (see Format of Additionnal Input Files).
In the absence of a specification of the file format, a single column (SingleColumn) format is assumed.
Velocities are defined at grid interfaces (not at the center of the nodes), so NX+1 X velocities are needed in the X coordinate direction, NY+1 Y velocities in the Y coordinate direction, and NZ+1 Z velocities in the Z coordinate direction.
CrunchFlow uses the convention that velocities for a particular grid cell are defined at the right hand interface (the so-called “right hand rule”).
So, for example, the velocity, qx(2,1,1) is assumed to be the X velocity at the interface between grid cell 2 and 3 (Figure XX).
Accordingly, the first entry in the velocityfile.vx file will be the X velocity defined at the interface between grid cell (0,1,1) and (1,1,1).
So, in the case of where ContinuousRead is specified as the file format, the code would read:
```
READ(52,*) (qx(jx,jy,jz),jx=0,nx),jy=1,ny),jz=1,nz)
READ(53,*) (qy(jx,jy,jz),jx=1,nx),jy=0,ny),jz=1,nz)
READ(54,*) (qz(jx,jy,jz),jx=1,nx),jy=1,ny),jz=0,nz)
```

#### Example

```
read_VelocityFile  copperleach   SingleColumn
```




## Keywords for saturated flow

### initialize_hydrostatic

#### Syntax
```
initialize_hydrostatic [logical]
```
[logical] = true or yes or false or no (Default=false) 

#### Explanation
The keyword parameter initialize_hydrostatic is followed by true or false (or yes or no) and is used to select, when true, the option to initialize the pressure field to hydrostatic conditions before time stepping begins.
When set to true, hydrostatic conditions will be initialized no matter what other pressure boundary conditions or source/sink terms take effect once time stepping begins.
However, a coordinate direction for the gravity vector must be selected with the gravity keyword and appropriate boundary conditions (typically just a fixed pressure at one end of the domain).
#### Example

```
initialize_hydrostatic true
```

### pressure

#### Syntax
```
pressure [value] zone  [jxbegin-jxend   jybegin-jyend   jzbegin-jzend] [fix]
```
or
```
pressure [value] default
```
[value] is a positive real number  (Default: none). 
[jxbegin-jxend   jybegin-jyend   jzbegin-jzend] are X, Y, and Z coordinates in the form of a beginning JX, and ending JX, a beginning JY and an ending JY, and a beginning JZ and an ending JZ, in each case separated by a hyphen (-).

#### Explanation
The pressure keyword uses an input format similar to other parameters specified by zone (e.g., tortuosity, permeability).
The keyword is followed by the value of the fluid pressure in units of Pascals, which is then followed by either the label zone or default. 
If zone is specified, then the code expects the X, Y, and Z coordinates in the form of a beginning JX, and ending JX, a beginning JY and an ending JY, and a beginning JZ and an ending JZ, in each case separated by a hyphen (-).
If default is specified instead of zone, the code will initialize the entire spatial domain (i.e., all grid cells) to the value.
More than one specification of pressure can be (and normally is) provided, for example
```
pressure    1000.0  default
pressure    0.0     zone  1-20 21-21 1-1 fix
```
Zones specified later in the sequence overwrite zones specified earlier—in the example above, the default specification sets the fluid pressure to 1000 Pa initially in all of the grid cells, while the next keyword sets the pressure at JY = 21 to zero.
In the case of a 20 by 20 X-Y simulation, the grid cell 21 corresponds to node NY+1 and is therefore a ghost cell.
This is the normal way in which fixed pressure boundary conditions are set.
Following the zone coordinates with fix instructs the code to lock the value in for the course of the simulation—not including fix would mean that the pressure is only set as an initial condition, after which it is free to evolve according to the physics of the problem.
The convention adopted here is that all three dimensions must be specified, even in a one-dimensional problem.
To create a hydrostatic pressure gradient in the Y direction, with depth increasing with increasing JY, one could use:
```
pressure    1000.0  default
pressure    0.0     zone  1-20 0-0 1-1 fix
gravity     90.0  0.0  90.0  up
```
With this input, the code would time step until the hydrostatic pressure gradient was achieved in the absence of other forcings.  To initialize the system to hydrostatic pressure before time stepping begins, use:
```
pressure    1000.0  default
pressure    0.0     zone  1-20 0-0 1-1 fix
gravity     90.0  0.0  90.0  up
```


#### Example

```
pressure    0.0     zone  1-20 0-0 1-1 fix
```


## Keywords for unsaturated flow

### Richards

#### Syntax
```
Richards [logical]
```
[logical] = true or false (Default=false) 

#### Explanation
This keyword indicates whether flow will be calculated according to the one-dimensional Richards equation.
Currently, the solver only supports the one-dimensional Richards equation.
Once the logical is set to true, the code will look for inputs on the soil hydraulic parameters and the initial and boundary conditions used to solve the one-dimensional Richards equation.

#### Example

```
Richards true
```


### Richards_steady

#### Syntax
```
Richards_steady [logical]
```
[logical] = true or false (Default=false) 

#### Explanation
This keyword indicates whether the initial condition of unsaturated flow will be calculated by solving the steady-state Richards equation.
Once the logical is set to true, the code will look for inputs on the boundary conditions used to solve the one-dimensional steady-state Richards equation.
If the user provides the initial condition, the solver uses it as the initial guess of the solution.
Otherwise, the initial guess is the zero water potential in the whole spatial domain.

#### Example

```
Richards_steady true
```

### Richards_print

#### Syntax
```
Richards_print [logical]
```
[logical] = true or false (Default=false) 

#### Explanation
This keyword indicates whether the screen shows print statements from the Richards solver (e.g., number of Newton iterations and line searches etc.).
This should be turned on only when you want to diagnose the convergence.

#### Example

```
Richards_print true
```

### richards_ic

#### Syntax
```
richards_ic [value]
```
[value] is the initial water potential used for both steady-state and time-dependent Richards solver.

#### Explanation
This keyword sets the initial condition for the water potential for the Richards solver to a constant value throughout the spatial domain.
This keyword is read only when the keyword 'read_richards_ic_file' is not provided.

#### Example
``` 
richards_ic -1.0
```


### read_richards_ic_file

#### Syntax
```
read_richards_ic_file [filename] [format]
```
[filename] gives the name of the file (up to 132 characters) containing initial condition (water potential or head) values for the Richards equation over the entire spatial domain.
[format]= SingleColumn (Default: SingleColumn)

#### Explanation
Read the initial condition for the Richards equation from a file.
If this keyword is used, the keyword 'richards_ic' is ignored.

#### Example
```
read_richards_ic_file initial_test_infiltration.csv SingleColumn
```



### psi_is_head

#### Syntax
```
psi_is_head [logical]
```
[logical] = true or false (Default=true) 

#### Explanation
This keyword indicates whether the primary variable psi in the Richards equation used in the input files is pressure head or not.
If false is selected, the input values (initial and boundary conditions, and vg_alpha) are interpreted as pressure [Pa].
In that case, the input pressure value is converted to pressure head by 
$$\text{pressure head} = \frac{\text{pressure [Pa]} - \text{air pressure [Pa]}}{998.23 \text{[kg m-3]} * 9.80665 \text{[m s-2]}}$$
, where pressure_air = 101325 Pa.
This keyward was added to be compatible with PFLOTRAN input files.
Also, when false is selected, the unit for space must be meters in the FLOW block.

#### Example

```
psi_is_head true
```




### theta_r_is_S_r

#### Syntax
```
theta_r_is_S_r [logical]
```
[logical] = true or false (Default=false) 

#### Explanation
This keyword indicates whether the input to the $\theta_r$ parameter (“residual water content”) in the van Genuchten model is residual saturation, or not.
If true is selected, the input value (keyword vg_theta_r) is interpreted as the residual saturation and converted into the residual water content by residual_water_content = $\theta_r \theta_s$.

#### Example

```
theta_r_is_S_r true
```

### theta_s_is_porosity

#### Syntax
```
theta_s_is_porosity [logical]
```
[logical] = true or false (Default=true) 

#### Explanation
This keyword indicates whether the $\theta_s$ parameter in the van Genuchten model is the same as the porosity value, or not.
If false is selected, users need to provide information on the $\theta_s$ parameter by the keyword vg_theta_s.
Otherwise, porosity value is used for the $theta_s$ parameter in the van Genuchten model, so the user does not need to provide $\theta_s$ value.

#### Example

```
theta_s_is_porosity true
```


### vg_is_n

#### Syntax
```
vg_is_n [logical]
```
[logical] = true or false (Default=true) 

#### Explanation
This keyword indicates whether the van Genuchten parameter $n$ is used in the input (vg_n), or not.
If false is selected, the input (vg_n) is interpreted as the $m$ parameter in the van Genuchten model.
In that case, the input value is converted to n by the formula: $n = 1/(1-m)$.
This keyward was added to be compatible with PFLOTRAN input file.


#### Example
```
vg_is_n true
```


### vg_alpha

#### Syntax
```
vg_alpha [#cells spacing] [#cells spacing] & 
                     ...[#cells spacing]
```
#cells is an integer and spacing is a real number

#### Explanation
This keyword parameter specifies the $\alpha$ parameter used in the van Genuchten model over the spatial domain.
This parameter is given by specifying pairs of integers and real numbers which give the number of cells of a given parameter value.
The following example would result in 10 cells of $\alpha$ = 1.00, 10 cells of $\alpha$ = 1.5, and 20 cells of $\alpha$ = 2.00.
The unit of the parameter is in length$^{-1}$ and follows the unit of the Flow block.
If this keyward and the values are not provide, you will get an error unless you provide the keyward read_vg_alpha.


#### Example
```
vg_alpha  10 1.00 10 1.5 20 2.0
```


### vg_n

#### Syntax
```
vg_n [#cells spacing] [#cells spacing] & 
                     ...[#cells spacing]
```
#cells is an integer and spacing is a real number

#### Explanation
This keyword parameter specifies the $n$ parameter used in the van Genuchten model over the spatial domain.
This parameter is given by specifying pairs of integers and real numbers which give the number of cells of a given parameter value.
The following example would result in 10 cells of $n$ = 1.00, 10 cells of $n$ = 1.5, and 20 cells of $n$ = 2.00.
If this keyward and the values are not provide, you will get an error unless you provide the keyward read_vg_n.


#### Example
```
vg_n  10 1.00 10 1.5 20 2.0
```

### vg_theta_r

#### Syntax
```
vg_theta_r [#cells spacing] [#cells spacing] & 
                     ...[#cells spacing] 
```
#cells is an integer and spacing is a real number

#### Explanation
This keyword parameter specifies the residual water content $\theta_r$ used in the van Genuchten model over the spatial domain.
This parameter is given by specifying pairs of integers and real numbers which give the number of cells of a given parameter value.
The following example would result in 10 cells of $\theta_r$ = 0.05, 10 cells of $\theta_r$ = 0.08, and 20 cells of $\theta_r$ = 0.1.
If this keyward and the values are not provide, you will get an error unless you provide the keyward read_vg_theta_r.


#### Example
```
vg_theta_r  10 0.05 10 0.08 20 0.1
```

### vg_theta_s

#### Syntax
```
vg_theta_s [#cells spacing] [#cells spacing] & 
                     ...[#cells spacing] 
```
#cells is an integer and spacing is a real number

#### Explanation
This keyword parameter specifies the residual water content $\theta_s$ used in the van Genuchten model over the spatial domain.
This parameter is given by specifying pairs of integers and real numbers which give the number of cells of a given parameter value.
The following example would result in 10 cells of theta_s = 0.40, 10 cells of theta_s = 0.30, and 20 cells of theta_s = 0.35.
If this keyward and the values are not provide, you will get an error.
The saturated water content is supposed to be equal to the porosity.
If these numbers are not equal, you will get a warning.


#### Example
```
vg_theta_s  10 0.40 10 0.30 20 0.35
```


### read_vg_alpha

#### Syntax
```
read_vg_alpha
```

#### Explanation
Read alpha parameter in the van Genuchten model from a file.

#### Example


### read_vg_n

#### Syntax
```
read_vg_n
```
#### Explanation
Read n parameter in the van Genuchten model from a file.

#### Example


### read_vg_theta_r 

#### Syntax
```
read_vg_theta_r 
```
#### Explanation
Read theta_r parameter in the van Genuchten model from a file.

#### Example


### x_begin_bc_type 

#### Syntax
```
x_begin_bc_type [BC_type] [value or file name] [number of time steps in the file]
```
#### Explanation
This keyword set the type and the type of boundary condition applied to the x_begin boundary (i.e., between jx = 0 and jx = 1).
Currently, the following boundary conditions are available:
- constant_Dirichlet
- constant_neumann
- constant_flux
- constant_atomosphere
- variable_Dirichlet
- variable_neumann
- variable_flux
- variable_atomosphere

Dirichlet boundary conditions enforce the water potential at the boundary to be the specified value.
Neumann boundary conditions enforce the gradient of the water potential at the boundary to be the specified value.
Flux boundary conditions enforce the water flux at the boundary to be the specified value.
Atomosphere boundary conditions swtich between Dirichlet and flux boundary conditions depending on the water potentail at the surface.
Until the boundary water potential is above threshold value $\psi_0$ ($-10^4$ m by default but can be chanaged by the keyword set_psi_0), flux boundary condition is used.
Once the boundary water potential is below the threshold (i.e., the surface soil is extremely dry), Dirichlet boundary condition with $psi_0$ is applied. 
Thus, the value for the flux boundary condition needs to be provided here.
When selecting constant boundary conditions, the value for the constnat boundary condition needs to be provided.
When selecting variable boundary conditions, the file name for the transient boundary condition data needs to be provided.
The format of the file is explained in the example below.

#### Example
```
x_begin_bc_type constant_flux 0.2
```
for constant flux boundary condition.

```
x_begin_bc_type variable_flux water_upper_BC.dat 12
```
for transient flux boundary condition. Here, the water_upper_BC.dat file is

```
0.000000000000000000e+00 1.000000000000000056e-01
2.500000000000000000e-01 1.000000000000000056e-01
2.510000000000000009e-01 0.000000000000000000e+00
3.000000000000000000e+00 0.000000000000000000e+00
3.009999999999999787e+00 1.000000000000000056e-01
3.250000000000000000e+00 1.000000000000000056e-01
3.250999999999999890e+00 0.000000000000000000e+00
6.000000000000000000e+00 0.000000000000000000e+00
6.009999999999999787e+00 1.000000000000000056e-01
6.250000000000000000e+00 1.000000000000000056e-01
6.251000000000000334e+00 0.000000000000000000e+00
9.000000000000000000e+00 0.000000000000000000e+00
```
, where the first column is the time, and the second column is the value. 
These data are interpolated for during time stepping to get a value used for each step.


### x_end_bc_type 

#### Syntax
```
x_end_bc_type [BC_type] [value or file name] [number of time steps in the file]
```
#### Explanation
This keyword set the type and the type of boundary condition applied to the x_end boundary (i.e., between jx = nx and jx = nx+1).
Currently, the following boundary conditions are available:
- constant_Dirichlet
- constant_neumann
- constant_flux
- variable_Dirichlet
- variable_neumann
- variable_flux

Dirichlet boundary conditions enforce the water potential at the boundary to be the specified value.
Neumann boundary conditions enforce the gradient of the water potential at the boundary to be the specified value.
Flux boundary conditions enforce the water flux at the boundary to be the specified value.
When selecting constant boundary conditions, the value for the constnat boundary condition needs to be provided.
When selecting variable boundary conditions, the file name for the transient boundary condition data needs to be provided.
The format of the file is explained in the example below.

#### Example
```
x_end_bc_type constant_neumann 0.0
```
for constant Neumann boundary condition.

### x_begin_bc_type_steady

#### Syntax
```
x_begin_bc_type_steady [BC_type] [value]
```
#### Explanation
This keyword set the type and the type of boundary condition applied to the x_begin boundary (i.e., between jx = 0 and jx = 1) for steady-state Richards equation.
Currently, the following boundary conditions are available:
- constant_Dirichlet
- constant_neumann
- constant_flux

Dirichlet boundary conditions enforce the water potential at the boundary to be the specified value.
Neumann boundary conditions enforce the gradient of the water potential at the boundary to be the specified value.
Flux boundary conditions enforce the water flux at the boundary to be the specified value.
When selecting constant boundary conditions, the value for the constnat boundary condition needs to be provided.

#### Example
```
x_begin_bc_type_steady constant_flux 0.0
```
for constant flux boundary condition for the steady-state Richards equation to get the intiial condition.

### x_end_bc_type_steady

#### Syntax
```
x_end_bc_type_steady [BC_type] [value]
```
#### Explanation
This keyword set the type and the type of boundary condition applied to the x_end boundary (i.e., between jx =  and jx = nx+1) for steady-state Richards equation.
Currently, the following boundary conditions are available:
- constant_Dirichlet
- constant_neumann
- constant_flux

Dirichlet boundary conditions enforce the water potential at the boundary to be the specified value.
Neumann boundary conditions enforce the gradient of the water potential at the boundary to be the specified value.
Flux boundary conditions enforce the water flux at the boundary to be the specified value.
When selecting constant boundary conditions, the value for the constnat boundary condition needs to be provided.

#### Example
```
x_end_bc_type_steady constant_Dirichlet 0.0
```
for constant Dirichlet boundary condition for the steady-state Richards equation to get the intiial condition.


### set_evaporation_boundary 
#### Syntax
```
set_evaporation_boundary [BC_location]
```

#### Explanation
This keyword sets the boundary, where no chemcial transport due to advection is applied.
This keyword is used to simulate the accumulation of chemical species near the soil surface due to evaporation.

#### Example
```
set_evaporation_boundary x_begin
```


### set_psi_0
#### Syntax
```
set_psi_0 [value]
```

#### Explanation
This keyword sets minimum water potentail allowed at the atomospheric boundary.
This value corresponds to the water potetial equal to the vapor in the air.
The default value is $\psi_0 = -10^4$ m.
The unit follows the unit in the flow block.

#### Example
```
set_psi_0 -1.0d3
```


### set_dpsi_max 
#### Syntax
```
set_dpsi_max [value]
```

#### Explanation
This keyword sets the maximum water potentail change allowed during Newton iterations.
The default value is $\psi_0 = -10^3$ m.
The unit follows the unit in the flow block.
When the Richards solver is not converged, a smaller value of this value may help the convergence (not always though).

#### Example
```
set_dpsi_max 1.0d2
```


## Keywords for gas transport

### constant_gasflow

#### Syntax
```
constant_gasflow [qxgas] [qygas] [qzgas]
```
[qxgas] [qygas] [qzgas] are real numbers giving the X, Y, and Z components of the gas flux (Default: 0.0 m3/m2/yr). 

#### Explanation
This keyword is used to set a constant gas flux.  As a vector, it includes X, Y, and Z components.
In a one-dimensional problem in X, the Y and Z velocities are optional.
Units are changed inside the FLOW block with the space_units and time_units keywords.
#### Example

```
constant_gasflow 5 0.0 3.0
```


### gaspump

#### Syntax
```
gaspump [rate] [label] [JX] [JY] [JZ]
```
[rate] is a real number giving pumping rate (Default: none).
[label] gives the name of the geochemical condition that applies to the well.
[JX] [JY] [JZ] are three integers setting the coordinates of the well.

#### Explanation
This keyword beginning with gaspump sets a pumping rate (positive for injection, negative for withdrawal) in units of L/s and assigns a geochemical condition to it.
This is followed by the coordinates of the well, with JY and JZ being optional in the case of a one-dimensional problem.
The geochemical condition is identified by its label, although the condition is only used in the case of an injection well.

#### Example

```
gaspump -2 WELL2 2 4 1
```


### read_GasVelocityFile

#### Syntax
```
read_GasVelocityFile [filename] [format]
```
[filename] gives the name of the file (up to 132 characters) containing gas velocity vector values over the entire spatial domain (Default: none).
The extensions vx, vy, and vz are assumed to indicate the files containing the X, Y, and Z velocities, respectively.
[format]= SingleColumn, ContinuousRead, FullFormat, or Unformatted (Default: SingleColumn) 

#### Explanation
This keyword provides a filename prefix for ASCII files containing velocity vector values over the entire spatial domain.
The extensions vx, vy, and vz are assumed to indicate the files containing the X, Y, and Z velocities, respectively.
For example, the syntax
```
read_GasVelocityFile  co2stream  ContinuousRead
```
would open and attempt to read the files:
```
co2stream.vx
co2stream.vy
co2stream.vz
```
using the SingleColumn format.
The format of the file depends on the optional format specified (see Format of Additionnal Input Files).
Units are set with the space_units and time_units keywords—otherwise, default units are m3/m2/yr.
The format of the file depends on the optional format specified (see Format of Additionnal Input Files). 
In the absence of a specification of the file format, a single column (SingleColumn) format is assumed.
Velocities are defined at grid interfaces (not at the center of the nodes), so NX+1 X velocities are needed in the X coordinate direction, NY+1 Y velocities in the Y coordinate direction, and NZ+1 Z velocities in the Z coordinate direction.
CrunchFlow uses the convention that velocities for a particular grid cell are defined at the right hand interface (the so-called “right hand rule”).
So, for example, the velocity, qx(2,1,1) is assumed to be the X velocity at the interface between grid cell 2 and 3 (Figure XX).
Accordingly, the first entry in the gasvelocityfile.vx file will be the X velocity defined at the interface between grid cell (0,1,1) and (1,1,1).
So, in the case of where ContinuousRead is specified as the file format, the code would read:
```
READ(52,*) (qx(jx,jy,jz),jx=0,nx),jy=1,ny),jz=1,nz)
READ(53,*) (qy(jx,jy,jz),jx=1,nx),jy=0,ny),jz=1,nz)
READ(54,*) (qz(jx,jy,jz),jx=1,nx),jy=1,ny),jz=0,nz)
```

#### Example

```
read_GasVelocityFile  co2stream  ContinuousRead
```