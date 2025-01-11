# Boundary Conditions Block

Boundary conditions are required for any problem involving transport.
The input currently involves specifying boundary conditions, corresponding 
essentially to "ghost cells" outside of the domain. Two possible approaches are possible:  1) grid cell by grid cell specific zones along boundary (recommended), or 2) the legacy approach in which a single boundary condition is applied across the entire face. 

The preferred way to specify boundary conditions is
```
     boundarycondition   condition_name  zone [jxbegin-jxend   jybegin-jyend jzbegin-jzend]  boundary_condition_type 
     boundarycondition   condition_name  zone [jxbegin-jxend   jybegin-jyend   jzbegin-jzend]  boundary_condition_type 
```
<u> Default:</u> &nbsp; None. All grid cells along a face must be specified.

Alternatively, the legacy approach can be used, but this is deprecated since it does not allow for heterogeneity along the boundary.

    x_begin  condition_name  boundary_condition_type
    x_end  condition_name  boundary_condition_type
    [y_begin  condition_name  boundary_condition_type]
    [y_end  condition_name  boundary_condition_type]
    [z_begin  condition_name  boundary_condition_type]
    [z_end  condition_name  boundary_condition_type]

<u> Default:</u> &nbsp; None.

You cannot mix the these two ways to specify boundary conditions.

## Boundary Condition Types

For fixed flux problems, it may be easier to make an injection well with the correct flux at the grid cell along the boundary. The volumetric source term in the input file in units of l/s is converted to m$^3$/yr and then to a Darcy flux by dividing by the cross-sectional area of the face, thus yielding units of m$^3$/m$^2$/yr. If this approach is used, then the actual boundary should be made to be no flux so that the fluid entering the cell (e.g., grid cell JX = 1) comes only from the volumetric source term.

<u>*Flux:*</u> &nbsp; In the case of a flux specification, the code uses the value and direction of the advective flux that is provided. If the flux is into the domain, then the upstream ghost cell will be used, if out of the domain (the downstream node), then there is no effect. For a pure diffusion problem (i.e., with no advective flux via flow across the boundary of the domain), a flux specification will result in a no-flux boundary.

<u>*Dirichlet:*</u> &nbsp; In the case of a Dirichlet boundary in GIMRT mode, this will cause both advective and diffusive flux to be calculated.
To use a Dirichlet condition in OS3D mode requires that an initial condition be fixed just inside the boundary of the problem.

## Example
For a 1D domain with 100 grid cells (NX=100), the specification for the preferred
option would be:
```
    boundarycondition   ambient  zone 0-0 1-1 1-1  flux 
    boundarycondition   injection  zone 101-101 1-1  1-1  dirichlet 
```
In this case, the ghost cells are jx=0 and jx=101.