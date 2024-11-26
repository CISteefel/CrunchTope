# Boundary Conditions Block


[Syntax](#syntax)

[Explanation](#explanation)

Boundary conditions are required for any problem involving transport.
The input currently involves specifying boundary conditions, corresponding essentially to “ghost cells” outside of the domain, at the 


## Syntax

```
x_begin  condition_name  boundary_condition_type
x_end  condition_name  boundary_condition_type
[y_begin  condition_name  boundary_condition_type]
[y_end  condition_name  boundary_condition_type]
[z_begin  condition_name  boundary_condition_type]
[z_end  condition_name  boundary_condition_type]
```
Default:  None.

Alternatively, you could specify boundary conditions via
```
boundarycondition   condition_name  zone [jxbegin-jxend   jybegin-jyend   jzbegin-jzend]  boundary_condition_type 
boundarycondition   condition_name  zone [jxbegin-jxend   jybegin-jyend   jzbegin-jzend]  boundary_condition_type 
```
You cannot mix the these two ways to specify boundary conditions.

## Explanation

Boundary conditions are set for ghost cells outside the faces of the 3D domain.
Currently, uniform boundary conditions are specified over the entire face, although this is being revised to provide more flexibility.
For additional flexibility, the user may fix initial conditions at the first or last grid cell in a coordinate direction, thus making it possible to specify multiple geochemical conditions on any one face.
For fixed flux problems, it may be easier to make an injection well with the correct flux at the grid cell along the boundary.
Boundary condition types are Dirichlet (used only in GIMRT mode) and flux.
In the case of a flux specification, the code uses the value and direction of the advective flux that is provided—if the flux is into the domain, then the upstream ghost cell will be used, if out of the domain, then there is no effect.
For a pure diffusion problem (i.e., with no advective flux via flow across the boundary of the domain), a flux specification will result in a no-flux boundary.
In the case of a Dirichlet boundary in GIMRT mode, this will cause both advective and diffusive flux to be calculated.
To use a Dirichlet condition in OS3D mode requires that an initial condition be fixed just inside the boundary of the problem.


## Example

```
boundarycondition   ambient  zone 0-0 1-1 1-1  flux 
boundarycondition   injection  zone 101-101 1-1  1-1  dirichlet 
```
