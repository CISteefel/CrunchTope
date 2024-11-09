MODULE Richards_module
  
USE crunchtype
USE hydraulic_function_module

IMPLICIT NONE
  
PRIVATE
  
  
TYPE :: RichardsBase ! store information used in the Richards solver, including discretization information
  CHARACTER (LEN=60) :: spatial_domain ! shape of the spatial domain (only used in 2D Richards)
  REAL(DP), ALLOCATABLE :: x(:)
  REAL(DP), ALLOCATABLE :: y(:)
  REAL(DP), ALLOCATABLE :: dx(:)
  REAL(DP), ALLOCATABLE :: dy(:)
  INTEGER(I4B) :: n_inner_cells
  INTEGER(I4B) :: n_faces
  INTEGER(I4B) :: n_bfaces
  INTEGER(I4B) :: n_total_cells
  INTEGER(I4B), ALLOCATABLE :: cell_to_face(:,:) ! link function from global cell number to global face number
  INTEGER(I4B), ALLOCATABLE :: face_to_cell(:,:) ! link function from global face number to global cell number
  INTEGER(I4B), ALLOCATABLE :: bface_to_face(:) ! link function from global boundary cell number to global face number
  INTEGER(I4B), ALLOCATABLE :: cell_to_coordinate(:,:) ! link function from cell number to spatial coordinate
  REAL(DP) :: psi_0 ! minimum water potential allowed when selecting 'atomosphere' boundary condition
END TYPE RichardsBase
  

TYPE :: RichardsOptions ! store options for Richards solver
  LOGICAL(LGT) :: is_steady = .FALSE. ! True when solving the steady-state Richards equation to get the initial condition
  LOGICAL(LGT) :: is_print = .FALSE. ! True if you want print statements from the Richards solver
  CHARACTER (LEN=60) :: hydraulic_function ! type of hydraulic function
  LOGICAL(LGT) :: vg_is_n = .TRUE. ! True if the input to vg_n is the n parameter in the van Genuchten model, otherwise, the input value is interpreted as the m parameter
  LOGICAL(LGT) :: psi_is_head = .TRUE. ! True if the primary variable psi in the Richards equation is pressure head [L] or not. If false, the input values for the initial and boundary conditions, and vg_alpha are interpreted as in terms of pressure [Pa].  
  LOGICAL(LGT) :: theta_s_is_porosity = .TRUE. ! True if the input to theta_s is the same as the porosity
  LOGICAL(LGT) :: theta_r_is_S_r = .FALSE. ! True if the input to theta_r is the residual saturation
END TYPE RichardsOptions
  
    
TYPE :: RichardsSolver ! store information on the nonlinear solver used in Richards solver
  REAL(DP) :: dpsi_max ! maximum update for water potential during the Newton iteration

END TYPE RichardsSolver

TYPE :: RichardsState ! store information on the nonlinear solver used in Richards solver
  REAL(DP), ALLOCATABLE :: psi(:,:,:) ! water potential [L]
  REAL(DP), ALLOCATABLE :: theta(:,:,:) ! volumetric water content [L3 L-3]
  REAL(DP), ALLOCATABLE :: theta_prev(:,:,:) ! volumetric water content from previous time step [L3 L-3]
  REAL(DP), ALLOCATABLE :: head(:,:,:) ! pressure head [L]
  REAL(DP), ALLOCATABLE :: kr(:,:,:) ! relative permeability [-]
  REAL(DP), ALLOCATABLE :: rho_water(:,:,:) ! the density of water [kg m-3]; this is temperature dependent
  REAL(DP), ALLOCATABLE :: mu_water(:,:,:) ! dynamics viscosity of water [Pa year]; this is temperature dependent
  REAL(DP), ALLOCATABLE :: xi(:,:,:) ! physical constant used in the Richards equation
  REAL(DP), ALLOCATABLE :: kr_faces(:) ! relative permeability at cell faces [-]
  REAL(DP), ALLOCATABLE :: xi_faces(:) ! physical constant in the Richards equation at faces    
END TYPE RichardsState

TYPE(RichardsBase), PUBLIC :: Richards_Base
TYPE(RichardsOptions), PUBLIC :: Richards_Options
TYPE(RichardsSolver), PUBLIC :: Richards_Solver
TYPE(RichardsState), PUBLIC :: Richards_State
  
PUBLIC RichardsAllocate, &
       RichardsDiscretize, &
       RichardsUpdateFluid
       
  
  CONTAINS

! ************************************************************************** !
SUBROUTINE AllocateArray(array, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
USE crunchtype
IMPLICIT NONE

REAL(DP), ALLOCATABLE, INTENT(INOUT) :: array(:,:,:)
INTEGER(I4B), INTENT(IN) :: lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z
  
! Check if the array is allocated, deallocate if necessary, then allocate with the new dimensions
IF (ALLOCATED(array)) THEN
  DEALLOCATE(array)
  ALLOCATE(array(lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z))
ELSE
  ALLOCATE(array(lbound_x:ubound_x, lbound_y:ubound_y, lbound_z:ubound_z))
END IF
   
END SUBROUTINE AllocateArray

! ************************************************************************** !
SUBROUTINE AllocateArray_1D(array, lbound_x, ubound_x)
USE crunchtype
IMPLICIT NONE

REAL(DP), ALLOCATABLE, INTENT(INOUT) :: array(:)
INTEGER(I4B), INTENT(IN) :: lbound_x, ubound_x
  
! Check if the array is allocated, deallocate if necessary, then allocate with the new dimensions
IF (ALLOCATED(array)) THEN
  DEALLOCATE(array)
  ALLOCATE(array(lbound_x:ubound_x))
ELSE
  ALLOCATE(array(lbound_x:ubound_x))
END IF
   
END SUBROUTINE AllocateArray_1D

! ************************************************************************** !
SUBROUTINE CreateLinkFunction(nx, ny, nz)
USE crunchtype

INTEGER(I4B), INTENT(IN) :: nx, ny, nz
INTEGER(I4B) :: i, j

IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
  IF (Richards_Base%spatial_domain == 'regular') THEN
    ! cell to face function
    ALLOCATE(Richards_Base%cell_to_face(Richards_Base%n_inner_cells, 4)) ! because the cell is limited to rectangular mesh, the column number is 4
        
    DO i = 1, Richards_Base%n_inner_cells
      Richards_Base%cell_to_face(i, 1) = i ! south face
      Richards_Base%cell_to_face(i, 2) = i + nx ! north face
      Richards_Base%cell_to_face(i, 3) = nx*(ny+1) + (nx+1)*((i-1)/nx) + MOD(i-1, nx) + 1 ! west face
      Richards_Base%cell_to_face(i, 4) = nx*(ny+1) + (nx+1)*((i-1)/nx) + MOD(i-1, nx) + 2 ! east face
    END DO 
        
    ! boundary face to face function
    ALLOCATE(Richards_Base%bface_to_face(Richards_Base%n_bfaces))
        
    !! bottom boundary
    DO i = 1, nx
      Richards_Base%bface_to_face(i) = i
    END DO
        
    !! right boundary
    DO i = 1, ny
      Richards_Base%bface_to_face(nx+i) = nx * (ny + 1) + (nx + 1) * i
    END DO
        
    !! top boundary
    DO i = 1, nx
      Richards_Base%bface_to_face(nx+ny+i) = nx*ny + (nx + 1 - i)
    END DO
        
    !! left boundary
    DO i = 1, ny
      Richards_Base%bface_to_face(2*nx+ny+i) = nx * (ny + 1) + (nx+1) * (ny - i) + 1
    END DO
         
    ! from cell number to coordinate numbers
    ALLOCATE(Richards_Base%cell_to_coordinate(Richards_Base%n_total_cells, 3))
    
    Richards_Base%cell_to_coordinate(:,3) = 1 ! z-coodinate is 1 for 2D problem
    
    DO i = 1, Richards_Base%n_inner_cells
      Richards_Base%cell_to_coordinate(i, 1) = MOD(i-1, nx)+1
      Richards_Base%cell_to_coordinate(i, 2) = (i-1)/nx + 1
    END DO
        
    ! loop over ghost cells
    !! bottom boundary
    DO i = 1, nx
      Richards_Base%cell_to_coordinate(Richards_Base%n_inner_cells + i, 1) = i
      Richards_Base%cell_to_coordinate(Richards_Base%n_inner_cells + i, 2) = 0
    END DO
        
    !! right boundary
    DO i = 1, ny
      Richards_Base%cell_to_coordinate(Richards_Base%n_inner_cells + nx + i, 1) = nx+1
      Richards_Base%cell_to_coordinate(Richards_Base%n_inner_cells + nx + i, 2) = i
    END DO
        
    !! top boundary
    DO i = 1, nx
      Richards_Base%cell_to_coordinate(Richards_Base%n_inner_cells + nx + ny + i, 1) = nx - i + 1
      Richards_Base%cell_to_coordinate(Richards_Base%n_inner_cells + nx + ny + i, 2) = ny + 1
    END DO
        
    !! left boundary
    DO i = 1, ny
      Richards_Base%cell_to_coordinate(Richards_Base%n_inner_cells + 2 * nx + ny + i, 1) = 0
      Richards_Base%cell_to_coordinate(Richards_Base%n_inner_cells + 2 * nx + ny + i, 2) = ny - i + 1
    END DO
        
        
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Currently, the spatial domain ', Richards_Base%spatial_domain, ' is supported.'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
    
  ! face to cell function
  ALLOCATE(Richards_Base%face_to_cell(Richards_Base%n_faces, 2))
        
  DO i = 1, Richards_Base%n_faces
    ! loop over internal cells
    DO j = 1, Richards_Base%n_inner_cells
      IF (Richards_Base%cell_to_face(j, 1) == i) THEN
        ! south face
        Richards_Base%face_to_cell(i, 2) = j
      ELSE IF (Richards_Base%cell_to_face(j, 2) == i) THEN
        ! north face
        Richards_Base%face_to_cell(i, 1) = j
      ELSE IF (Richards_Base%cell_to_face(j, 3) == i) THEN
        ! west face
        Richards_Base%face_to_cell(i, 2) = j
      ELSE IF (Richards_Base%cell_to_face(j, 4) == i) THEN
        ! east face
        Richards_Base%face_to_cell(i, 1) = j
      END IF
    END DO
    ! loop over ghost cells
    DO j = 1, Richards_Base%n_bfaces
      IF (Richards_Base%bface_to_face(j) == i) THEN
        IF (Richards_Base%face_to_cell(i, 1) == 0) THEN
          Richards_Base%face_to_cell(i, 1) = j + Richards_Base%n_inner_cells
        ELSE
          Richards_Base%face_to_cell(i, 2) = j + Richards_Base%n_inner_cells
        END IF
      END IF
    END DO
  END DO
  
ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Currently, three-dimensional Richards solver is supported.'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
      
END SUBROUTINE CreateLinkFunction
! ************************************************************************** !
SUBROUTINE RichardsAllocate(nx, ny, nz)
USE crunchtype
IMPLICIT NONE

INTEGER(I4B), INTENT(IN) :: nx, ny, nz
INTEGER(I4B) :: lbound_x = 1
INTEGER(I4B) :: ubound_x = 1
INTEGER(I4B) :: lbound_y = 1
INTEGER(I4B) :: ubound_y = 1
INTEGER(I4B) :: lbound_z = 1
INTEGER(I4B) :: ubound_z = 1

IF (nx > 1 .AND. ny == 1 .AND. nz == 1) THEN ! one-dimensional problem
  lbound_x = 0
  ubound_x = nx + 1

ELSE IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
  lbound_x = 0
  ubound_x = nx + 1
  lbound_y = 0
  ubound_y = ny + 1

END IF

! allocate state variables
CALL AllocateArray(Richards_State%psi, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(Richards_State%theta, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(Richards_State%theta_prev, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(Richards_State%head, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(Richards_State%kr, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(Richards_State%rho_water, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(Richards_State%mu_water, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(Richards_State%xi, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)

! allocate state variables at faces
CALL AllocateArray_1D(Richards_State%kr_faces, 1, Richards_Base%n_faces)
CALL AllocateArray_1D(Richards_State%xi_faces, 1, Richards_Base%n_faces)

! allocate hydraulic parameters
CALL AllocateArray(VGM_parameters%theta_r, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(VGM_parameters%theta_s, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(VGM_parameters%alpha, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(VGM_parameters%n, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
  
END SUBROUTINE RichardsAllocate
! ************************************************************************** !
  
SUBROUTINE RichardsDiscretize(nx, ny, nz)
USE crunchtype
USE medium
    
IMPLICIT NONE

INTEGER(I4B), INTENT(IN) :: nx, ny, nz
    
! call internal subroutine for each dimension to get cell coordinate and the width
CALL Discretize(nx, Richards_Base%x, Richards_Base%dx, x, dxx)
CALL Discretize(ny, Richards_Base%y, Richards_Base%dy, y, dyy)
   
! compute numbers of faces, cells etc.
IF (nx > 1 .AND. ny == 1 .AND. nz == 1) THEN ! one-dimensional problem
  Richards_Base%n_inner_cells = nx
  Richards_Base%n_faces = nx + 1
  Richards_Base%n_bfaces = 2
  Richards_Base%n_total_cells = Richards_Base%n_inner_cells + Richards_Base%n_bfaces

ELSE IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
  
  SELECT CASE (Richards_Base%spatial_domain)
  CASE('regular')
    Richards_Base%n_inner_cells = nx * ny
    Richards_Base%n_faces = nx*(ny+1) + (nx+1)*ny
    Richards_Base%n_bfaces = 2*(nx+ny)
    Richards_Base%n_total_cells = Richards_Base%n_inner_cells + Richards_Base%n_bfaces
  CASE DEFAULT
    WRITE(*,*)
    WRITE(*,*) ' The spatial domain ', Richards_Base%spatial_domain, ' is not supported. '
    WRITE(*,*)
    READ(*,*)
    STOP
  END SELECT
    
ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Currently, three-dimensional Richards solver is supported.'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

! create link functions
CALL CreateLinkFunction(nx, ny, nz)

CONTAINS
  SUBROUTINE Discretize(n, x_out, dx_out, x_in, dx_in)
  INTEGER(I4B), INTENT(IN) :: n
  REAL(DP), INTENT(IN) :: x_in(0:n), dx_in(-1:n+2)
  REAL(DP), ALLOCATABLE, INTENT(INOUT) :: x_out(:), dx_out(:)
  INTEGER(I4B) :: i
        
  IF (n > 1) THEN
    ! get the discretization for one-dimension
    
    CALL AllocateArray_1D(x_out, 0, n+1)
    CALL AllocateArray_1D(dx_out, 0, n+1)
    
    DO i = 1, n
      x_out(i) = x_in(i)
      dx_out(i) = dx_in(i)
    END DO
    
    ! fill the ghost cells (the width is the same as the neighboring cell)
    dx_out(0) = dx_out(1)
    x_out(0) = x_out(1) - dx_out(0)
      
    dx_out(n+1) = dx_out(n)
    x_out(n+1) = x_out(n) + dx_out(n+1)
    
  END IF
        
  END SUBROUTINE Discretize
      
END SUBROUTINE RichardsDiscretize

! ************************************************************************** !
SUBROUTINE RichardsUpdateFluid(temp)
USE crunchtype
IMPLICIT NONE

REAL(DP), ALLOCATABLE, INTENT(IN) :: temp(:,:,:) ! temperature in C

Richards_State%rho_water = 0.99823d0 * 1.0E3
Richards_State%mu_water = 0.0010005* 86400.0d0 * 365.0d0 ! dynamic viscosity of water
!Richards_State%mu_water = 10.0d0**(-4.5318d0 - 220.57d0/(149.39 - temp - 273.15d0)) * 86400.0d0 * 365.0d0 ! dynamic viscosity of water
      
END SUBROUTINE RichardsUpdateFluid
! ************************************************************************** !
  


END MODULE Richards_module