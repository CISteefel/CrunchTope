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
  REAL(DP) :: gravity_vector(3)
END TYPE RichardsBase
  
TYPE :: RichardsBC ! store information on the boundary condition for each face
  REAL(DP) :: BC_value
  REAL(DP), ALLOCATABLE :: BC_values(:,:) ! transient boundary condition data
  INTEGER(I4B) :: BC_type
  LOGICAL(LGT) :: is_atmosphere = .FALSE.
END TYPE RichardsBC

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
  REAL(DP) :: tau_a  = 1.0d-8 ! absolute tolerance for the Newton solver
  REAL(DP) :: max_Newton = 100 ! maximum number of Newton iterations
  REAL(DP) :: max_line_search = 30 ! maximum number of line searches
END TYPE RichardsSolver


TYPE :: RichardsState ! store information on the nonlinear solver used in Richards solver
  REAL(DP), ALLOCATABLE :: psi(:,:,:) ! water potential [L]
  REAL(DP), ALLOCATABLE :: theta(:,:,:) ! volumetric water content [L3 L-3]
  REAL(DP), ALLOCATABLE :: theta_prev(:,:,:) ! volumetric water content from previous time step [L3 L-3]
  REAL(DP), ALLOCATABLE :: head(:,:,:) ! pressure head [L]
  REAL(DP), ALLOCATABLE :: kr(:,:,:) ! relative permeability [-]
  REAL(DP), ALLOCATABLE :: dtheta(:,:,:) ! the derivative of the water content wirh respect to the water potential [-]
  REAL(DP), ALLOCATABLE :: dkr(:,:,:) ! the derivative of the relative permeability  wirh respect to the water potential [-]
  REAL(DP), ALLOCATABLE :: xi(:,:,:) ! physical constant used in the Richards equation
  REAL(DP), ALLOCATABLE :: kr_faces(:) ! relative permeability at cell faces [-]
  REAL(DP), ALLOCATABLE :: xi_faces(:) ! physical constant in the Richards equation at faces    
END TYPE RichardsState

TYPE(RichardsBase), PUBLIC :: Richards_Base
TYPE(RichardsBC), ALLOCATABLE, PUBLIC :: Richards_BCs(:)
TYPE(RichardsBC), ALLOCATABLE, PUBLIC :: Richards_BCs_steady(:)
TYPE(RichardsOptions), PUBLIC :: Richards_Options
TYPE(RichardsSolver), PUBLIC :: Richards_Solver
TYPE(RichardsState), PUBLIC :: Richards_State
  
PUBLIC RichardsAllocate, &
       RichardsDiscretize, &
       RichardsFlux, &
       RichardsSolve, &
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
SUBROUTINE AllocateArray_2D(array, lbound_x, ubound_x, lbound_y, ubound_y)
USE crunchtype
IMPLICIT NONE

REAL(DP), ALLOCATABLE, INTENT(INOUT) :: array(:,:)
INTEGER(I4B), INTENT(IN) :: lbound_x, ubound_x, lbound_y, ubound_y
  
! Check if the array is allocated, deallocate if necessary, then allocate with the new dimensions
IF (ALLOCATED(array)) THEN
  DEALLOCATE(array)
  ALLOCATE(array(lbound_x:ubound_x, lbound_y:ubound_y))
ELSE
  ALLOCATE(array(lbound_x:ubound_x, lbound_y:ubound_y))
END IF
   
END SUBROUTINE AllocateArray_2D
! ************************************************************************** !
REAL(DP) FUNCTION ArithmaticMean(S_1, S_2, dx_1, dx_2)
  
REAL(DP), INTENT(IN) :: S_1, S_2, dx_1, dx_2
  
ArithmaticMean = (S_1*dx_2 + S_2*dx_1)/(dx_1 + dx_2)
   
END FUNCTION ArithmaticMean
! ************************************************************************** !
SUBROUTINE CreateLinkFunction(nx, ny, nz)
USE crunchtype

IMPLICIT NONE

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
CALL AllocateArray(Richards_State%dtheta, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
CALL AllocateArray(Richards_State%dkr, lbound_x, ubound_x, lbound_y, ubound_y, lbound_z, ubound_z)
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
SUBROUTINE RichardsFlux(nx, ny, nz)
USE crunchtype
USE flow, ONLY: harx, hary, harz
USE transport, ONLY: qx, qy, qz

IMPLICIT NONE

INTEGER(I4B), INTENT(IN) :: nx, ny, nz
INTEGER(I4B) :: i, jx, jy, jz
INTEGER(I4B) :: jx_west, jy_west, jx_east, jy_east
INTEGER(I4B) :: cell_ID_west, cell_ID_east
REAL(DP) :: temp_theta, temp_kr, temp_dtheta, temp_dkr
REAL(DP) :: dx_west, dx_east, delta_x
REAL(DP) :: gravity
REAL(DP) :: psi_grad ! gradient of water potential
REAL(DP) :: q_diff ! diffusion flow
REAL(DP) :: q_grav ! gravitational flow
REAL(DP) :: K_face ! permeability at face

IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
! apply van Genuchten model to all grid cells
  jz = 1
  DO jy = 0, ny + 1
    DO jx = 0, nx + 1
      CALL VGM_Model(Richards_State%psi(jx, jy, jz), &
                     VGM_parameters%theta_r(jx, jy, jz), &
                     VGM_parameters%theta_s(jx, jy, jz), &
                     VGM_parameters%alpha(jx, jy, jz), &
                     VGM_parameters%n(jx, jy, jz), &
                     temp_theta, temp_kr, temp_dtheta, temp_dkr, &
                     .FALSE.)
      
      ! Assign temporary results back to derived type
      Richards_State%theta(jx, jy, jz) = temp_theta
      Richards_State%kr(jx, jy, jz) = temp_kr
    
      Richards_State%head(jx, jy, jz) = Richards_State%psi(jx, jy, jz) + & 
                                        Richards_Base%gravity_vector(1)*Richards_Base%x(jx) + &
                                        Richards_Base%gravity_vector(2)*Richards_Base%y(jy)
    END DO
  END DO
  
  !compute flux by looping over faces
  DO i = 1, Richards_Base%n_faces
    cell_ID_west = Richards_Base%face_to_cell(i, 1) ! south or west cell
    cell_ID_east = Richards_Base%face_to_cell(i, 2) ! south or west cell
    
    jx_west = Richards_Base%cell_to_coordinate(cell_ID_west, 1) ! x-coordinate of the west (or south) cell
    jy_west = Richards_Base%cell_to_coordinate(cell_ID_west, 2) ! y-coordinate of the west (or south) cell
    jx_east = Richards_Base%cell_to_coordinate(cell_ID_east, 1) ! x-coordinate of the east (or north) cell
    jy_east = Richards_Base%cell_to_coordinate(cell_ID_east, 2) ! y-coordinate of the east (or north) cell

    IF (jx_west == jx_east) THEN
      ! faces paralell to x-axis
      dx_west = Richards_Base%dy(jy_west)
      dx_east = Richards_Base%dy(jy_east)
      delta_x = Richards_Base%y(jy_east) - Richards_Base%y(jy_west)
      K_face = hary(jx_west, jy_west, jz)
      gravity = Richards_Base%gravity_vector(2)
    ELSE
      ! faces paralell to y-axis
      dx_west = Richards_Base%dx(jx_west)
      dx_east = Richards_Base%dx(jx_east)
      delta_x = Richards_Base%x(jx_east) - Richards_Base%x(jx_west)
      K_face = harx(jx_west, jy_west, jz)
      gravity = Richards_Base%gravity_vector(1)
    END IF

    ! distance weighted arithmatic mean for face values
    Richards_State%xi_faces(i) = ArithmaticMean(Richards_State%xi(jx_west, jy_west, jz),&
                                                Richards_State%xi(jx_east, jy_east, jz),&
                                                dx_west, dx_east)
    
    Richards_State%kr_faces(i) = ArithmaticMean(Richards_State%kr(jx_west, jy_west, jz),&
                                                Richards_State%kr(jx_east, jy_east, jz),&
                                                dx_west, dx_east)
    
    ! compute fluxes
    
    psi_grad = (Richards_State%psi(jx_east, jy_east, jz) - Richards_State%psi(jx_west, jy_west, jz))/delta_x
    q_diff = -Richards_State%xi_faces(i)*K_face*Richards_State%kr_faces(i)*psi_grad
    q_grav = -gravity*K_face* &
              MERGE(Richards_State%kr(jx_east, jy_east, jz)*Richards_State%xi(jx_east, jy_east, jz), &
                    Richards_State%kr(jx_west, jy_west, jz)*Richards_State%xi(jx_west, jy_west, jz), &
                    Richards_State%head(jx_east, jy_east, jz) - Richards_State%head(jx_west, jy_west, jz) >= 0.0D0)
    
    ! store flux in the flux variable
    IF (jx_west == jx_east) THEN
      ! faces paralell to x-axis
      qy(jx_west, jy_west, jz) = q_diff + q_grav
    ELSE
      ! faces paralell to y-axis
      qx(jx_west, jy_west, jz) = q_diff + q_grav
    END IF
        
  END DO
  
END IF

END SUBROUTINE RichardsFlux
! ************************************************************************** !
SUBROUTINE RichardsJacobian(nx, ny, nz, dt, J)
USE crunchtype
USE flow, ONLY: harx, hary, harz

INTEGER(I4B), INTENT(IN) :: nx, ny, nz
REAL(DP), INTENT(IN) :: dt ! time step

REAL(DP), ALLOCATABLE, INTENT(INOUT) :: J(:,:)
REAL(DP) :: temp_theta, temp_kr, temp_dtheta, temp_dkr
INTEGER(I4B) :: i
INTEGER(I4B) :: jx, jy, jz
INTEGER(I4B) :: jx_inner, jy_inner, jx_outer, jy_outer
INTEGER(I4B) :: face_ID
INTEGER(I4B) :: face_IDs(4)
INTEGER(I4B) :: face_ID_south
INTEGER(I4B) :: face_ID_north
INTEGER(I4B) :: face_ID_west
INTEGER(I4B) :: face_ID_east
  
INTEGER(I4B) :: cell_ID_south
INTEGER(I4B) :: cell_ID_north
INTEGER(I4B) :: cell_ID_west
INTEGER(I4B) :: cell_ID_east
INTEGER(I4B) :: cell_ID_inner
INTEGER(I4B) :: cell_ID_outer

REAL(DP) :: delta_x
REAL(DP) :: gravity
REAL(DP) :: K_face ! permeability at face    
    
! allocate Jacobian
CALL AllocateArray_2D(J, 1, Richards_Base%n_total_cells, 1, Richards_Base%n_total_cells)

! initialize
J= 0.0d0

IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
  
  ! allocate arrays for derivative
  
  
  ! apply van Genuchten model to all grid cells
  jz = 1
  DO jy = 0, ny + 1
    DO jx = 0, nx + 1
      CALL VGM_Model(Richards_State%psi(jx, jy, jz), &
                     VGM_parameters%theta_r(jx, jy, jz), &
                     VGM_parameters%theta_s(jx, jy, jz), &
                     VGM_parameters%alpha(jx, jy, jz), &
                     VGM_parameters%n(jx, jy, jz), &
                     temp_theta, temp_kr, temp_dtheta, temp_dkr, &
                     .TRUE.)
      
      ! Assign temporary results back to derived type
      Richards_State%dtheta(jx, jy, jz) = temp_dtheta
      Richards_State%dkr(jx, jy, jz) = temp_dkr
    

    END DO
  END DO
  
  ! internal cells
  DO i = 1, Richards_Base%n_inner_cells
    face_IDs = Richards_Base%cell_to_face(i, :)
    
    face_ID_south = face_IDs(1)
    face_ID_north = face_IDs(2)
    face_ID_west = face_IDs(3)
    face_ID_east = face_IDs(4)
    
    cell_ID_south = Richards_Base%face_to_cell(face_ID_south, 1)
    cell_ID_north = Richards_Base%face_to_cell(face_ID_north, 2)
    cell_ID_west = Richards_Base%face_to_cell(face_ID_west, 1)
    cell_ID_east = Richards_Base%face_to_cell(face_ID_east, 2)
    
    jx = Richards_Base%cell_to_coordinate(i, 1)
    jy = Richards_Base%cell_to_coordinate(i, 2)
    
    ! add temporal derivative part
    J(i, i) = J(i, i) + Richards_State%dtheta(jx, jy, jz)
    
    ! add diffusive flux part
        
    !! center cell
    !!! due to q at west face
    J(i, i) = J(i, i) + dt/Richards_Base%dx(jx) * harx(jx-1, jy, jz) * Richards_State%xi_faces(face_ID_west) / (Richards_Base%x(jx) - Richards_Base%x(jx-1)) * &
                  (Richards_State%dkr(jx, jy, jz)*Richards_Base%dx(jx-1)/(Richards_Base%dx(jx) + Richards_Base%dx(jx-1))*(Richards_State%psi(jx, jy, jz) - Richards_State%psi(jx-1, jy, jz)) + Richards_State%kr_faces(face_ID_west))
  
    !!! due to q at east face
    J(i, i) = J(i, i) - dt/Richards_Base%dx(jx) * harx(jx, jy, jz) * Richards_State%xi_faces(face_ID_east) / (Richards_Base%x(jx+1) - Richards_Base%x(jx)) * &
                  (Richards_State%dkr(jx, jy, jz)*Richards_Base%dx(jx+1)/(Richards_Base%dx(jx+1) + Richards_Base%dx(jx))*(Richards_State%psi(jx+1, jy, jz) - Richards_State%psi(jx, jy, jz)) - Richards_State%kr_faces(face_ID_east))
  
    !!! due to q at south face
    J(i, i) = J(i, i) + dt/Richards_Base%dy(jy) * hary(jx, jy-1, jz) * Richards_State%xi_faces(face_ID_south) / (Richards_Base%y(jy) - Richards_Base%y(jy-1)) * &
                  (Richards_State%dkr(jx, jy, jz)*Richards_Base%dy(jy-1)/(Richards_Base%dy(jy) + Richards_Base%dy(jy-1))*(Richards_State%psi(jx, jy, jz) - Richards_State%psi(jx, jy-1, jz)) + Richards_State%kr_faces(face_ID_south))
    
    !!! due to q at east face
    J(i, i) = J(i, i) - dt/Richards_Base%dy(jy) * hary(jx, jy, jz) * Richards_State%xi_faces(face_ID_north) / (Richards_Base%y(jy+1) - Richards_Base%y(jy)) * &
                  (Richards_State%dkr(jx, jy, jz)*Richards_Base%dy(jy+1)/(Richards_Base%dy(jy+1) + Richards_Base%dy(jy))*(Richards_State%psi(jx, jy+1, jz) - Richards_State%psi(jx, jy, jz)) - Richards_State%kr_faces(face_ID_north))
    
    !! west cell due to q at west face
    J(i, cell_ID_west) = J(i, cell_ID_west) + dt/Richards_Base%dx(jx) * harx(jx-1, jy, jz) * Richards_State%xi_faces(face_ID_west) / (Richards_Base%x(jx) - Richards_Base%x(jx-1)) * &
                  (Richards_State%dkr(jx-1, jy, jz)*Richards_Base%dx(jx)/(Richards_Base%dx(jx) + Richards_Base%dx(jx-1))*(Richards_State%psi(jx, jy, jz) - Richards_State%psi(jx-1, jy, jz)) - Richards_State%kr_faces(face_ID_west))
    
    !! east cell due to q at east face
    J(i, cell_ID_east) = J(i, cell_ID_east) - dt/Richards_Base%dx(jx) * harx(jx, jy, jz) * Richards_State%xi_faces(face_ID_east) / (Richards_Base%x(jx+1) - Richards_Base%x(jx)) * &
                  (Richards_State%dkr(jx+1, jy, jz)*Richards_Base%dx(jx)/(Richards_Base%dx(jx+1) + Richards_Base%dx(jx))*(Richards_State%psi(jx+1, jy, jz) - Richards_State%psi(jx, jy, jz)) + Richards_State%kr_faces(face_ID_east))
    
    !! south cell due to q at south face
    J(i, cell_ID_south) = J(i, cell_ID_south) + dt/Richards_Base%dy(jy) * hary(jx, jy-1, jz) * Richards_State%xi_faces(face_ID_south) / (Richards_Base%y(jy) - Richards_Base%y(jy-1)) * &
                  (Richards_State%dkr(jx, jy-1, jz)*Richards_Base%dy(jy)/(Richards_Base%dy(jy) + Richards_Base%dy(jy-1))*(Richards_State%psi(jx, jy, jz) - Richards_State%psi(jx, jy-1, jz)) - Richards_State%kr_faces(face_ID_south))
    
    !! north cell due to q at north face
    J(i, cell_ID_north) = J(i, cell_ID_north) - dt/Richards_Base%dy(jy) * hary(jx, jy, jz) * Richards_State%xi_faces(face_ID_north) / (Richards_Base%y(jy+1) - Richards_Base%y(jy)) * &
                  (Richards_State%dkr(jx, jy+1, jz)*Richards_Base%dy(jy)/(Richards_Base%dy(jy+1) + Richards_Base%dy(jy))*(Richards_State%psi(jx, jy+1, jz) - Richards_State%psi(jx, jy, jz)) + Richards_State%kr_faces(face_ID_north))
  
    ! add gravitational flux part
    
    !! center cell
    !!! due to q at west face
    J(i, i) = J(i, i) + dt/Richards_Base%dx(jx) * Richards_Base%gravity_vector(1)*harx(jx-1, jy, jz)*MERGE(Richards_State%dkr(jx, jy, jz) * Richards_State%xi(jx, jy, jz), 0.0d0, Richards_State%head(jx, jy, jz) - Richards_State%head(jx-1, jy, jz) >= 0.0d0)
  
    !!! due to q at east face
    J(i, i) = J(i, i) - dt/Richards_Base%dx(jx) * Richards_Base%gravity_vector(1)*harx(jx, jy, jz)*MERGE(0.0d0, Richards_State%dkr(jx, jy, jz) * Richards_State%xi(jx, jy, jz), Richards_State%head(jx+1, jy, jz) - Richards_State%head(jx, jy, jz) >= 0.0d0)
  
    !!! due to q at south face
    J(i, i) = J(i, i) + dt/Richards_Base%dy(jy) * Richards_Base%gravity_vector(2)*hary(jx, jy-1, jz)*MERGE(Richards_State%dkr(jx, jy, jz) * Richards_State%xi(jx, jy, jz), 0.0d0, Richards_State%head(jx, jy, jz) - Richards_State%head(jx, jy-1, jz) >= 0.0d0)
  
    !!! due to q at north face
    J(i, i) = J(i, i) - dt/Richards_Base%dy(jy) * Richards_Base%gravity_vector(2)*hary(jx, jy, jz)*MERGE(0.0d0, Richards_State%dkr(jx, jy, jz) * Richards_State%xi(jx, jy, jz), Richards_State%head(jx, jy+1, jz) - Richards_State%head(jx, jy, jz) >= 0.0d0)
  
    !! west cell due to q at west face
    J(i, cell_ID_west) = J(i, cell_ID_west) + dt/Richards_Base%dx(jx) * Richards_Base%gravity_vector(1)*harx(jx-1, jy, jz)*MERGE(0.0d0, Richards_State%dkr(jx-1, jy, jz) * Richards_State%xi(jx-1, jy, jz), Richards_State%head(jx, jy, jz) - Richards_State%head(jx-1, jy, jz) >= 0.0d0)
  
    !! east cell due to q at east face
    J(i, cell_ID_east) = J(i, cell_ID_east) - dt/Richards_Base%dx(jx) * Richards_Base%gravity_vector(1)*harx(jx, jy, jz)*MERGE(Richards_State%dkr(jx+1, jy, jz) * Richards_State%xi(jx+1, jy, jz), 0.0d0, Richards_State%head(jx+1, jy, jz) - Richards_State%head(jx, jy, jz) >= 0.0d0)
  
    !! south cell due to q at south face
    J(i, cell_ID_south) = J(i, cell_ID_south) + dt/Richards_Base%dy(jy) * Richards_Base%gravity_vector(2)*hary(jx, jy-1, jz)*MERGE(0.0d0, Richards_State%dkr(jx, jy-1, jz) * Richards_State%xi(jx, jy-1, jz), Richards_State%head(jx, jy, jz) - Richards_State%head(jx, jy-1, jz) >= 0.0d0)
    
    !! north cell due to q at north face
    J(i, cell_ID_north) = J(i, cell_ID_north) - dt/Richards_Base%dy(jy) * Richards_Base%gravity_vector(2)*hary(jx, jy, jz)*MERGE(Richards_State%dkr(jx, jy+1, jz) * Richards_State%xi(jx, jy+1, jz), 0.0d0, Richards_State%head(jx, jy+1, jz) - Richards_State%head(jx, jy, jz) >= 0.0d0)
  
  END DO
  
  ! evaluate residual for boundary conditions
  DO i = 1, Richards_Base%n_bfaces
    face_ID = Richards_Base%bface_to_face(i)
    
    cell_ID_inner = MINVAL(Richards_Base%face_to_cell(face_ID, :)) ! inner cell
    cell_ID_outer = MAXVAL(Richards_Base%face_to_cell(face_ID, :)) ! ghost cell
    
    jx_inner = Richards_Base%cell_to_coordinate(cell_ID_inner, 1) ! x-coordinate of the inner cell
    jy_inner = Richards_Base%cell_to_coordinate(cell_ID_inner, 2) ! y-coordinate of the inner cell
    jx_outer = Richards_Base%cell_to_coordinate(cell_ID_outer, 1) ! x-coordinate of the ghost cell
    jy_outer = Richards_Base%cell_to_coordinate(cell_ID_outer, 2) ! y-coordinate of the ghost cell
      
    SELECT CASE (Richards_BCs(i)%BC_type)
    CASE (1) ! Dirichlet
      J(cell_ID_outer, cell_ID_outer) = 0.5d0
      J(cell_ID_outer, cell_ID_inner) = 0.5d0
  
    CASE (2) ! Neumann
      IF (jx_inner == jx_outer) THEN
        ! faces paralell to x-axis
        delta_x = Richards_Base%y(jy_outer) - Richards_Base%y(jy_inner)
      ELSE
        ! faces paralell to y-axis
        delta_x = Richards_Base%x(jx_outer) - Richards_Base%x(jx_inner)
      END IF
      J(cell_ID_outer, cell_ID_outer) = 1.0d0/delta_x
      J(cell_ID_outer, cell_ID_inner) = -1.0d0/delta_x
  
    CASE (3) ! flux
      IF (jx_inner == jx_outer) THEN
        ! faces paralell to x-axis
        delta_x = Richards_Base%y(jy_outer) - Richards_Base%y(jy_inner)
        gravity = Richards_Base%gravity_vector(2)
        IF (jy_inner < jy_outer) THEN
          K_face = hary(jx_inner, jy_inner, jz)
        ELSE
          K_face = hary(jx_outer, jy_outer, jz)    
        END IF
      ELSE
        ! faces paralell to y-axis
        delta_x = Richards_Base%x(jx_outer) - Richards_Base%x(jx_inner)
        gravity = Richards_Base%gravity_vector(1)
        IF (jx_inner < jx_outer) THEN
          K_face = harx(jx_inner, jy_inner, jz)
        ELSE
          K_face = harx(jx_outer, jy_outer, jz)    
        END IF
      END IF
   
      ! add diffusive flux part
  
      J(cell_ID_outer, cell_ID_outer) = J(cell_ID_outer, cell_ID_outer) - K_face * Richards_State%xi_faces(face_ID) / delta_x * &
                    (Richards_State%dkr(jx_outer, jy_outer, jz)*0.5d0*(Richards_State%psi(jx_outer, jy_outer, jz) - Richards_State%psi(jx_inner, jy_inner, jz)) + Richards_State%kr_faces(face_ID))
      
      J(cell_ID_outer, cell_ID_inner) = J(cell_ID_outer, cell_ID_inner) - K_face * Richards_State%xi_faces(face_ID) / delta_x * &
                    (Richards_State%dkr(jx_inner, jy_inner, jz)*0.5d0*(Richards_State%psi(jx_outer, jy_outer, jz) - Richards_State%psi(jx_inner, jy_inner, jz)) - Richards_State%kr_faces(face_ID))
  
      ! add gravitational flux part
      J(cell_ID_outer, cell_ID_outer) = J(cell_ID_outer, cell_ID_outer) - gravity*K_face*MERGE(Richards_State%dkr(jx_outer, jy_outer, jz) * Richards_State%xi(jx_outer, jy_outer, jz), 0.0d0, Richards_State%head(jx_outer, jy_outer, jz) - Richards_State%head(jx_inner, jy_inner, jz) >= 0.0d0)
  
      J(cell_ID_outer, cell_ID_inner) = J(cell_ID_outer, cell_ID_inner) - gravity*K_face*MERGE(0.0d0, Richards_State%dkr(jx_inner, jy_inner, jz) * Richards_State%xi(jx_inner, jy_inner, jz), Richards_State%head(jx_outer, jy_outer, jz) - Richards_State%head(jx_inner, jy_inner, jz) >= 0.0d0)
      
    CASE DEFAULT
      WRITE(*,*)
      WRITE(*,*) ' The boundary condition type ', Richards_BCs(i)%BC_type, ' is not recognized. '
      WRITE(*,*)
      READ(*,*)
      STOP
  
    END SELECT
  END DO
  
END IF

END SUBROUTINE RichardsJacobian
! ************************************************************************** !
SUBROUTINE RichardsResidual(nx, ny, nz, dt, F_residual)
USE crunchtype
USE transport, ONLY: qx, qy, qz

INTEGER(I4B), INTENT(IN) :: nx, ny, nz
REAL(DP), INTENT(IN) :: dt ! time step

REAL(DP), ALLOCATABLE, INTENT(INOUT) :: F_residual(:)

INTEGER(I4B) :: i
INTEGER(I4B) :: jx, jy, jz
REAL(DP) :: delta_x
REAL(DP) :: divergence
REAL(DP) :: psi_b ! water potential at a boundary (this is also used for the gradient)
INTEGER(I4B) :: face_ID
INTEGER(I4B) :: cell_ID_inner, cell_ID_outer
INTEGER(I4B) :: jx_inner, jy_inner, jx_outer, jy_outer

! allocate residual
CALL AllocateArray_1D(F_residual, 1, Richards_Base%n_total_cells)

! initialize
F_residual= 0.0d0

IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
  
  jz = 1
  
  ! internal cells
  DO i = 1, Richards_Base%n_inner_cells
    jx = Richards_Base%cell_to_coordinate(i, 1)
    jy = Richards_Base%cell_to_coordinate(i, 2)
        
    divergence = (qy(jx, jy, jz) - qy(jx, jy-1, jz)) / Richards_Base%dy(jy) + &
                 (qx(jx, jy, jz) - qx(jx-1, jy, jz)) / Richards_Base%dx(jx)
    
    F_residual(i) = Richards_State%theta(jx, jy, jz) - Richards_State%theta_prev(jx, jy, jz) + divergence*dt
    
  END DO
  
  ! evaluate residual for boundary conditions
  DO i = 1, Richards_Base%n_bfaces
    face_ID = Richards_Base%bface_to_face(i)
    
    cell_ID_inner = MINVAL(Richards_Base%face_to_cell(face_ID, :)) ! inner cell
    cell_ID_outer = MAXVAL(Richards_Base%face_to_cell(face_ID, :)) ! ghost cell
    
    jx_inner = Richards_Base%cell_to_coordinate(cell_ID_inner, 1) ! x-coordinate of the inner cell
    jy_inner = Richards_Base%cell_to_coordinate(cell_ID_inner, 2) ! y-coordinate of the inner cell
    jx_outer = Richards_Base%cell_to_coordinate(cell_ID_outer, 1) ! x-coordinate of the ghost cell
    jy_outer = Richards_Base%cell_to_coordinate(cell_ID_outer, 2) ! y-coordinate of the ghost cell
      
    SELECT CASE (Richards_BCs(i)%BC_type)
    CASE (1) ! Dirichlet
      ! water potential at the boundary
      psi_b = (Richards_State%psi(jx_inner, jy_inner, jz) + Richards_State%psi(jx_outer, jy_outer, jz))/2.0d0
      
      F_residual(Richards_Base%n_inner_cells + i) = psi_b - Richards_BCs(i)%BC_value
  
    CASE (2) ! Neumann
      IF (jx_inner == jx_outer) THEN
        ! faces paralell to x-axis
        delta_x = Richards_Base%y(jy_outer) - Richards_Base%y(jy_inner)
      ELSE
        ! faces paralell to y-axis
        delta_x = Richards_Base%x(jx_outer) - Richards_Base%x(jx_inner)
      END IF
      ! water potential gradient at the boundary
      psi_b = (Richards_State%psi(jx_outer, jy_outer, jz) - Richards_State%psi(jx_inner, jy_inner, jz))/delta_x
      F_residual(Richards_Base%n_inner_cells + i) = psi_b - Richards_BCs(i)%BC_value
  
    CASE (3) ! flux
      IF (jx_inner == jx_outer) THEN
        ! faces paralell to x-axis
        IF (jy_inner < jy_outer) THEN
          F_residual(Richards_Base%n_inner_cells + i) = qy(jx_inner, jy_inner, jz) - Richards_BCs(i)%BC_value    
        ELSE
          F_residual(Richards_Base%n_inner_cells + i) = qy(jx_outer, jy_outer, jz) - Richards_BCs(i)%BC_value    
        END IF
        
      ELSE
        ! faces paralell to y-axis
        IF (jx_inner < jx_outer) THEN
          F_residual(Richards_Base%n_inner_cells + i) = qx(jx_inner, jy_inner, jz) - Richards_BCs(i)%BC_value
        ELSE
          F_residual(Richards_Base%n_inner_cells + i) = qx(jx_outer, jy_outer, jz) - Richards_BCs(i)%BC_value
        END IF
          
      END IF
          
    CASE DEFAULT
      WRITE(*,*)
      WRITE(*,*) ' The boundary condition type ', Richards_BCs(i)%BC_type, ' is not recognized. '
      WRITE(*,*)
      READ(*,*)
      STOP
  
    END SELECT
  END DO
  
  
END IF

END SUBROUTINE RichardsResidual
! ************************************************************************** !
SUBROUTINE RichardsSolve(nx, ny, nz, dt)
USE crunchtype

IMPLICIT NONE
INTEGER(I4B), INTENT(IN) :: nx, ny, nz
REAL(DP), INTENT(IN) :: dt ! time step

INTEGER(I4B) :: i, jx, jy, jz

REAL(DP), ALLOCATABLE :: F_residual(:) ! residual
REAL(DP), ALLOCATABLE :: J(:, :) ! Jacobian matrix
REAL(DP), ALLOCATABLE :: dpsi_Newton(:) ! Newton step

REAL(DP), ALLOCATABLE :: psi_prev(:,:,:) ! old solution for line searches

! parameters for the linear solver
INTEGER(I4B) :: N, info, lda, ldb, nrhs
INTEGER(I4B), ALLOCATABLE :: ipiv(:)

! parameters for Newtons' method for forward problem
REAL(DP), PARAMETER                                        :: tau_a = 1.0d-8

INTEGER(I4B) :: maxitr_Newton = 100
INTEGER(I4B) :: maxitr_line_search = 30
! variables for line search
REAL(DP) :: error_old, tol, alpha_line, lam, error_new
INTEGER(I4B) :: no_backtrack, descent
INTEGER(I4B) :: iteration
INTEGER(I4B) :: total_line

! allocate arrays

CALL AllocateArray_2D(J, 1, Richards_Base%n_total_cells, 1, Richards_Base%n_total_cells)
CALL AllocateArray_1D(dpsi_Newton, 1, Richards_Base%n_total_cells)

IF (ALLOCATED(ipiv)) THEN
  DEALLOCATE(ipiv)
  ALLOCATE(ipiv(Richards_Base%n_total_cells))
ELSE
  ALLOCATE(ipiv(Richards_Base%n_total_cells))
END IF
  
! set parameters for linear solver
nrhs = 1
lda = Richards_Base%n_total_cells
ldb = Richards_Base%n_total_cells
N = Richards_Base%n_total_cells

! initialize parameters for Newton's method
iteration = 0
total_line = 0
error_old = 1.0d20 ! at least one Newton iteration will be conducted
tol = Richards_Solver%tau_a

! begin Newton's method
newton_loop: DO
  IF (error_old < tol) EXIT
  ! Evaluate the flux and the residual
  CALL RichardsFlux(nx, ny, nz)
  CALL RichardsResidual(nx, ny, nz, dt, F_residual)
  
  ! Evaluate the Jacobian matrix
  CALL RichardsJacobian(nx, ny, nz, dt, J)
  
  dpsi_Newton = -F_residual
  
  ! Solve the linear system
  CALL dgesv(N, nrhs, J, lda, ipiv, dpsi_Newton, ldb, info)
  
  ! Armijo line search
  alpha_line = 1.0d-4
  lam = 1.0d0
  descent = 0
  no_backtrack = 0
  psi_prev = Richards_State%psi
  
  ! line search
  line: DO
    IF (descent /= 0 .OR. no_backtrack > maxitr_line_search) EXIT line
    ! update water potential
    
    IF (nx > 1 .AND. ny > 1) THEN ! two-dimensional problem or three-dimensional problem
      DO i = 1, Richards_Base%n_total_cells
        jx = Richards_Base%cell_to_coordinate(i, 1)
        jy = Richards_Base%cell_to_coordinate(i, 2)
        jz = Richards_Base%cell_to_coordinate(i, 3)
        IF (ABS(lam * dpsi_Newton(i)) > Richards_Solver%dpsi_max) THEN
          Richards_State%psi(jx, jy, jz) = Richards_State%psi(jx, jy, jz) + SIGN(Richards_Solver%dpsi_max, dpsi_Newton(i))
        ELSE
          Richards_State%psi(jx, jy, jz) = Richards_State%psi(jx, jy, jz) + lam * dpsi_Newton(i)
        END IF
      END DO
    END IF
    
    ! evaluate the flux and the residual
    CALL RichardsFlux(nx, ny, nz)
    CALL RichardsResidual(nx, ny, nz, dt, F_residual)
  
    ! update tolerance
    error_new = MAXVAL(ABS(F_residual))
    ! Check if Armijo conditions are satisfied
    Armijo: IF (error_new < error_old - error_old*alpha_line*lam) THEN
      error_old = error_new
      descent = 1
    ELSE Armijo
      lam = 0.5d0 * lam
      no_backtrack = no_backtrack + 1
      Richards_State%psi = psi_prev
          
    END IF Armijo
        
  END DO line
  
  IF (iteration > Richards_Solver%max_Newton) THEN
    WRITE(*,*) ' The Newton method failed to converge in the Richards solver. '
    READ(*,*)
    STOP
  END IF
  
  total_line = total_line + no_backtrack
  iteration = iteration + 1
  
  IF (Richards_Options%is_print) THEN
    WRITE(*,120) iteration, tol, error_old
    120 FORMAT(1X, 'At the', I3, ' th Newton iteration, the tolerance is  ', ES14.4, ' , and the error is ', ES20.8)
  END IF
  
END DO newton_loop

END SUBROUTINE RichardsSolve
! ************************************************************************** !
ELEMENTAL SUBROUTINE RichardsUpdateFluid(temp, xi)
USE crunchtype
USE params, ONLY: g

IMPLICIT NONE

REAL(DP), INTENT(IN) :: temp ! temperature in C
REAL(DP), INTENT(INOUT) :: xi ! temperature in C
REAL(DP) :: rho_water
REAL(DP) :: mu_water

rho_water = 0.99823d0 * 1.0E3
mu_water = 0.0010005* 86400.0d0 * 365.0d0 ! dynamic viscosity of water
!mu_water = 10.0d0**(-4.5318d0 - 220.57d0/(149.39 - temp - 273.15d0)) * 86400.0d0 * 365.0d0 ! dynamic viscosity of water

xi = rho_water*g/mu_water

END SUBROUTINE RichardsUpdateFluid
! ************************************************************************** !
  


END MODULE Richards_module