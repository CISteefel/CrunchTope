SUBROUTINE residual_Richards(nx, ny, nz, dtflow, F_residual)
! This subroutine calculates the residual of the Richards equation
! F = d theta / dt + div q - S (dimensionless)
! S is the source/sink term (positive for source)
USE crunchtype
USE io
USE params
USE runtime
USE concentration
USE medium
USE flow
USE transport

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                   :: nx
INTEGER(I4B), INTENT(IN)                                   :: ny
INTEGER(I4B), INTENT(IN)                                   :: nz

INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz
INTEGER(I4B)                                               :: n_inner_cells
INTEGER(I4B)                                               :: n_total_cells

INTEGER(I4B)                                               :: n_xfaces
INTEGER(I4B)                                               :: n_yfaces
INTEGER(I4B)                                               :: n_bfaces

REAL(DP), INTENT(IN)                                       :: dtflow
REAL(DP), ALLOCATABLE, INTENT(INOUT)                                         :: F_residual(:) ! residual


REAL(DP)                                                   :: divergence
REAL(DP)                                                   :: q_x_1
REAL(DP)                                                   :: q_x_2
REAL(DP)                                                   :: q_y_1 
REAL(DP)                                                   :: q_y_2
REAL(DP)                                                   :: dx_1
REAL(DP)                                                   :: dx_2
INTEGER(I4B)                                 :: face_ID
INTEGER(I4B), DIMENSION(4)                                 :: face_IDs
INTEGER(I4B), DIMENSION(2)                                 :: coordinate_cell ! coordinate of cell 1 for 2D problem
INTEGER(I4B), DIMENSION(2)                                 :: coordinate_cell_1 ! coordinate of cell 1 for 2D problem
INTEGER(I4B), DIMENSION(2)                                 :: coordinate_cell_2 ! coordinate of cell 1 for 2D problem

INTEGER(I4B)                                               :: jx_1
INTEGER(I4B)                                               :: jx_2
INTEGER(I4B)                                               :: jy_1
INTEGER(I4B)                                               :: jy_2

REAL(DP)                                                   :: psi_b ! water potential at a boundary
REAL(DP)                                                   :: psi_grad_b ! water potential gradient at a boundary

!REAL(DP)                                                   :: water_balance ! water balance to prevent the cell from drying out
!REAL(DP)                                                   :: adjusted_extraction ! total water extraction (=evaporation + transpiration) adjusted to prevent the cell from drying out



!*********************************************************************
IF (nx > 1 .AND. ny == 1 .AND. nz == 1) THEN ! one-dimensional problem
  
  jy = 1
  jz = 1

  ! internal cells
  DO jx = 1, nx
    F_residual(jx) = theta(jx, jy, jz) - theta_prev(jx, jy, jz) + (qx(jx, jy, jz) - qx(jx-1, jy, jz))*dtflow/dxx_2(jx)
  END DO

  ! boundary condition at the inlet (begin boundary condition)

  SELECT CASE (x_begin_BC_type)
  CASE ('constant_dirichlet', 'variable_dirichlet')
    psi_b = (psi(0, jy, jz)*dxx_2(1) + psi(1, jy, jz)*dxx_2(0))/(dxx_2(0) + dxx_2(1))
    F_residual(0) = psi_b - value_x_begin_BC
  
    CONTINUE
  
  CASE ('constant_neumann', 'variable_neumann')
    psi_grad_b = (psi(1, jy, jz) - psi(0, jy, jz))/(x_2(1) - x_2(0))
    F_residual(0) = psi_grad_b - value_x_begin_BC
  
  CASE ('constant_flux', 'variable_flux')
    F_residual(0) = qx(0, jy, jz) - value_x_begin_BC

  CASE ('constant_atomosphere', 'variable_atomosphere')
    psi_b = (psi(0, jy, jz)*dxx_2(1) + psi(1, jy, jz)*dxx_2(0))/(dxx_2(0) + dxx_2(1))
    IF (psi_b < psi_0) THEN
      ! the boundary water potential is below psi_0, so switch to Dirichlet BC
      F_residual(0) = psi_b - psi_0
    ELSE
      F_residual(0) = qx(0, jy, jz) - value_x_begin_BC
    END IF
  !CASE ('environmental_forcing')
  !  CONTINUE
  
  CASE DEFAULT
    WRITE(*,*)
    WRITE(*,*) ' The boundary condition type ', x_begin_BC_type, ' is not supported. '
    WRITE(*,*)
    READ(*,*)
    STOP
  
  END SELECT


  ! boundary condition at the inlet (end boundary condition)
  SELECT CASE (x_end_BC_type)
  CASE ('constant_dirichlet', 'variable_dirichlet')
    psi_b = (psi(nx, jy, jz)*dxx_2(nx+1) + psi(nx+1, jy, jz)*dxx_2(nx))/(dxx_2(nx) + dxx_2(nx+1))
    F_residual(nx+1) = psi_b - value_x_end_BC
  
  CASE ('constant_neumann', 'variable_neumann')
    psi_grad_b = (psi(nx+1, jy, jz) - psi(nx, jy, jz))/(x_2(nx+1) - x_2(nx))
    F_residual(nx+1) = psi_grad_b - value_x_end_BC
  
  CASE ('constant_flux', 'variable_flux')
    F_residual(nx+1) = qx(nx, jy, jz) - value_x_end_BC
  
  CASE DEFAULT
    WRITE(*,*)
    WRITE(*,*) ' The boundary condition type ', x_end_BC_type, ' is not supported. '
    WRITE(*,*)
    READ(*,*)
    STOP
  
  END SELECT

!*********************************************************************
ELSE IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
  
  IF (domain_shape_flow == 'regular') THEN
    n_inner_cells = nx * ny
    n_xfaces = nx*(ny+1) ! number of faces
    n_yfaces =  (nx+1)*ny
    n_bfaces = 2*nx + 2*ny
    n_total_cells = n_inner_cells + n_bfaces
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Currently, only regular spatial domain is supported.'
    WRITE(*,*)
    READ(*,*)
    STOP
        
  END IF
  

  IF (ALLOCATED(F_residual)) THEN
    DEALLOCATE(F_residual)
    ALLOCATE(F_residual(n_total_cells))
  ELSE
    ALLOCATE(F_residual(n_total_cells))
  END IF
  
  F_residual= 0.0
  
  ! internal cells
  DO i = 1, n_inner_cells
    face_IDs = cell_to_face(i, :)
    jx = cell_to_coordinate(i, 1)
    jy = cell_to_coordinate(i, 2)
    jz = 1
    
    q_y_1 = q_Richards(face_IDs(1))
    q_y_2 = q_Richards(face_IDs(2))
    q_x_1 = q_Richards(face_IDs(3))
    q_x_2 = q_Richards(face_IDs(4))
    divergence = (q_y_2 - q_y_1) / dyy_2(jy) + (q_x_1 - q_x_1) / dxx_2(jx)
    F_residual(i) = theta(jx, jy, jz) - theta_prev(jx, jy, jz) + divergence*dtflow
  END DO
  
  
  ! evaluate residual for boundary conditions
  DO i = 1, n_bfaces
    face_ID = bface_to_face(i)

    coordinate_cell_1 = cell_to_coordinate(face_to_cell(face_ID, 1), :)
    coordinate_cell_2 = cell_to_coordinate(face_to_cell(face_ID, 2), :)
    
    jx_1 = coordinate_cell_1(1)
    jx_2 = coordinate_cell_2(1)
    jy_1 = coordinate_cell_1(2)
    jy_2 = coordinate_cell_2(2)
    
    jz = 1
    
    SELECT CASE (BC_type_Richards(i))
    CASE (1) ! Dirichlet
      IF (jx_1 == jx_2) THEN
        ! faces paralell to x-axis
        dx_1 = dyy_2(jy_1)
        dx_2 = dyy_2(jy_2)
      ELSE
        ! faces paralell to y-axis
        dx_1 = dxx_2(jx_1)
        dx_2 = dxx_2(jx_2)
      END IF
      psi_b = (psi(jx_1, jy_1, jz)*dx_2 + psi(jx_2, jy_2, jz)*dx_1)/(dx_1 + dx_2)
      F_residual(n_inner_cells + i) = psi_b - BC_value_Richards(i)
  
    CASE (2) ! Neumann
      IF (jx_1 == jx_2) THEN
        ! faces paralell to x-axis
        dx_1 = y_2(jy_2) - y_2(jy_1) 
      ELSE
        ! faces paralell to y-axis
        dx_1 = x_2(jx_2) - x_2(jx_1)
      END IF
      psi_grad_b = (psi(jx_2, jy_2, jz) - psi(jx_1, jy_1, jz))/dx_1
      F_residual(n_inner_cells + i) = psi_grad_b - BC_value_Richards(i)
  
    CASE (3) ! flux
      F_residual(n_inner_cells + i) = q_Richards(face_ID) - BC_value_Richards(i)
  
    CASE DEFAULT
      WRITE(*,*)
      WRITE(*,*) ' The boundary condition type ', x_end_BC_type, ' is not supported. '
      WRITE(*,*)
      READ(*,*)
      STOP
  
    END SELECT
  END DO
!*********************************************************************
ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Currently, three-dimensional Richards solver is supported.'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
END SUBROUTINE residual_Richards