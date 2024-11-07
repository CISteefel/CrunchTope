SUBROUTINE Jacobian_Richards(nx, ny, nz, dtflow, low_bound_F, high_bound_F, J)
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
INTEGER(I4B)                                               :: n_inner_cells
INTEGER(I4B)                                               :: n_xfaces
INTEGER(I4B)                                               :: n_yfaces
INTEGER(I4B)                                               :: n_bfaces
INTEGER(I4B), INTENT(IN)                                   :: low_bound_F
INTEGER(I4B), INTENT(IN)                                   :: high_bound_F

REAL(DP), INTENT(OUT), DIMENSION(low_bound_F:high_bound_F, low_bound_F:high_bound_F)                   :: J ! Jacobian matrix
INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

INTEGER(I4B)                                               :: jx_1
INTEGER(I4B)                                               :: jy_1
INTEGER(I4B)                                               :: jx_2
INTEGER(I4B)                                               :: jy_2

REAL(DP)                                                   :: dx_1
REAL(DP)                                                   :: dx_2
REAL(DP)                                                   :: gravity_effect
INTEGER(I4B)                                 :: face_ID
INTEGER(I4B)                                 :: face_ID_south
INTEGER(I4B)                                 :: face_ID_north
INTEGER(I4B)                                 :: face_ID_west
INTEGER(I4B)                                 :: face_ID_east
  
INTEGER(I4B)                                 :: cell_ID_south
INTEGER(I4B)                                 :: cell_ID_north
INTEGER(I4B)                                 :: cell_ID_west
INTEGER(I4B)                                 :: cell_ID_east

INTEGER(I4B), DIMENSION(4)                                 :: face_IDs
INTEGER(I4B), DIMENSION(2)                                 :: cell_IDs
INTEGER(I4B)                                 :: cell_ID_inner
INTEGER(I4B)                                 :: cell_ID_outer
INTEGER(I4B), DIMENSION(2)                                 :: coordinate_cell ! coordinate of cell 1 for 2D problem
INTEGER(I4B), DIMENSION(2)                                 :: coordinate_cell_inner ! coordinate of cell 1 for 2D problem
INTEGER(I4B), DIMENSION(2)                                 :: coordinate_cell_outer ! coordinate of cell 1 for 2D problem

REAL(DP), INTENT(IN)                                       :: dtflow
REAL(DP)                                                   :: psi_b ! water potential at a boundary

xi_2 = rho_water_2*g/mu_water
J = 0.0 ! initialize Jacobian matrix

!*********************************************************************
IF (nx > 1 .AND. ny == 1 .AND. nz == 1) THEN ! one-dimensional problem
  jy = 1
  jz = 1

  ! evaluate Jacobian matrix
  ! interior cells
  inner: DO jx = 1, nx  
    ! add diffusive flux part
    !! q at jx-1 part
    J(jx, jx - 1) = J(jx, jx - 1) + dtflow/dxx_2(jx) * K_faces(jx-1) * xi_2_faces(jx-1) / (x_2(jx) - x_2(jx-1)) * &
                  (dkr(jx-1, jy, jz)*dxx_2(jx)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) - kr_faces(jx-1))
  
    !! q at jx-1 part
    J(jx, jx) = J(jx, jx) + dtflow/dxx_2(jx) * K_faces(jx-1) * xi_2_faces(jx-1) / (x_2(jx) - x_2(jx-1)) * &
                  (dkr(jx, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) + kr_faces(jx-1))
  
    !! q at jx part
    J(jx, jx) = J(jx, jx) - dtflow/dxx_2(jx) * K_faces(jx) * xi_2_faces(jx) / (x_2(jx+1) - x_2(jx)) * &
                  (dkr(jx, jy, jz)*dxx_2(jx+1)/(dxx_2(jx+1) + dxx_2(jx))*(psi(jx+1, jy, jz) - psi(jx, jy, jz)) - kr_faces(jx))
  
    !! q at jx part
    J(jx, jx + 1) = J(jx, jx + 1) - dtflow/dxx_2(jx) * K_faces(jx) * xi_2_faces(jx) / (x_2(jx+1) - x_2(jx)) * &
                  (dkr(jx+1, jy, jz)*dxx_2(jx)/(dxx_2(jx+1) + dxx_2(jx))*(psi(jx+1, jy, jz) - psi(jx, jy, jz)) + kr_faces(jx))
  
    ! add gravitational flux part
    !! q at jx-1 part
    J(jx, jx - 1) = J(jx, jx - 1) + dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces(jx-1)*MERGE(0.0d0, dkr(jx-1, jy, jz) * xi_2(jx-1, jy, jz), head(jx, jy, jz) - head(jx-1, jy, jz) >= 0.0d0)
  
    !! q at jx-1 part
    J(jx, jx) = J(jx, jx) + dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces(jx-1)*MERGE(dkr(jx, jy, jz) * xi_2(jx, jy, jz), 0.0d0, head(jx, jy, jz) - head(jx-1, jy, jz) >= 0.0d0)
  
    !! q at jx part
    J(jx, jx) = J(jx, jx) - dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces(jx)*MERGE(0.0d0, dkr(jx, jy, jz) * xi_2(jx, jy, jz), head(jx+1, jy, jz) - head(jx, jy, jz) >= 0.0d0)
  
    !! q at jx part
    J(jx, jx + 1) = J(jx, jx + 1) - dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces(jx)*MERGE(dkr(jx+1, jy, jz) * xi_2(jx+1, jy, jz), 0.0d0, head(jx+1, jy, jz) - head(jx, jy, jz) >= 0.0d0)
  
    ! add temporal derivative part
    J(jx, jx) = J(jx, jx) + dtheta(jx, jy, jz)
  
  END DO inner

  ! boundary condition at the inlet (begin boundary condition)
  SELECT CASE (x_begin_BC_type)
  CASE ('constant_dirichlet', 'variable_dirichlet')
    J(0, 0) = dxx_2(1)/(dxx_2(0) + dxx_2(1))
    J(0, 1) = dxx_2(0)/(dxx_2(0) + dxx_2(1))
  CASE ('constant_neumann', 'variable_neumann')
    J(0, 0) = 1.0d0/(x_2(1) - x_2(0))
    J(0, 1) = -1.0d0/(x_2(1) - x_2(0))
  
  CASE ('constant_flux', 'variable_flux')
    ! add diffusive flux part
    jx = 1
    J(0, 0) = J(0, 0) -  K_faces(jx-1) * xi_2_faces(jx-1) / (x_2(jx) - x_2(jx-1)) * &
                  (dkr(jx-1, jy, jz)*dxx_2(jx)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) - kr_faces(jx-1))
  
    J(0, 1) = J(0, 1) - K_faces(jx-1) * xi_2_faces(jx-1) / (x_2(jx) - x_2(jx-1)) * &
                  (dkr(jx, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) + kr_faces(jx-1))
  
    ! add gravitational flux part
    jx = 1
    J(0, 0) = J(0, 0) - SignGravity*COSD(x_angle)*K_faces(jx-1)*MERGE(0.0d0, dkr(jx-1, jy, jz) * xi_2(jx-1, jy, jz), head(jx, jy, jz) - head(jx-1, jy, jz) >= 0.0d0)
  
    J(0, 1) = J(0, 1) - SignGravity*COSD(x_angle)*K_faces(jx-1)*MERGE(dkr(jx, jy, jz) * xi_2(jx, jy, jz), 0.0d0, head(jx, jy, jz) - head(jx-1, jy, jz) >= 0.0d0)
  
  CASE ('constant_atomosphere', 'variable_atomosphere')
    psi_b = (psi(0, jy, jz)*dxx_2(1) + psi(1, jy, jz)*dxx_2(0))/(dxx_2(0) + dxx_2(1))
    IF (psi_b < psi_0) THEN
      ! the boundary water potential is below psi_0, so switch to Dirichlet BC
      J(0, 0) = dxx_2(1)/(dxx_2(0) + dxx_2(1))
      J(0, 1) = dxx_2(0)/(dxx_2(0) + dxx_2(1))
    ELSE
      ! add diffusive flux part
      jx = 1
      J(0, 0) = J(0, 0) -  K_faces(jx-1) * xi_2_faces(jx-1) / (x_2(jx) - x_2(jx-1)) * &
                    (dkr(jx-1, jy, jz)*dxx_2(jx)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) - kr_faces(jx-1))
  
      J(0, 1) = J(0, 1) - K_faces(jx-1) * xi_2_faces(jx-1) / (x_2(jx) - x_2(jx-1)) * &
                    (dkr(jx, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) + kr_faces(jx-1))
  
      ! add gravitational flux part
      jx = 1
      J(0, 0) = J(0, 0) - SignGravity*COSD(x_angle)*K_faces(jx-1)*MERGE(0.0d0, dkr(jx-1, jy, jz) * xi_2(jx-1, jy, jz), head(jx, jy, jz) - head(jx-1, jy, jz) >= 0.0d0)
  
      J(0, 1) = J(0, 1) - SignGravity*COSD(x_angle)*K_faces(jx-1)*MERGE(dkr(jx, jy, jz) * xi_2(jx, jy, jz), 0.0d0, head(jx, jy, jz) - head(jx-1, jy, jz) >= 0.0d0)
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
    J(nx+1, nx) = dxx_2(nx+1)/(dxx_2(nx) + dxx_2(nx+1))
    J(nx+1, nx+1) = dxx_2(nx)/(dxx_2(nx) + dxx_2(nx+1))
  
  CASE ('constant_neumann', 'variable_neumann')
    J(nx+1, nx) = -1.0d0 / (x_2(nx+1) - x_2(nx))
    J(nx+1, nx+1) = 1.0d0 / (x_2(nx+1) - x_2(nx))
  
  CASE ('constant_flux', 'variable_flux')
    ! add diffusive flux part
    jx = nx + 1
    J(nx+1, nx) = J(nx+1, nx) - K_faces(jx-1) * xi_2_faces(jx-1) / (x_2(jx) - x_2(jx-1)) * &
                  (dkr(jx-1, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) - kr_faces(jx-1))
  
    J(nx+1, nx+1) = J(nx+1, nx+1) - K_faces(jx-1) * xi_2_faces(jx-1) / (x_2(jx) - x_2(jx-1)) * &
                  (dkr(jx, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) + kr_faces(jx-1))
  
    ! add gravitational flux part
    J(nx+1, nx) = J(nx+1, nx) - SignGravity*COSD(x_angle)*K_faces(jx-1)*MERGE(0.0d0, dkr(jx-1, jy, jz) * xi_2(jx-1, jy, jz), head(jx, jy, jz) - head(jx-1, jy, jz) >= 0.0d0)
  
    J(nx+1, nx+1) = J(nx+1, nx+1) - SignGravity*COSD(x_angle)*K_faces(jx-1)*MERGE(dkr(jx, jy, jz) * xi_2(jx, jy, jz), 0.0d0, head(jx, jy, jz) - head(jx-1, jy, jz) >= 0.0d0)
  
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
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Currently, only regular spatial domain is supported.'
    WRITE(*,*)
    READ(*,*)
    STOP
        
  END IF
  
  
  ! internal cells
  DO i = 1, n_inner_cells
    face_IDs = cell_to_face(i, :)
    face_ID_south = face_IDs(1)
    face_ID_north = face_IDs(2)
    face_ID_west = face_IDs(3)
    face_ID_east = face_IDs(4)
    
    cell_ID_south = face_to_cell(face_ID_south, 1)
    cell_ID_north = face_to_cell(face_ID_north, 2)
    cell_ID_west = face_to_cell(face_ID_west, 1)
    cell_ID_east = face_to_cell(face_ID_east, 2)
    
    jx = cell_to_coordinate(i, 1)
    jy = cell_to_coordinate(i, 2)
    jz = 1
    
    ! add temporal derivative part
    J(i, i) = J(i, i) + dtheta(jx, jy, jz)
    
    ! add diffusive flux part
    
    !! center cell
    !!! due to q at west face
    J(i, i) = J(i, i) + dtflow/dxx_2(jx) * K_faces(face_ID_west) * xi_2_faces(face_ID_west) / (x_2(jx) - x_2(jx-1)) * &
                  (dkr(jx, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) + kr_faces(face_ID_west))
  
    !!! due to q at east face
    J(i, i) = J(i, i) - dtflow/dxx_2(jx) * K_faces(face_ID_east) * xi_2_faces(face_ID_east) / (x_2(jx+1) - x_2(jx)) * &
                  (dkr(jx, jy, jz)*dxx_2(jx+1)/(dxx_2(jx+1) + dxx_2(jx))*(psi(jx+1, jy, jz) - psi(jx, jy, jz)) - kr_faces(face_ID_east))
  
    !!! due to q at south face
    J(i, i) = J(i, i) + dtflow/dyy_2(jy) * K_faces(face_ID_south) * xi_2_faces(face_ID_south) / (y_2(jy) - y_2(jy-1)) * &
                  (dkr(jx, jy, jz)*dyy_2(jy-1)/(dyy_2(jy) + dyy_2(jy-1))*(psi(jx, jy, jz) - psi(jx, jy-1, jz)) + kr_faces(face_ID_south))
    
    !!! due to q at east face
    J(i, i) = J(i, i) - dtflow/dyy_2(jy) * K_faces(face_ID_north) * xi_2_faces(face_ID_north) / (y_2(jy+1) - y_2(jy)) * &
                  (dkr(jx, jy, jz)*dyy_2(jy+1)/(dyy_2(jy+1) + dyy_2(jy))*(psi(jx, jy+1, jz) - psi(jx, jy, jz)) - kr_faces(face_ID_north))
    
    !! west cell due to q at west face
    J(i, cell_ID_west) = J(i, cell_ID_west) + dtflow/dxx_2(jx) * K_faces(face_ID_west) * xi_2_faces(face_ID_west) / (x_2(jx) - x_2(jx-1)) * &
                  (dkr(jx-1, jy, jz)*dxx_2(jx)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) - kr_faces(face_ID_west))
    
    !! east cell due to q at east face
    J(i, cell_ID_east) = J(i, cell_ID_east) - dtflow/dxx_2(jx) * K_faces(face_ID_east) * xi_2_faces(face_ID_east) / (x_2(jx+1) - x_2(jx)) * &
                  (dkr(jx+1, jy, jz)*dxx_2(jx)/(dxx_2(jx+1) + dxx_2(jx))*(psi(jx+1, jy, jz) - psi(jx, jy, jz)) + kr_faces(face_ID_east))
    
    !! south cell due to q at south face
    J(i, cell_ID_south) = J(i, face_ID_south) + dtflow/dyy_2(jy) * K_faces(face_ID_west) * xi_2_faces(face_ID_south) / (y_2(jy) - y_2(jy-1)) * &
                  (dkr(jx, jy-1, jz)*dyy_2(jy)/(dyy_2(jy) + dyy_2(jy-1))*(psi(jx, jy, jz) - psi(jx, jy-1, jz)) - kr_faces(face_ID_south))
    
    !! north cell due to q at north face
    J(i, cell_ID_north) = J(i, cell_ID_north) - dtflow/dyy_2(jy) * K_faces(cell_ID_north) * xi_2_faces(cell_ID_north) / (y_2(jy+1) - y_2(jy)) * &
                  (dkr(jx, jy+1, jz)*dyy_2(jy)/(dyy_2(jy+1) + dyy_2(jy))*(psi(jx, jy+1, jz) - psi(jx, jy, jz)) + kr_faces(cell_ID_north))
  
    ! add gravitational flux part
    
    !! center cell
    !!! due to q at west face
    J(i, i) = J(i, i) + dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces(face_ID_west)*MERGE(dkr(jx, jy, jz) * xi_2(jx, jy, jz), 0.0d0, head(jx, jy, jz) - head(jx-1, jy, jz) >= 0.0d0)
  
    !!! due to q at east face
    J(i, i) = J(i, i) - dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces(face_ID_east)*MERGE(0.0d0, dkr(jx, jy, jz) * xi_2(jx, jy, jz), head(jx+1, jy, jz) - head(jx, jy, jz) >= 0.0d0)
  
    !!! due to q at south face
    J(i, i) = J(i, i) + dtflow/dyy_2(jy) * SignGravity*COSD(y_angle)*K_faces(face_ID_south)*MERGE(dkr(jx, jy, jz) * xi_2(jx, jy, jz), 0.0d0, head(jx, jy, jz) - head(jx, jy-1, jz) >= 0.0d0)
  
    !!! due to q at north face
    J(i, i) = J(i, i) - dtflow/dyy_2(jy) * SignGravity*COSD(y_angle)*K_faces(face_ID_north)*MERGE(0.0d0, dkr(jx, jy, jz) * xi_2(jx, jy, jz), head(jx, jy+1, jz) - head(jx, jy, jz) >= 0.0d0)
  
    !! west cell due to q at west face
    J(i, cell_ID_west) = J(i, cell_ID_west) + dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces(face_ID_west)*MERGE(0.0d0, dkr(jx-1, jy, jz) * xi_2(jx-1, jy, jz), head(jx, jy, jz) - head(jx-1, jy, jz) >= 0.0d0)
  
    
    !! east cell due to q at east face
    J(i, cell_ID_east) = J(i, cell_ID_east) - dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces(face_ID_east)*MERGE(dkr(jx+1, jy, jz) * xi_2(jx+1, jy, jz), 0.0d0, head(jx+1, jy, jz) - head(jx, jy, jz) >= 0.0d0)
  
    !! south cell due to q at south face
    J(i, cell_ID_south) = J(i, cell_ID_south) + dtflow/dyy_2(jy) * SignGravity*COSD(y_angle)*K_faces(face_ID_south)*MERGE(0.0d0, dkr(jx, jy-1, jz) * xi_2(jx, jy-1, jz), head(jx, jy, jz) - head(jx, jy-1, jz) >= 0.0d0)
    
    !! north cell due to q at north face
    J(i, cell_ID_north) = J(i, cell_ID_north) - dtflow/dyy_2(jy) * SignGravity*COSD(y_angle)*K_faces(face_ID_north)*MERGE(dkr(jx, jy+1, jz) * xi_2(jx, jy+1, jz), 0.0d0, head(jx, jy+1, jz) - head(jx, jy, jz) >= 0.0d0)
  
  END DO
  
  
  ! evaluate residual for boundary conditions
  DO i = 1, n_bfaces
    face_ID = bface_to_face(i)
  
    cell_IDs = face_to_cell(face_ID, :)
    cell_ID_inner = minval(cell_IDs)
    cell_ID_outer = maxval(cell_IDs)
    coordinate_cell_inner = cell_to_coordinate(cell_ID_inner, :)
    coordinate_cell_outer = cell_to_coordinate(cell_ID_outer, :)
    
    jx_1 = coordinate_cell_inner(1)
    jx_2 = coordinate_cell_outer(1)
    jy_1 = coordinate_cell_inner(2)
    jy_2 = coordinate_cell_outer(2)
    
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
      
      J(cell_ID_outer, cell_ID_outer) = dx_1/(dx_1 + dx_2)
      J(cell_ID_outer, cell_ID_inner) = dx_1/(dx_1 + dx_2)  
    
    CASE (2) ! Neumann
      IF (jx_1 == jx_2) THEN
        ! faces paralell to x-axis
        dx_1 = y_2(jy_2) - y_2(jy_1) 
      ELSE
        ! faces paralell to y-axis
        dx_1 = x_2(jx_2) - x_2(jx_1)
      END IF
      J(cell_ID_outer, cell_ID_outer) = 1.0d0/dx_1
      J(cell_ID_outer, cell_ID_inner) = -1.0d0/dx_1
        
    CASE (3) ! flux
      IF (jx_1 == jx_2) THEN
        ! faces paralell to x-axis
        dx_1 = y_2(jy_2) - y_2(jy_1)
        dx_2 = dyy_2(jy_2) + dyy_2(jy_1)
        gravity_effect = SignGravity*COSD(y_angle)
      ELSE
        ! faces paralell to y-axis
        dx_1 = x_2(jx_2) - x_2(jx_1)
        dx_2 = dxx_2(jx_2) + dxx_2(jx_1)
        gravity_effect = SignGravity*COSD(x_angle)
      END IF
      ! add diffusive flux part
  
      J(cell_ID_outer, cell_ID_outer) = J(cell_ID_outer, cell_ID_outer) - K_faces(face_ID) * xi_2_faces(face_ID) / dx_1 * &
                    (dkr(jx_2, jy_2, jz)*0.5d0*dx_2/dx_2*(psi(jx_2, jy_2, jz) - psi(jx_1, jy_1, jz)) + kr_faces(face_ID))
      
      J(cell_ID_outer, cell_ID_inner) = J(cell_ID_outer, cell_ID_inner) - K_faces(face_ID) * xi_2_faces(face_ID) / dx_1 * &
                    (dkr(jx_1, jy_1, jz)*0.5d0*dx_2/dx_2*(psi(jx_2, jy_2, jz) - psi(jx_1, jy_1, jz)) - kr_faces(face_ID))
  
      ! add gravitational flux part
      J(cell_ID_outer, cell_ID_inner) = J(cell_ID_outer, cell_ID_inner) - SignGravity*COSD(x_angle)*K_faces(face_ID)*MERGE(0.0d0, dkr(jx_1, jy_1, jz) * xi_2(jx_1, jy_1, jz), head(jx_2, jy_2, jz) - head(jx_1, jy_1, jz) >= 0.0d0)
  
      J(cell_ID_outer, cell_ID_outer) = J(cell_ID_outer, cell_ID_outer) - SignGravity*COSD(x_angle)*K_faces(face_ID)*MERGE(dkr(jx_2, jy_2, jz) * xi_2(jx_2, jy_2, jz), 0.0d0, head(jx_2, jy_2, jz) - head(jx_1, jy_1, jz) >= 0.0d0)
  
    
    
    
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

END SUBROUTINE Jacobian_Richards