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
INTEGER(I4B), INTENT(IN)                                   :: low_bound_F
INTEGER(I4B), INTENT(IN)                                   :: high_bound_F

REAL(DP), INTENT(OUT), DIMENSION(low_bound_F:high_bound_F, low_bound_F:high_bound_F)                   :: J ! Jacobian matrix
INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

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
    jx = cell_to_coordinate(i, 1)
    jy = cell_to_coordinate(i, 2)
    jz = 1
    
    q_y_1 = q_Richards(face_IDs(1))
    q_y_2 = q_Richards(face_IDs(2))
    q_x_1 = q_Richards(face_IDs(3))
    q_x_2 = q_Richards(face_IDs(4))
    divergence = (q_y_2 - q_y_1) / dyy_2(jy) + (q_x_2 - q_x_1) / dxx_2(jx)
    F_residual(i) = theta(jx, jy, jz) - theta_prev(jx, jy, jz) + divergence*dtflow
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