SUBROUTINE Jacobian_Richards(nx, ny, nz, dtflow, J)
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
REAL(DP), INTENT(OUT), DIMENSION(0:nx+1, 0:nx+1)                   :: J ! Jacobian matrix
INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

REAL(DP), INTENT(IN)                                       :: dtflow
REAL(DP)                                                   :: psi_b ! water potential at a boundary

xi_2 = rho_water_2*g/mu_water

jy = 1
jz = 1
DO jx = 0, nx+1
  CALL vanGenuchten_model(psi(jx, jy, jz), theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), &
                          theta(jx, jy, jz), kr(jx, jy, jz), dtheta(jx, jy, jz), dkr(jx, jy, jz))
  head(jx, jy, jz) = psi(jx, jy, jz) + SignGravity * COSD(x_angle) * x_2(jx)
END DO

! compute xi and kr at faces
DO jx = 0, nx
  xi_2_faces(jx, jy, jz) = (xi_2(jx, jy, jz)*dxx_2(jx+1) + xi_2(jx+1, jy, jz)*dxx_2(jx))/(dxx_2(jx+1) + dxx_2(jx))
  kr_faces(jx, jy, jz) = (kr(jx, jy, jz)*dxx_2(jx+1) + kr(jx+1, jy, jz)*dxx_2(jx))/(dxx_2(jx+1) + dxx_2(jx))
END DO


! evaluate Jacobian matrix
J = 0.0 ! initialize Jacobian matrix
  
! interior cells
inner: DO jx = 1, nx  
  ! add diffusive flux part
  !! q at jx-1 part
  J(jx, jx - 1) = J(jx, jx - 1) + dtflow/dxx_2(jx) * K_faces_x(jx-1, jy, jz) * xi_2_faces(jx-1, jy, jz) / (x_2(jx) - x_2(jx-1)) * &
                (dkr(jx-1, jy, jz)*dxx_2(jx)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) - kr_faces(jx-1, jy, jz))
  
  !! q at jx-1 part
  J(jx, jx) = J(jx, jx) + dtflow/dxx_2(jx) * K_faces_x(jx-1, jy, jz) * xi_2_faces(jx-1, jy, jz) / (x_2(jx) - x_2(jx-1)) * &
                (dkr(jx, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) + kr_faces(jx-1, jy, jz))
  
  !! q at jx part
  J(jx, jx) = J(jx, jx) - dtflow/dxx_2(jx) * K_faces_x(jx, jy, jz) * xi_2_faces(jx, jy, jz) / (x_2(jx+1) - x_2(jx)) * &
                (dkr(jx, jy, jz)*dxx_2(jx+1)/(dxx_2(jx+1) + dxx_2(jx))*(psi(jx+1, jy, jz) - psi(jx, jy, jz)) - kr_faces(jx, jy, jz))
  
  !! q at jx part
  J(jx, jx + 1) = J(jx, jx + 1) - dtflow/dxx_2(jx) * K_faces_x(jx, jy, jz) * xi_2_faces(jx, jy, jz) / (x_2(jx+1) - x_2(jx)) * &
                (dkr(jx+1, jy, jz)*dxx_2(jx)/(dxx_2(jx+1) + dxx_2(jx))*(psi(jx+1, jy, jz) - psi(jx, jy, jz)) + kr_faces(jx, jy, jz))
  
  ! add gravitational flux part
  !! q at jx-1 part
  J(jx, jx - 1) = J(jx, jx - 1) + dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces_x(jx-1, jy, jz)*MERGE(0.0d0, dkr(jx-1, jy, jz) * xi_2(jx-1, jy, jz), SignGravity * COSD(x_angle) * (x_2(jx) - x_2(jx-1)) >= 0.0d0)
  
  !! q at jx-1 part
  J(jx, jx) = J(jx, jx) + dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces_x(jx-1, jy, jz)*MERGE(dkr(jx, jy, jz) * xi_2(jx, jy, jz), 0.0d0, SignGravity * COSD(x_angle) * (x_2(jx) - x_2(jx-1)) >= 0.0d0)
  
  !! q at jx part
  J(jx, jx) = J(jx, jx) - dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces_x(jx, jy, jz)*MERGE(0.0d0, dkr(jx, jy, jz) * xi_2(jx, jy, jz), SignGravity * COSD(x_angle) * (x_2(jx+1) - x_2(jx)) >= 0.0d0)
  
  !! q at jx part
  J(jx, jx + 1) = J(jx, jx + 1) - dtflow/dxx_2(jx) * SignGravity*COSD(x_angle)*K_faces_x(jx, jy, jz)*MERGE(dkr(jx+1, jy, jz) * xi_2(jx+1, jy, jz), 0.0d0, SignGravity * COSD(x_angle) * (x_2(jx+1) - x_2(jx)) >= 0.0d0)
  
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
  J(0, 0) = J(0, 0) -  K_faces_x(jx-1, jy, jz) * xi_2_faces(jx-1, jy, jz) / (x_2(jx) - x_2(jx-1)) * &
                (dkr(jx-1, jy, jz)*dxx_2(jx)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) - kr_faces(jx-1, jy, jz))
  
  J(0, 1) = J(0, 1) - K_faces_x(jx-1, jy, jz) * xi_2_faces(jx-1, jy, jz) / (x_2(jx) - x_2(jx-1)) * &
                (dkr(jx, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) + kr_faces(jx-1, jy, jz))
  
  ! add gravitational flux part
  jx = 1
  J(0, 0) = J(0, 0) - SignGravity*COSD(x_angle)*K_faces_x(jx-1, jy, jz)*MERGE(0.0d0, dkr(jx-1, jy, jz) * xi_2(jx-1, jy, jz), SignGravity * COSD(x_angle) * (x_2(jx) - x_2(jx-1)) >= 0.0d0)
  
  J(0, 1) = J(0, 1) - SignGravity*COSD(x_angle)*K_faces_x(jx-1, jy, jz)*MERGE(dkr(jx, jy, jz) * xi_2(jx, jy, jz), 0.0d0, SignGravity * COSD(x_angle) * (x_2(jx) - x_2(jx-1)) >= 0.0d0)
  
CASE ('constant_atomosphere', 'variable_atomosphere')
  psi_b = (psi(0, jy, jz)*dxx_2(1) + psi(1, jy, jz)*dxx_2(0))/(dxx_2(0) + dxx_2(1))
  IF (psi_b < psi_0) THEN
    ! the boundary water potential is below psi_0, so switch to Dirichlet BC
    J(0, 0) = dxx_2(1)/(dxx_2(0) + dxx_2(1))
    J(0, 1) = dxx_2(0)/(dxx_2(0) + dxx_2(1))
  ELSE
    ! add diffusive flux part
    jx = 1
    J(0, 0) = J(0, 0) -  K_faces_x(jx-1, jy, jz) * xi_2_faces(jx-1, jy, jz) / (x_2(jx) - x_2(jx-1)) * &
                  (dkr(jx-1, jy, jz)*dxx_2(jx)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) - kr_faces(jx-1, jy, jz))
  
    J(0, 1) = J(0, 1) - K_faces_x(jx-1, jy, jz) * xi_2_faces(jx-1, jy, jz) / (x_2(jx) - x_2(jx-1)) * &
                  (dkr(jx, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) + kr_faces(jx-1, jy, jz))
  
    ! add gravitational flux part
    jx = 1
    J(0, 0) = J(0, 0) - SignGravity*COSD(x_angle)*K_faces_x(jx-1, jy, jz)*MERGE(0.0d0, dkr(jx-1, jy, jz) * xi_2(jx-1, jy, jz), SignGravity * COSD(x_angle) * (x_2(jx) - x_2(jx-1)) >= 0.0d0)
  
    J(0, 1) = J(0, 1) - SignGravity*COSD(x_angle)*K_faces_x(jx-1, jy, jz)*MERGE(dkr(jx, jy, jz) * xi_2(jx, jy, jz), 0.0d0, SignGravity * COSD(x_angle) * (x_2(jx) - x_2(jx-1)) >= 0.0d0)
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
  J(nx+1, nx) = J(nx+1, nx) - K_faces_x(jx-1, jy, jz) * xi_2_faces(jx-1, jy, jz) / (x_2(jx) - x_2(jx-1)) * &
                (dkr(jx-1, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) - kr_faces(jx-1, jy, jz))
  
  J(nx+1, nx+1) = J(nx+1, nx+1) - K_faces_x(jx-1, jy, jz) * xi_2_faces(jx-1, jy, jz) / (x_2(jx) - x_2(jx-1)) * &
                (dkr(jx, jy, jz)*dxx_2(jx-1)/(dxx_2(jx) + dxx_2(jx-1))*(psi(jx, jy, jz) - psi(jx-1, jy, jz)) + kr_faces(jx-1, jy, jz))
  
  ! add gravitational flux part
  J(nx+1, nx) = J(nx+1, nx) - SignGravity*COSD(x_angle)*K_faces_x(jx-1, jy, jz)*MERGE(0.0d0, dkr(jx-1, jy, jz) * xi_2(jx-1, jy, jz), SignGravity * COSD(x_angle) * (x_2(jx) - x_2(jx-1)) >= 0.0d0)
  
  J(nx+1, nx+1) = J(nx+1, nx+1) - SignGravity*COSD(x_angle)*K_faces_x(jx-1, jy, jz)*MERGE(dkr(jx, jy, jz) * xi_2(jx, jy, jz), 0.0d0, SignGravity * COSD(x_angle) * (x_2(jx) - x_2(jx-1)) >= 0.0d0)
  
CASE DEFAULT
  WRITE(*,*)
  WRITE(*,*) ' The boundary condition type ', x_end_BC_type, ' is not supported. '
  WRITE(*,*)
  READ(*,*)
  STOP
  
END SELECT

END SUBROUTINE Jacobian_Richards