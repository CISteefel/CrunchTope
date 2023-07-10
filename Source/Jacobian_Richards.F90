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
REAL(DP), INTENT(OUT), DIMENSION(nx, nx)                   :: J ! Jacobian matrix

INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

REAL(DP), INTENT(IN)                                       :: dtflow

! variables not declared in CrunchTope
REAL(DP)                                                   :: psi_lb
REAL(DP)                                                   :: head_lb
REAL(DP)                                                   :: theta_lb
REAL(DP)                                                   :: dtheta_lb
REAL(DP)                                                   :: kr_lb
REAL(DP)                                                   :: dkr_lb
REAL(DP)                                                   :: psi_ub
REAL(DP)                                                   :: head_ub
REAL(DP)                                                   :: theta_ub
REAL(DP)                                                   :: dtheta_ub
REAL(DP)                                                   :: kr_ub
REAL(DP)                                                   :: dkr_ub
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: xi ! physical constant

xi = rho_water2*g/mu_water

jy = 1
jz = 1
DO jx = 1, nx
  CALL vanGenuchten_model(psi(jx, jy, jz), theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), &
                          theta(jx, jy, jz), kr(jx, jy, jz), dtheta(jx, jy, jz), dkr(jx, jy, jz))
  head(jx, jy, jz) = psi(jx, jy, jz) + x(jx)
END DO

! evaluate Jacobian matrix
J = 0.0 ! initialize Jacobian matrix

! interior cells
inner: DO jx = 2, nx - 1
  J(jx, jx - 1) = 1.0d0/dxx(jx) * xi(jx, jy, jz) * K_faces_x(jx-1, jy, jz) / (x(jx) - x(jx - 1))* &
              MERGE(-kr(jx, jy, jz), dkr(jx-1, jy, jz)*(head(jx, jy, jz) - head(jx-1, jy, jz)) - kr(jx-1, jy, jz), head(jx, jy, jz) - head(jx-1, jy, jz) >= 0)
  J(jx, jx) = dtheta(jx, jy, jz)/dtflow + &
          1.0d0/dxx(jx) * xi(jx, jy, jz) * K_faces_x(jx-1, jy, jz) / (x(jx) - x(jx - 1)) &
          *MERGE(dkr(jx, jy, jz)*(head(jx, jy, jz) - head(jx-1, jy, jz)) + kr(jx, jy, jz), kr(jx-1, jy, jz), head(jx, jy, jz) - head(jx-1, jy, jz) >= 0) &    
          - 1.0d0/dxx(jx) * xi(jx, jy, jz) * K_faces_x(jx, jy, jz) / (x(jx + 1) - x(jx)) &
          * MERGE(-kr(jx+1, jy, jz), dkr(jx, jy, jz)*(head(jx+1, jy, jz) - head(jx, jy, jz)) - kr(jx, jy, jz), head(jx+1, jy, jz) - head(jx, jy, jz) >= 0)
  J(jx, jx + 1) = - 1.0d0/dxx(jx) * xi(jx, jy, jz) * K_faces_x(jx, jy, jz) / (x(jx + 1) - x(jx)) &
              * MERGE(dkr(jx+1, jy, jz)*(head(jx+1, jy, jz) - head(jx, jy, jz)) + kr(jx+1, jy, jz), kr(jx, jy, jz), head(jx+1, jy, jz) - head(jx, jy, jz) >= 0)
END DO inner

! lower boundary (the first cell)
jx = 1
J(1, 1) = dtheta(1, jy, jz)/dtflow &
        - 1.0d0/dxx(1) * xi(jx, jy, jz)*K_faces_x(1, jy, jz)/(x(2) - x(1))*MERGE(-kr(2, jy, jz), dkr(1, jy, jz)*(head(2, jy, jz) - head(1, jy, jz)) - kr(1, jy, jz), head(2, jy, jz) - head(1, jy, jz) >= 0)

J(1, 2) = - 1.0d0/dxx(1) * xi(jx, jy, jz) * K_faces_x(1, jy, jz) / (x(2) - x(1)) * MERGE(dkr(2, jy, jz)*(head(2, jy, jz) - head(1, jy, jz)) + kr(2, jy, jz), kr(1, jy, jz), head(2, jy, jz) - head(1, jy, jz) >= 0)

! the derivative of the lower boundary flux with respect to the first cell psi depends on boundary condition
SELECT CASE (lower_BC_type)
CASE ('constant_dirichlet', 'variable_dirichlet')
  head_lb = value_lower_BC + (x(jx) - 0.5d0 * dxx(jx))
  CALL vanGenuchten_model_kr(value_lower_BC, theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), kr_lb)
  J(1, 1) = J(1, 1) + 1.0d0/dxx(1) * xi(jx, jy, jz)*K_faces_x(0, jy, jz)/(dxx(1)*0.5d0)*MERGE(dkr(1, jy, jz)*(head(1, jy, jz) - head_lb) + kr(1, jy, jz), kr_lb, head(1, jy, jz) - head_lb >= 0)
CASE ('constant_neumann', 'variable_neumann')
  psi_lb = psi(1, jy, jz) - value_lower_BC*(0.5d0 * dxx(jx))
  head_lb = psi_lb + (x(jx) - 0.5d0 * dxx(jx))
  CALL vanGenuchten_model(psi_lb, theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), &
                          theta_lb, kr_lb, dtheta_lb, dkr_lb)
  ! head(1, jy, jz) - head_lb does not depend on psi(1, jy, jz), so the derivative is simple
  ! derivative of k_lb with respect to psi_1 is dkr_lb * d psi_lb/d psi_1 = dkr_lb
  J(1, 1) = J(1, 1) + 1.0d0/dxx(1) * xi(jx, jy, jz)*K_faces_x(0, jy, jz)/(dxx(1)*0.5d0)*MERGE(dkr(1, jy, jz)*(head(1, jy, jz) - head_lb), dkr_lb*(head(1, jy, jz) - head_lb), head(1, jy, jz) - head_lb >= 0)
  CONTINUE
CASE ('constant_flux', 'variable_flux')
  ! the derivative of the lower boundary flux with respect to the first cell psi is zero, so do nothing
  CONTINUE
CASE DEFAULT
  WRITE(*,*)
  WRITE(*,*) ' The boundary condition type ', lower_BC_type, ' is not supported. '
  WRITE(*,*)
  READ(*,*)
  STOP
  
END SELECT


! upper boundary (the last cell)
jx = nx
J(nx, nx - 1) = 1.0d0/dxx(nx) * xi(jx, jy, jz) * K_faces_x(nx-1, jy, jz) / (x(nx) - x(nx-1)) * &
                MERGE(-kr(nx, jy, jz), dkr(nx-1, jy, jz)*(head(nx, jy, jz) - head(nx-1, jy, jz)) - kr(nx-1, jy, jz), head(nx, jy, jz) - head(nx-1, jy, jz) >= 0)

J(nx, nx) = dtheta(nx, jy, jz)/dtflow + 1.0d0/dxx(nx) * xi(jx, jy, jz) * K_faces_x(nx-1, jy, jz) / (x(nx) - x(nx-1)) &
            *MERGE(dkr(nx, jy, jz)*(head(nx, jy, jz) - head(nx-1, jy, jz)) + kr(nx, jy, jz), kr(nx-1, jy, jz), head(nx, jy, jz) - head(nx-1, jy, jz) >= 0)

! the derivative of the upper boundary flux with respect to the lass cell psi depends on boundary condition
SELECT CASE (upper_BC_type)
CASE ('constant_dirichlet', 'variable_dirichlet')
  head_ub = value_upper_BC + (x(jx) + 0.5d0 * dxx(jx))
  CALL vanGenuchten_model_kr(value_upper_BC, theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), kr_ub)
  J(nx, nx) = J(nx, nx) - 1.0d0/dxx(nx) * xi(jx, jy, jz)*K_faces_x(nx, jy, jz)/(dxx(nx)*0.5d0)*MERGE(-kr_ub, dkr(nx, jy, jz)*(head_ub - head(nx, jy, jz)), head_ub - head(nx, jy, jz) >= 0)
  
CASE ('constant_neumann', 'variable_neumann')
  CONTINUE
CASE ('constant_flux', 'variable_flux')
  ! the derivative of the upper boundary flux with respect to the last cell psi is zero, so do nothing
  CONTINUE
  
CASE ('environmental_forcing')
  ! the derivative of the upper boundary flux with respect to the last cell psi is zero, so do nothing
  CONTINUE
CASE DEFAULT
  WRITE(*,*)
  WRITE(*,*) ' The boundary condition type ', upper_BC_type, ' is not supported. '
  WRITE(*,*)
  READ(*,*)
  STOP
  
END SELECT

END SUBROUTINE Jacobian_Richards