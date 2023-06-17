SUBROUTINE flux_Richards_steady(nx, ny, nz)
! This subroutine computes water flux based on the Buckingham Darcy's law
! The computed water flux is positive upward with respect to the x-direction
! This subroutine is used for the steady-state Richards equation solver (the only difference from the flux_Richards.F90 is boundary condition)
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

REAL(DP)                                                   :: psi_lb
REAL(DP)                                                   :: head_lb
REAL(DP)                                                   :: kr_lb
REAL(DP)                                                   :: psi_ub
REAL(DP)                                                   :: head_ub
REAL(DP)                                                   :: kr_ub
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: xi ! physical constant

xi = rho_water2*g/mu_water

! apply van Genuchten model to all grid cells
jy = 1
jz = 1
DO jx = 1, nx
  CALL vanGenuchten_model(psi(jx, jy, jz), theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), &
                 theta(jx, jy, jz), kr(jx, jy, jz), dtheta(jx, jy, jz), dkr(jx, jy, jz))
  head(jx, jy, jz) = psi(jx, jy, jz) + x(jx)
END DO

! internal faces
DO jx = 1, nx - 1
  qx(jx, jy, jz) = -xi(jx,jy,jz)*K_faces_x(jx, jy, jz)/(x(jx+1) - x(jx))*(kr(jx+1, jy, jz)*MAX(head(jx+1, jy, jz) - head(jx, jy, jz), 0.0d0) &
                      + kr(jx, jy, jz)*MIN(head(jx+1, jy, jz) - head(jx, jy, jz), 0.0d0))
END DO

! lower boundary face
jx = 1
SELECT CASE (lower_BC_type_steady)
CASE ('constant_dirichlet', 'variable_dirichlet')
  head_lb = value_lower_BC_steady + (x(jx) - 0.5d0 * dxx(jx))
  CALL vanGenuchten_model_kr(value_lower_BC_steady, theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), kr_lb)
  qx(0, jy, jz) = -xi(jx,jy,jz)*K_faces_x(0, jy, jz)/(dxx(1)*0.5d0)*((kr(1, jy, jz))*MAX(head(1, jy, jz) - head_lb, 0.0d0) + kr_lb*MIN(head(1, jy, jz) - head_lb, 0.0d0))
  
CASE ('constant_neumann', 'variable_neumann')
  ! compute the water potential at the lower boundary from the given gradient of the water potential
  ! this is based on the one-sided finite difference approximation at the boundary
  psi_lb = psi(1, jy, jz) - value_lower_BC_steady*(0.5d0 * dxx(jx))
  head_lb = psi_lb + (x(jx) - 0.5d0 * dxx(jx))
  CALL vanGenuchten_model_kr(psi_lb, theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), kr_lb)
  qx(0, jy, jz) = -xi(jx,jy,jz)*K_faces_x(0, jy, jz)/(dxx(1)*0.5d0)*((kr(1, jy, jz))*MAX(head(1, jy, jz) - head_lb, 0.0d0) + kr_lb*MIN(head(1, jy, jz) - head_lb, 0.0d0))
  
CASE ('constant_flux', 'variable_flux')
  qx(0, jy, jz) = value_lower_BC_steady
  
CASE DEFAULT
  WRITE(*,*)
  WRITE(*,*) ' The boundary condition type ', lower_BC_type_steady, ' is not supported. '
  WRITE(*,*)
  READ(*,*)
  STOP
  
END SELECT

! upper boundary face
jx = nx
SELECT CASE (upper_BC_type_steady)
CASE ('constant_dirichlet', 'variable_dirichlet')
  
  head_ub = value_upper_BC_steady + (x(jx) + 0.5d0 * dxx(jx))
  CALL vanGenuchten_model_kr(value_upper_BC_steady, theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), kr_ub)
  qx(nx, jy, jz) = -xi(jx,jy,jz)*K_faces_x(nx, jy, jz)/(dxx(nx)*0.5d0)*((kr_ub)*MAX(head_ub - head(nx, jy, jz), 0.0d0) + kr(nx, jy, jz)*MIN(head_ub - head(nx, jy, jz), 0.0d0))
  
CASE ('constant_neumann', 'variable_neumann')
  
  psi_ub = psi(1, jy, jz) + value_upper_BC_steady*(0.5d0 * dxx(jx))
  head_ub = psi_ub + (x(jx) + 0.5d0 * dxx(jx))
  qx(nx, jy, jz) = -xi(jx,jy,jz)*K_faces_x(nx, jy, jz)/(dxx(nx)*0.5d0)*((kr_ub)*MAX(head_ub - head(nx, jy, jz), 0.0d0) + kr(nx, jy, jz)*MIN(head_ub - head(nx, jy, jz), 0.0d0))
  
CASE ('constant_flux', 'variable_flux')
  qx(nx, jy, jz) = value_upper_BC_steady
  
CASE DEFAULT
  WRITE(*,*)
  WRITE(*,*) ' The boundary condition type ', upper_BC_type_steady, ' is not supported. '
  WRITE(*,*)
  READ(*,*)
  STOP
  
END SELECT


END SUBROUTINE flux_Richards_steady