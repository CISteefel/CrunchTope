SUBROUTINE flux_Richards(nx, ny, nz)
! This subroutine computes water flux based on the Buckingham Darcy's law
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

REAL(DP)                                                   :: psi_grad ! gradient of water potential
REAL(DP)                                                   :: q_diff ! diffusion flow
REAL(DP)                                                   :: q_grav ! gravitational flow

xi_2 = rho_water_2*g/mu_water


! apply van Genuchten model to all grid cells
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

! flux at faces (only gravitational flow is upwinded)
DO jx = 0, nx
  psi_grad = (psi(jx+1, jy, jz) - psi(jx, jy, jz))/(x_2(jx+1) - x_2(jx))
  q_diff = -xi_2_faces(jx, jy, jz)*K_faces(jx)*kr_faces(jx, jy, jz)*psi_grad
  q_grav = -SignGravity*COSD(x_angle)*K_faces(jx)*MERGE(kr(jx+1, jy, jz) * xi_2(jx+1, jy, jz), kr(jx, jy, jz) * xi_2(jx, jy, jz), head(jx+1, jy, jz) - head(jx, jy, jz) >= 0.0d0)
  qx(jx, jy, jz) = q_diff + q_grav
END DO

END SUBROUTINE flux_Richards