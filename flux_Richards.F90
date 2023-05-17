SUBROUTINE flux_Richards(nx, ny, nz, psi_lb_value, qx_ub_value)
! subroutine to compute water flux based on the Buckingham Darcy's law
USE crunchtype
USE io
USE params
USE runtime
USE concentration
USE medium
USE flow
USE transport

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                          :: nx
INTEGER(I4B), INTENT(IN)                                          :: ny
INTEGER(I4B), INTENT(IN)                                          :: nz

INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

REAL(DP), INTENT(IN)                                                   :: psi_lb_value
REAL(DP), INTENT(IN)                                                   :: qx_ub_value
REAL(DP)                                                   :: head_lb
REAL(DP)                                                   :: kr_lb

!**************************************************
! physical parameters for Richards solver added by Toshiyuki Bandai, 2023, May
REAL(DP), PARAMETER :: mu = 1.0016d0 * 1.0E-3 * (86400.0d0 * 365.0d0) ! dynamics viscosity of water [Pa year] at 20 degC
REAL(DP), PARAMETER :: rho = 0.99823d0 * 1.0E3 ! density of water [kg m-3] at 20 degC
REAL(DP), PARAMETER :: g = 9.80665d0 * (86400.0d0 * 365.0d0) ** 2 ! gravitational acceleration [m year-2]
REAL(DP), PARAMETER :: xi = rho*g/mu ! constant used to solve the Richards equation
! End of edits by Toshiyuki Bandai, 2023, May
!**************************************************

jy = 1
jz = 1
DO jx = 1, nx
  CALL vanGenuchten_model(psi(jx, jy, jz), theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), psi_s(jx, jy, jz), &
                 theta(jx, jy, jz), kr(jx, jy, jz), dtheta(jx, jy, jz), dkr(jx, jy, jz))
  head(jx, jy, jz) = psi(jx, jy, jz) + x(jx)
END DO

jx = 1
jy = 1
jz = 1
head_lb = psi_lb_value + (x(jx) - 0.5d0 * dxx(jx))
CALL vanGenuchten_model_kr(psi_lb_value, theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), psi_s(jx, jy, jz),&
                              kr_lb)

jy = 1
jz = 1
! lower boundary
qx(0, jy, jz) = -xi*K_faces_x(0, jy, jz)/(dxx(1)*0.5d0)*((kr(1, jy, jz))*MAX(head(1, jy, jz) - head_lb, 0.0d0) + kr_lb*MIN(head(1, jy, jz) - head_lb, 0.0))

! internal faces
DO jx = 1, nx - 1
  qx(jx, jy, jz) = -xi*K_faces_x(jx, jy, jz)/(x(jx+1) - x(jx))*(kr(jx+1, jy, jz)*MAX(head(jx+1, jy, jz) - head(jx, jy, jz), 0.0) &
                      + kr(jx, jy, jz)*MIN(head(jx+1, jy, jz) - head(jx, jy, jz), 0.0))
END DO

! upper boundary
qx(nx, jy, jz) = qx_ub_value

END SUBROUTINE flux_Richards