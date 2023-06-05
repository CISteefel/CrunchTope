SUBROUTINE Jacobian_Richards_steady(nx, ny, nz, psi_lb_value, J)
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
REAL(DP), INTENT(OUT), DIMENSION(nx, nx) :: J

INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

! variables not declared in CrunchTope
REAL(DP), INTENT(IN)                                                   :: psi_lb_value
REAL(DP)                                                   :: head_lb
REAL(DP)                                                   :: kr_lb

!**************************************************
! physical parameters for Richards solver added by Toshiyuki Bandai, 2023, May
REAL(DP), PARAMETER :: mu = 1.0016d0 * 1.0E-3 * 86400.0d0 * 365.0d0 ! dynamics viscosity of water [Pa day] at 20 degC
REAL(DP), PARAMETER :: rho = 0.99823d0 * 1.0E3 ! density of water [kg m-3] at 20 degC
REAL(DP), PARAMETER :: g = 9.80665d0 * (86400.0d0 * 365.0d0) ** 2 ! gravitational acceleration [m day-2]
REAL(DP), PARAMETER :: xi = rho*g/mu ! constant used to solve the Richards equation
! End of edits by Toshiyuki Bandai, 2023, May
!**************************************************

jy = 1
jz = 1
DO jx = 1, nx
  CALL vanGenuchten_model(psi(jx, jy, jz), theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), &
                          theta(jx, jy, jz), kr(jx, jy, jz), dtheta(jx, jy, jz), dkr(jx, jy, jz))
  head(jx, jy, jz) = psi(jx, jy, jz) + x(jx)
END DO

jx = 1
jy = 1
jz = 1
head_lb = psi_lb_value + (x(jx) - 0.5d0 * dxx(jx))
CALL vanGenuchten_model_kr(psi_lb_value, theta_r(jx, jy, jz), theta_s(jx, jy, jz), VG_alpha(jx, jy, jz), VG_n(jx, jy, jz), &
                              kr_lb)

! evaluate Jacobian matrix
J = 0.0 ! initialize Jacobian matrix

! lower boundary
J(1, 1) = xi*K_faces_x(0, jy, jz)/(dxx(1)*0.5d0)*MERGE(dkr(1, jy, jz)*(head(1, jy, jz) - head_lb) + kr(1, jy, jz), kr_lb, head(1, jy, jz) - head_lb >= 0)  &
        - xi*K_faces_x(1, jy, jz)/(x(2) - x(1))*MERGE(-kr(2, jy, jz), dkr(1, jy, jz)*(head(2, jy, jz) - head(1, jy, jz)) - kr(1, jy, jz), head(2, jy, jz) - head(1, jy, jz) >= 0)
! Jacobian computation
J(1, 2) = - xi * K_faces_x(1, jy, jz) / (x(2) - x(1)) * MERGE(dkr(2, jy, jz)*(head(2, jy, jz) - head(1, jy, jz)) + kr(2, jy, jz), kr(1, jy, jz), head(2, jy, jz) - head(1, jy, jz) >= 0)


! interior nodes
inner: DO jx = 2, nx - 1
  J(jx, jx - 1) = xi * K_faces_x(jx-1, jy, jz) / (x(jx) - x(jx - 1))* &
              MERGE(-kr(jx, jy, jz), dkr(jx-1, jy, jz)*(head(jx, jy, jz) - head(jx-1, jy, jz)) - kr(jx-1, jy, jz), head(jx, jy, jz) - head(jx-1, jy, jz) >= 0)
  J(jx, jx) = xi * K_faces_x(jx-1, jy, jz) / (x(jx) - x(jx - 1)) &
          *MERGE(dkr(jx, jy, jz)*(head(jx, jy, jz) - head(jx-1, jy, jz)) + kr(jx, jy, jz), kr(jx-1, jy, jz), head(jx, jy, jz) - head(jx-1, jy, jz) >= 0) &    
          - xi * K_faces_x(jx, jy, jz) / (x(jx + 1) - x(jx)) &
          * MERGE(-kr(jx+1, jy, jz), dkr(jx, jy, jz)*(head(jx+1, jy, jz) - head(jx, jy, jz)) - kr(jx, jy, jz), head(jx+1, jy, jz) - head(jx, jy, jz) >= 0)
  J(jx, jx + 1) = - xi * K_faces_x(jx, jy, jz) / (x(jx + 1) - x(jx)) &
              * MERGE(dkr(jx+1, jy, jz)*(head(jx+1, jy, jz) - head(jx, jy, jz)) + kr(jx+1, jy, jz), kr(jx, jy, jz), head(jx+1, jy, jz) - head(jx, jy, jz) >= 0)
END DO inner

! upper boundary
J(nx, nx - 1) = xi * K_faces_x(nx-1, jy, jz) / (x(nx) - x(nx-1)) * &
                MERGE(-kr(nx, jy, jz), dkr(nx-1, jy, jz)*(head(nx, jy, jz) - head(nx-1, jy, jz)) - kr(nx-1, jy, jz), head(nx, jy, jz) - head(nx-1, jy, jz) >= 0)
J(nx, nx) = xi * K_faces_x(nx-1, jy, jz) / (x(nx) - x(nx-1)) &
            *MERGE(dkr(nx, jy, jz)*(head(nx, jy, jz) - head(nx-1, jy, jz)) + kr(nx, jy, jz), kr(nx-1, jy, jz), head(nx, jy, jz) - head(nx-1, jy, jz) >= 0)

  
END SUBROUTINE Jacobian_Richards_steady