SUBROUTINE residual_Richards_steady(nx, ny, nz, dtflow, F_residual)
! This subroutine calculates the residual of the steady-state Richards equation
! F = div q - S (dimensionless)
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

REAL(DP), INTENT(IN)                                       :: dtflow
REAL(DP), DIMENSION(0:nx + 1), INTENT(OUT)                 :: F_residual

REAL(DP)                                                   :: psi_b ! water potential at a boundary
REAL(DP)                                                   :: psi_grad_b ! water potential gradient at a boundary

!REAL(DP)                                                   :: water_balance ! water balance to prevent the cell from drying out
!REAL(DP)                                                   :: adjusted_extraction ! total water extraction (=evaporation + transpiration) adjusted to prevent the cell from drying out


F_residual= 0.0

jy = 1
jz = 1

! internal cells
DO jx = 1, nx
  F_residual(jx) = (qx(jx, jy, jz) - qx(jx-1, jy, jz))*dtflow/dxx_2(jx)
END DO

! boundary condition at the inlet (begin boundary condition)

SELECT CASE (x_begin_BC_type_steady)
CASE ('constant_dirichlet')
  psi_b = (psi(0, jy, jz)*dxx_2(1) + psi(1, jy, jz)*dxx_2(0))/(dxx_2(0) + dxx_2(1))
  F_residual(0) = psi_b - value_x_begin_BC_steady
  
  CONTINUE
  
CASE ('constant_neumann')
  psi_grad_b = (psi(1, jy, jz) - psi(0, jy, jz))/(x_2(1) - x_2(0))
  F_residual(0) = psi_grad_b - value_x_begin_BC_steady
  
CASE ('constant_flux')
  F_residual(0) = qx(0, jy, jz) - value_x_begin_BC_steady
  
CASE DEFAULT
  WRITE(*,*)
  WRITE(*,*) ' The boundary condition type ', x_begin_BC_type_steady, ' is not supported. '
  WRITE(*,*)
  READ(*,*)
  STOP
  
END SELECT


! boundary condition at the inlet (end boundary condition)
SELECT CASE (x_end_BC_type_steady)
CASE ('constant_dirichlet')
  psi_b = (psi(nx, jy, jz)*dxx_2(nx+1) + psi(nx+1, jy, jz)*dxx_2(nx))/(dxx_2(nx) + dxx_2(nx+1))
  F_residual(nx+1) = psi_b - value_x_end_BC_steady
  
CASE ('constant_neumann')
  psi_grad_b = (psi(nx+1, jy, jz) - psi(nx, jy, jz))/(x_2(nx+1) - x_2(nx))
  F_residual(nx+1) = psi_grad_b - value_x_end_BC_steady
  
CASE ('constant_flux')
  F_residual(nx+1) = qx(nx, jy, jz) - value_x_end_BC_steady
  
CASE DEFAULT
  WRITE(*,*)
  WRITE(*,*) ' The boundary condition type ', x_end_BC_type_steady, ' is not supported. '
  WRITE(*,*)
  READ(*,*)
  STOP
  
END SELECT

END SUBROUTINE residual_Richards_steady