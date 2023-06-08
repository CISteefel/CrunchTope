SUBROUTINE residual_Richards_steady(nx, ny, nz, F_residual)
! This subroutine calculates the residual of the steady-state Richards equation
! F = d theta / dt + div q - S
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

REAL(DP), DIMENSION(nx), INTENT(OUT)                       :: F_residual

F_residual= 0.0

jy = 1
jz = 1
! lower boundary
F_residual(1) = qx(1, jy, jz) - qx(0, jy, jz)

! internal cells
DO jx = 2, nx-1
  F_residual(jx) = qx(jx, jy, jz) - qx(jx-1, jy, jz)
END DO

! upper boundary
F_residual(nx) = qx(nx, jy, jz) - qx(nx-1, jy, jz)

END SUBROUTINE residual_Richards_steady