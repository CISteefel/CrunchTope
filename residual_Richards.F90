SUBROUTINE residual_Richards(nx, ny, nz, dtflow, F_residual)
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

REAL(DP), INTENT(IN) :: dtflow
REAL(DP), DIMENSION(nx), INTENT(OUT) :: F_residual

F_residual= 0.0

jy = 1
jz = 1
! lower boundary
F_residual(1) = (dxx(1)/dtflow)*(theta(1, jy, jz) - theta_prev(1, jy, jz)) - qx(0, jy, jz) + qx(1, jy, jz)

! internal cells
DO jx = 2, nx-1
  F_residual(jx) = (dxx(jx)/dtflow)*(theta(jx, jy, jz) - theta_prev(jx, jy, jz)) - qx(jx-1, jy, jz) + qx(jx, jy, jz)
END DO

! upper boundary
F_residual(nx) = (dxx(nx)/dtflow)*(theta(nx, jy, jz) - theta_prev(nx, jy, jz)) - qx(nx-1, jy, jz) + qx(nx, jy, jz)

! add source/sink terms
IF (transpitimeseries) THEN
  DO i = 1, transpicells
    F_residual(nx+1-i) = F_residual(nx+1-i) - transpirate
  END DO
END IF

END SUBROUTINE residual_Richards