SUBROUTINE residual_Richards(nx, ny, nz, dtflow, F_residual)
! This subroutine calculates the residual of the Richards equation
! F = d theta / dt + div q - S (dimensionless)
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
REAL(DP), DIMENSION(nx), INTENT(OUT)                       :: F_residual

REAL(DP)                                                   :: upper_water_balance ! water balance in the top cell to prevent the cell from drying out

F_residual= 0.0

jy = 1
jz = 1
! lower boundary
F_residual(1) = theta(1, jy, jz) - theta_prev(1, jy, jz) + (qx(1, jy, jz) - qx(0, jy, jz))*dtflow/dxx(1)

! internal cells
DO jx = 2, nx-1
  F_residual(jx) = theta(jx, jy, jz) - theta_prev(jx, jy, jz) + (qx(jx, jy, jz) - qx(jx-1, jy, jz))*dtflow/dxx(jx)
END DO

! upper boundary
SELECT CASE (upper_BC_type)
CASE ('constant_flux', 'variable_flux')
  ! check if the prescribed flux (evaporation) is smaller than the remaining water in the top cell
  ! note that the flux is positive for evaporation
  upper_water_balance = (theta_prev(nx, jy, jz) - theta_r(nx, jy, jz)) - qx(nx, jy, jz)*dtflow/dxx(nx)  ! the unit is dimensionless
  IF (upper_water_balance <= 0.001d0) THEN
    ! the top cell is extremely dry, update the flux into the top cell to prevent the cell from drying out
    qx(nx, jy, jz) = (theta_prev(nx, jy, jz) - theta_r(nx, jy, jz) - 0.001d0)*dxx(nx)/dtflow
  END IF
  F_residual(nx) = theta(nx, jy, jz) - theta_prev(nx, jy, jz) + (qx(nx, jy, jz) - qx(nx-1, jy, jz))*dtflow/dxx(nx)
    
CASE DEFAULT
  F_residual(nx) = theta(nx, jy, jz) - theta_prev(nx, jy, jz) + (qx(nx, jy, jz) - qx(nx-1, jy, jz))*dtflow/dxx(nx)
END SELECT

! add source/sink terms
!IF (transpitimeseries) THEN
!  DO i = 1, transpicells
!    F_residual(nx+1-i) = F_residual(nx+1-i) - transpirate
!  END DO
!END IF

END SUBROUTINE residual_Richards