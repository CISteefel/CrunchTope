SUBROUTINE solve_Richards(nx, ny, nz, dtflow)
! This subroutine solves the Richards equation using Newton's method with Armijo line search.
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

REAL(DP), INTENT(IN)                                       :: dtflow ! time step for flow

REAL(DP), DIMENSION(nx)                                    :: F_residual ! residual
REAL(DP), DIMENSION(nx, nx)                                :: J ! Jacobian matrix
REAL(DP), DIMENSION(nx)                                    :: dpsi_Newton ! Newton step

! parameters for the linear solver
INTEGER(I4B)                                               :: info, lda, ldb, nrhs
INTEGER(I4B), DIMENSION(nx)                                :: ipiv

! parameters for Newtons' method for forward problem
REAL(DP), PARAMETER                                        :: tau_a = 1.0d-8
REAL(DP), PARAMETER                                        :: tau_r = 1.0d-7
INTEGER(I4B), PARAMETER                                    :: maxitr = 1000

! variables for line search
REAL(DP)                                                   :: error_old, tol, alpha_line, lam, error_new
REAL(DP), PARAMETER                                        :: tau_line_saerch = 1.0d-3 ! below this tolerance, the line search stops
INTEGER(I4B)                                               :: no_backtrack, descent
INTEGER(I4B)                                               :: iteration
INTEGER(I4B)                                               :: total_line

! variables for checking water mass balance
REAL(DP)                                                   :: water_mass
REAL(DP)                                                   :: water_mass_error

! initialize parameters for linear solver
nrhs = 1
lda = nx
ldb = nx

! initialize parameters for Newton's method
iteration = 0
total_line = 0

! Evaluate the flux and the residual
CALL flux_Richards(nx, ny, nz)

CALL residual_Richards(nx, ny, nz, dtflow, F_residual)

! update tolerance
error_old = MAXVAL(ABS(F_residual))
!tol = tau_r * error_old + tau_a
tol = tau_a ! this tolerance should be used for problems that need very high accuracy

! begin Newton's method
newton_loop: DO
  IF (error_old < tol .AND. iteration > 2) EXIT
  !IF (error_old < tol) EXIT
  ! Evaluate the Jacobian matrix
  CALL Jacobian_Richards(nx, ny, nz, dtflow, J)
  
  dpsi_Newton = -F_residual
  ! Solve the linear system
  CALL dgesv(nx, nrhs, j, lda, ipiv, dpsi_Newton, ldb, info)
  
  ! Armijo line search
  alpha_line = 1.0d-4
  lam = 1.0d0
  descent = 0
  no_backtrack = 0
  psi_prev = psi
  
  ! line search
  line: DO
    IF (descent /= 0 .OR. no_backtrack > 100) EXIT line
    ! update water potential
    DO jx = 1, nx
      psi(jx, ny, nz) = psi(jx, ny, nz) + lam * dpsi_Newton(jx)
    END DO
    ! evaluate the flux and the residual
    CALL flux_Richards(nx, ny, nz)
    CALL residual_Richards(nx, ny, nz, dtflow, F_residual)
  
    ! update tolerance
    error_new = MAXVAL(ABS(F_residual))
    ! Check if Armijo conditions are satisfied
    !Armijo: IF (error_new < error_old - error_old*alpha_line*lam) THEN
    ! line search stops when the error is below the tolerance (tau_line_saerch), this is because line search is not necessary when the error is small (and in fact line search may fail when the residual is very small)
    Armijo: IF (error_new < error_old - error_old*alpha_line*lam .OR. error_new < tau_line_saerch) THEN
      error_old = error_new
      descent = 1
    ELSE Armijo
      lam = 0.5d0 * lam
      no_backtrack = no_backtrack + 1
      psi = psi_prev
          
    END IF Armijo
        
  END DO line
  
  IF (no_backtrack > 100) THEN
    WRITE(*,*) ' Line search failed in the Richards solver. '
    READ(*,*)
    STOP
  END IF
  
  IF (iteration > maxitr) THEN
    WRITE(*,*) ' The Newton method failed to converge in the Richards solver. '
    READ(*,*)
    STOP
  END IF
  
  total_line = total_line + no_backtrack
  iteration = iteration + 1
  
  IF (Richards_print) THEN
  WRITE(*,120) iteration, tol, error_old
  120 FORMAT(1X, 'At the', I3, ' th Newton iteration, the tolerance is  ', ES14.4, ' , and the error is ', ES20.8)
  END IF
END DO newton_loop

IF (Richards_print) THEN
  WRITE(*,100) iteration, total_line
  100 FORMAT(1X, 'The Newton method needed ', I5, ' iterations with ', I3, ' line searches in the Richards solver. ')
END IF

!***********************************************************************************************************************************************
! evaluate water mass balance (this is not compatible with evaporaiton and transpiration for the "enviornmental_forcing" boundary condition.)
IF (Richards_print) THEN
  jy = 1
  jz = 1
  water_mass = 0.0d0
  DO jx = 1, nx
    water_mass = water_mass + dxx(jx)*(theta(jx, jy, jz) - theta_prev(jx, jy, jz)) ! how much water content is increased in this time step
  END DO
  
  SELECT CASE (upper_BC_type)
  CASE ('environmental_forcing')
    water_mass_error = water_mass - dtflow*(qx(0, jy, jz) - (infiltration_rate + evaporate))
    DO i = 1, transpicells
      water_mass_error = water_mass_error + dtflow*transpirate_cell(nx - i + 1)
    END DO
    water_mass_error = 100.0d0*water_mass_error/water_mass ! in percent
  CASE DEFAULT
    water_mass_error = 100.0d0*(water_mass - dtflow*(qx(0, jy, jz) - qx(nx, jy, jz)))/water_mass ! in percent
  END SELECT
  
  WRITE(*,110) water_mass, water_mass_error
110 FORMAT(1X, 'The water mass increase is ', ES14.4, ' m, and the water mass balance error is ', ES14.4, '%.')
    
! check water mass balance in each cell
!  DO jx = 1, nx - transpicells
!    water_mass = 0.0d0
!    water_mass = water_mass + dxx(jx)*(theta(jx, jy, jz) - theta_prev(jx, jy, jz)) ! how much water content is increased in this time step
!    water_mass_error = 100.0d0*(water_mass - dtflow*(qx(jx-1, jy, jz) - qx(jx, jy, jz)))/water_mass ! in percent
!    WRITE(*,130) water_mass, water_mass_error, jx
!130 FORMAT(1X, 'The water mass increase is ', ES14.4, ' m, and the water mass balance error is ', ES14.4, '% at the cell ', I4)
!  END DO
!  
  !READ(*,*)
END IF
!***********************************************************************************************************************************************

END SUBROUTINE solve_Richards