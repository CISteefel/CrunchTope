SUBROUTINE solve_Richards_steady(nx, ny, nz)
! This subroutine solves the steady-state Richards equation using Newton's method with Armijo line search.
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

REAL(DP), DIMENSION(0:nx+1)                                    :: F_residual ! residual
REAL(DP), DIMENSION(0:nx+1, 0:nx+1)                                :: J ! Jacobian matrix
REAL(DP), DIMENSION(0:nx+1)                                    :: dpsi_Newton ! Newton step

! parameters for the linear solver
INTEGER(I4B)                                               :: info, lda, ldb, nrhs
INTEGER(I4B), DIMENSION(0:nx+1)                                :: ipiv

! parameters for Newtons' method for forward problem
REAL(DP), PARAMETER                                        :: tau_a = 1.0d-8
!REAL(DP), PARAMETER                                        :: tau_r = 1.0d-7
INTEGER(I4B), PARAMETER                                    :: maxitr = 1000

! variables for line search
REAL(DP)                                                   :: error_old, tol, alpha_line, lam, error_new
INTEGER(I4B)                                               :: no_backtrack, descent
INTEGER(I4B)                                               :: iteration
INTEGER(I4B)                                               :: total_line

! variables for checking water mass balance
REAL(DP)                                                   :: water_mass
REAL(DP)                                                   :: water_mass_error

! initialize parameters for linear solver
nrhs = 1
lda = nx + 2
ldb = nx + 2

! initialize parameters for Newton's method
iteration = 0
total_line = 0

! update tolerance
error_old = 1.0d20 ! at least one Newton iteration will be conducted
!tol = tau_r * error_old + tau_a
tol = tau_a ! this tolerance should be used for problems that need very high accuracy

! begin Newton's method
newton_loop: DO
  IF (error_old < tol) EXIT
  ! Evaluate the flux and the residual
  CALL flux_Richards_steady(nx, ny, nz)
  CALL residual_Richards_steady(nx, ny, nz, F_residual)
  ! Evaluate the Jacobian matrix
  CALL Jacobian_Richards_steady(nx, ny, nz, J)
  
  dpsi_Newton = -F_residual
  ! Solve the linear system
  CALL dgesv(nx+2, nrhs, j, lda, ipiv, dpsi_Newton, ldb, info)
  
  ! Armijo line search
  alpha_line = 1.0d-4
  lam = 1.0d0
  descent = 0
  no_backtrack = 0
  psi_prev = psi
  
  ! line search
  line: DO
    IF (descent /= 0 .OR. no_backtrack > 30) EXIT line
    ! update water potential
    DO jx = 0, nx+1
      psi(jx, ny, nz) = psi(jx, ny, nz) + lam * dpsi_Newton(jx)
    END DO
    ! evaluate the flux and the residual
    CALL flux_Richards_steady(nx, ny, nz)
    CALL residual_Richards_steady(nx, ny, nz, F_residual)
  
    ! update tolerance
    error_new = MAXVAL(ABS(F_residual))
    ! Check if Armijo conditions are satisfied
    Armijo: IF (error_new < error_old - error_old*alpha_line*lam) THEN
      error_old = error_new
      descent = 1
    ELSE Armijo
      lam = 0.5d0 * lam
      no_backtrack = no_backtrack + 1
      psi = psi_prev
          
    END IF Armijo
        
  END DO line
  
  IF (no_backtrack > 30) THEN
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

END SUBROUTINE solve_Richards_steady