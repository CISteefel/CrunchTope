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

REAL(DP), DIMENSION(nx)                                    :: F_residual ! residual
REAL(DP), DIMENSION(nx, nx)                                :: J ! Jacobian matrix
REAL(DP), DIMENSION(nx)                                    :: dpsi_Newton ! Newton step

! parameters for the linear solver
INTEGER(I4B)                                               :: info, lda, ldb, nrhs
INTEGER(I4B), DIMENSION(nx)                                :: ipiv

! parameters for Newtons' method for forward problem
! Because the residual has a unit of m year^-1, the meaning of the tolerance for the residual is different from the time-dependent case (dimensionless residula is used).
REAL(DP), PARAMETER                                        :: tau_a = 1.0d-7 ! tolerance for the residual
REAL(DP), PARAMETER                                        :: tau_r = 1.0d-7
INTEGER(I4B), PARAMETER                                    :: maxitr = 1000

! variables for line search
REAL(DP)                                                   :: error_old, tol, alpha_line, lam, error_new
INTEGER(I4B)                                               :: no_backtrack, descent
INTEGER(I4B)                                               :: iteration
INTEGER(I4B)                                               :: total_line

! initialize parameters for linear solver
nrhs = 1
lda = nx
ldb = nx

! initialize parameters for Newton's method
iteration = 0
total_line = 0

! Evaluate the flux and the residual
CALL flux_Richards_steady(nx, ny, nz)

CALL residual_Richards_steady(nx, ny, nz, F_residual)

! update tolerance
error_old = MAXVAL(ABS(F_residual))
tol = tau_r * error_old + tau_a

! begin Newton's method
newton_loop: DO
  IF (error_old < tol) EXIT
  ! Evaluate the Jacobian matrix
  CALL Jacobian_Richards_steady(nx, ny, nz, J)
  
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
  
  IF (no_backtrack > 100) THEN
    WRITE(*,*) ' Line search failed in the steady-state Richards solver. '
    READ(*,*)
    STOP
  END IF
  
  total_line = total_line + no_backtrack
  iteration = iteration + 1
END DO newton_loop

IF (iteration > maxitr) THEN
  WRITE(*,*) ' The Newton method failed to converge in the steady-state Richards solver. '
  READ(*,*)
  STOP
END IF

IF (Richards_print) THEN
  WRITE(*,100) iteration, total_line
  100 FORMAT(1X, 'The Newton method needed ', I3, ' iterations with ', I3, ' line searches in the steady-state Richards solver. ')
END IF

END SUBROUTINE solve_Richards_steady