SUBROUTINE solve_Richards(nx, ny, nz, qx_lb_value, qx_ub_value, dtflow)
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

! variables not declared in CrunchTope
REAL(DP), INTENT(IN)                                                   :: qx_lb_value
REAL(DP), INTENT(IN)                                                   :: qx_ub_value

REAL(DP), DIMENSION(nx) :: F_residual ! residual
REAL(DP), DIMENSION(nx, nx) :: J ! Jacobian matrix
REAL(DP), DIMENSION(nx) :: dpsi_Newton ! Newton step

! linear solver
INTEGER(I4B) :: info, lda, ldb, nrhs
INTEGER(I4B), DIMENSION(nx) :: ipiv

! parameters for Newtons' method for forward problem
REAL(DP), PARAMETER :: tau_a = 1.0d-7
REAL(DP), PARAMETER :: tau_r = 1.0d-7
INTEGER(I4B), PARAMETER :: maxitr = 1000

! parameters for line search
REAL(DP) :: error_old, tol, alpha_line, lam, error_new
INTEGER(I4B) :: no_backtrack, descent
INTEGER(I4B)  :: iteration
INTEGER(I4B) :: total_line

! initialize parameters for linear solver
nrhs = 1
lda = nx
ldb = nx

! initialize parameters for Newton's method
iteration = 0
total_line = 0

! Evaluate the residual
CALL flux_Richards_noflow(nx, ny, nz, qx_lb_value, qx_ub_value)

CALL residual_Richards(nx, ny, nz, dtflow, F_residual)

! update tolerance
error_old = MAXVAL(ABS(F_residual))
tol = tau_r * error_old + tau_a

newton_loop: DO
  IF (error_old < tol) EXIT
  ! Evaluated the Jacobian matrix
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
    ! evaluate the residual
    CALL flux_Richards_noflow(nx, ny, nz, qx_lb_value, qx_ub_value)
    CALL residual_Richards(nx, ny, nz, dtflow, F_residual)
  
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
    WRITE(*, 100)
100 FORMAT(1X, 'Line search failed.')
    READ(*,*)
  END IF
  
  total_line = total_line + no_backtrack
  iteration = iteration + 1
END DO newton_loop

IF (iteration > maxitr) THEN
  WRITE(*, 110)
  110 FORMAT(1X, 'The Newton method did not converge.')
END IF


WRITE(*,120) iteration, total_line
120 FORMAT(1X, 'The Newton method needed ', I3, ' iterations with ', I3, ' line searches.')
END SUBROUTINE solve_Richards