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

interface
    subroutine residual_Richards(nx, ny, nz, dtflow, F_residual)
        USE crunchtype
        INTEGER(I4B), INTENT(IN)                                   :: nx
        INTEGER(I4B), INTENT(IN)                                   :: ny
        INTEGER(I4B), INTENT(IN)                                   :: nz
        REAL(DP), INTENT(IN)                                       :: dtflow
        REAL(DP), ALLOCATABLE, INTENT(INOUT)                                         :: F_residual(:) ! residual
    end subroutine residual_Richards
end interface

interface
    subroutine Jacobian_Richards(nx, ny, nz, dtflow, J)
        USE crunchtype
        INTEGER(I4B), INTENT(IN)                                   :: nx
        INTEGER(I4B), INTENT(IN)                                   :: ny
        INTEGER(I4B), INTENT(IN)                                   :: nz
        REAL(DP), INTENT(IN)                                       :: dtflow
        REAL(DP), ALLOCATABLE, INTENT(INOUT)                                         :: J(:, :) ! residual
    end subroutine Jacobian_Richards
end interface

INTEGER(I4B), INTENT(IN)                                   :: nx
INTEGER(I4B), INTENT(IN)                                   :: ny
INTEGER(I4B), INTENT(IN)                                   :: nz

INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz
INTEGER(I4B)                                               :: n_inner_cells
INTEGER(I4B)                                               :: n_total_cells
INTEGER(I4B)                                               :: n_bfaces
INTEGER(I4B)                                               :: n_xfaces
INTEGER(I4B)                                               :: n_yfaces

REAL(DP), INTENT(IN)                                       :: dtflow ! time step for flow

REAL(DP), ALLOCATABLE                                      :: F_residual(:) ! residual
REAL(DP), ALLOCATABLE                                      :: J(:, :) ! Jacobian matrix
REAL(DP), ALLOCATABLE                                      :: dpsi_Newton(:) ! Newton step

! parameters for the linear solver
INTEGER(I4B)                                               :: info, lda, ldb, nrhs
INTEGER(I4B), ALLOCATABLE                         :: ipiv(:)

! parameters for Newtons' method for forward problem
REAL(DP), PARAMETER                                        :: tau_a = 1.0d-8
!REAL(DP), PARAMETER                                        :: tau_r = 1.0d-7
INTEGER(I4B), PARAMETER                                    :: maxitr_Newton = 100
INTEGER(I4B), PARAMETER                                    :: maxitr_line_search = 30
! variables for line search
REAL(DP)                                                   :: error_old, tol, alpha_line, lam, error_new
INTEGER(I4B)                                               :: no_backtrack, descent
INTEGER(I4B)                                               :: iteration
INTEGER(I4B)                                               :: total_line

!*********************************************************************
! allocate arrays
IF (nx > 1 .AND. ny == 1 .AND. nz == 1) THEN ! one-dimensional problem
  IF (ALLOCATED(F_residual)) THEN
    DEALLOCATE(F_residual)
    ALLOCATE(F_residual(0:nx+1))
  ELSE
    ALLOCATE(F_residual(0:nx+1))
  END IF

  IF (ALLOCATED(J)) THEN
    DEALLOCATE(J)
    ALLOCATE(J(0:nx+1, 0:nx+1))
  ELSE
    ALLOCATE(J(0:nx+1, 0:nx+1))
  END IF

  IF (ALLOCATED(dpsi_Newton)) THEN
    DEALLOCATE(dpsi_Newton)
    ALLOCATE(dpsi_Newton(0:nx+1))
  ELSE
    ALLOCATE(dpsi_Newton(0:nx+1))
  END IF
  
  
  ! initialize parameters for linear solver
  nrhs = 1
  lda = nx + 2
  ldb = nx + 2


!*********************************************************************
ELSE IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
  IF (domain_shape_flow == 'regular') THEN
    n_inner_cells = nx * ny
    n_xfaces = nx*(ny+1) ! number of faces
    n_yfaces =  (nx+1)*ny
    n_bfaces = 2*nx + 2*ny
    n_total_cells = n_inner_cells + n_bfaces
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Currently, only regular spatial domain is supported.'
    WRITE(*,*)
    READ(*,*)
    STOP
        
  END IF
  
  
  
  IF (ALLOCATED(F_residual)) THEN
    DEALLOCATE(F_residual)
    ALLOCATE(F_residual(n_total_cells))
  ELSE
    ALLOCATE(F_residual(n_total_cells))
  END IF

  IF (ALLOCATED(J)) THEN
    DEALLOCATE(J)
    ALLOCATE(J(n_total_cells, n_total_cells))
  ELSE
    ALLOCATE(J(n_total_cells, n_total_cells))
  END IF

  IF (ALLOCATED(dpsi_Newton)) THEN
    DEALLOCATE(dpsi_Newton)
    ALLOCATE(dpsi_Newton(n_total_cells))
  ELSE
    ALLOCATE(dpsi_Newton(n_total_cells))
  END IF
  
  IF (ALLOCATED(ipiv)) THEN
    DEALLOCATE(ipiv)
    ALLOCATE(ipiv(n_total_cells))
  ELSE
    ALLOCATE(ipiv(n_total_cells))
  END IF
  
  nrhs = 1
  lda = n_total_cells
  ldb = n_total_cells


ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Currently, three-dimensional Richards solver is supported.'
  WRITE(*,*)
  READ(*,*)
  STOP
  
END IF

! initialize parameters for Newton's method
iteration = 0
total_line = 0

! update tolerance
error_old = 1.0d20 ! at least one Newton iteration will be conducted
!tol = tau_r * error_old + tau_a
tol = tau_a

! begin Newton's method
newton_loop: DO
  IF (error_old < tol) EXIT
  ! Evaluate the flux and the residual
  CALL flux_Richards(nx, ny, nz)
  
  IF (Richards_steady) THEN
    CONTINUE
    !CALL residual_Richards_steady(nx, ny, nz, dtflow, F_residual)
    ! Evaluate the Jacobian matrix
    !CALL Jacobian_Richards_steady(nx, ny, nz, dtflow, J)
  ELSE
    CALL residual_Richards(nx, ny, nz, dtflow, F_residual)
    ! Evaluate the Jacobian matrix
    CALL Jacobian_Richards(nx, ny, nz, dtflow, J)
  END IF
  
  dpsi_Newton = -F_residual
  ! Solve the linear system
  CALL dgesv(n_total_cells, nrhs, J, lda, ipiv, dpsi_Newton, ldb, info)
  
  ! Armijo line search
  alpha_line = 1.0d-4
  lam = 1.0d0
  descent = 0
  no_backtrack = 0
  psi_prev = psi
  
  ! line search
  line: DO
    IF (descent /= 0 .OR. no_backtrack > maxitr_line_search) EXIT line
    ! update water potential
    !IF (nx > 1 .AND. ny == 1 .AND. nz == 1) THEN ! one-dimensional problem
    !  
    !  DO jx = 0, nx+1
    !    IF (ABS(lam * dpsi_Newton(jx)) > dpsi_max) THEN
    !      psi(jx, ny, nz) = psi(jx, ny, nz) + SIGN(dpsi_max, dpsi_Newton(jx))
    !    ELSE
    !      psi(jx, ny, nz) = psi(jx, ny, nz) + lam * dpsi_Newton(jx)
    !    END IF
    !  
    !  END DO
    IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
      DO i = 1, n_total_cells
        jx = cell_to_coordinate(i, 1)
        jy = cell_to_coordinate(i, 2)
        IF (ABS(lam * dpsi_Newton(i)) > dpsi_max) THEN
          psi(jx, jy, nz) = psi(jx, jy, nz) + SIGN(dpsi_max, dpsi_Newton(i))
        ELSE
          psi(jx, jy, nz) = psi(jx, jy, nz) + lam * dpsi_Newton(i)
        END IF
      END DO
    END IF
    
    ! evaluate the flux and the residual
    CALL flux_Richards(nx, ny, nz)
    IF (Richards_steady) THEN
      CONTINUE
      !CALL residual_Richards_steady(nx, ny, nz, dtflow, F_residual)
    ELSE
      CALL residual_Richards(nx, ny, nz, dtflow, F_residual)
    END IF
  
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
  
  !IF (no_backtrack > maxitr_line_search) THEN
  !  WRITE(*,*) ' Line search failed in the Richards solver. '
  !  READ(*,*)
  !  STOP
  !END IF
  
  IF (iteration > maxitr_Newton) THEN
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

END SUBROUTINE solve_Richards