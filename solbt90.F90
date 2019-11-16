

SUBROUTINE solbt90(m,n)
USE crunchtype
USE solver

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                       :: m
INTEGER(I4B), INTENT(IN)                                       :: n
!  *********************************************
interface
SUBROUTINE lubksb90(a,indx,b,n)
USE crunchtype
IMPLICIT NONE
REAL(DP),  DIMENSION(:,:), INTENT(IN)                          :: a
REAL(DP),  DIMENSION(:), INTENT(IN OUT)                        :: b
INTEGER(I4B),  DIMENSION(:),INTENT(IN)                         :: indx
INTEGER(I4B), INTENT(IN)                                       :: n
END SUBROUTINE lubksb90
END interface

!      dimension a(ndim,ndim,mx),b(ndim,ndim,mx),c(ndim,ndim,mx)
!      dimension y(ndim,mx)
!      dimension ip(ndim,mx)

!-----------------------------------------------------------------------bt001120
! solution of block-tridiagonal linear system.                          bt001130
! coefficient matrix must have been previously processed by decbt.      bt001140
! m, n, a, b, c, and ip  must not have been changed since call to decbt.bt001150
! written by a. c. hindmarsh.                                           bt001160
! input..                                                               bt001170
!     m = order of each block.                                          bt001180
!     n = number of blocks in each direction of matrix.                 bt001190
! a,b,c = m by m by n arrays containing block lu decomposition          bt001200
!         of coefficient matrix from decbt.                             bt001210
!    ip = m by n integer array of pivot information from decbt.         bt001220
!     y = array of length m*n containg the right-hand side vector       bt001230
!         (treated as an m by n array here).                            bt001240
! output..                                                              bt001250
!     y = solution vector, of length m*n.                               bt001260
! solbt makes calls to subroutine sol(m,m0,a,y,ip)                      bt001270
! for solution of m by m linear systems.                                bt001280
!-----------------------------------------------------------------------bt001290
!                                                                       bt001300

INTEGER(I4B)                                                   :: nm1
INTEGER(I4B)                                                   :: nm2
INTEGER(I4B)                                                   :: km1
INTEGER(I4B)                                                   :: i
INTEGER(I4B)                                                   :: j
INTEGER(I4B)                                                   :: k
INTEGER(I4B)                                                   :: kb
INTEGER(I4B)                                                   :: kp1

REAL(DP)                                                       :: dpr

nm1 = n - 1
nm2 = n - 2
! forward solution sweep. ----------------------------------------------bt001350
CALL lubksb90(aah(:,:,1),indexx(:,1),yh(:,1),m)
DO k = 2,nm1
  km1 = k - 1
  DO i = 1,m
    dpr = 0.
    DO j = 1,m
      dpr = dpr + cch(i,j,k)*yh(j,km1)
    END DO
    yh(i,k) = yh(i,k) - dpr
  END DO
  CALL lubksb90(aah(:,:,k),indexx(:,k),yh(:,k),m)
END DO
DO i = 1,m
  dpr = 0.
  DO j = 1,m
    dpr = dpr + cch(i,j,n)*yh(j,nm1) + bbh(i,j,n)*yh(j,nm2)
  END DO
  yh(i,n) = yh(i,n) - dpr
END DO
CALL lubksb90(aah(:,:,n),indexx(:,n),yh(:,n),m)
! backward solution sweep. ---------------------------------------------bt001540
DO kb = 1,nm1
  k = n - kb
  kp1 = k + 1
  DO i = 1,m
    dpr = 0.
    DO j = 1,m
      dpr = dpr + bbh(i,j,k)*yh(j,kp1)
    END DO
    yh(i,k) = yh(i,k) - dpr
  END DO
END DO
DO i = 1,m
  dpr = 0.
  DO j = 1,m
    dpr = dpr + cch(i,j,1)*yh(j,3)
  END DO
  yh(i,1) = yh(i,1) - dpr
END DO
RETURN
!-----------------------  end of subroutine solbt  ---------------------bt001720
END SUBROUTINE solbt90

