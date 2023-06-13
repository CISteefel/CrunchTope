
 
SUBROUTINE decbt90 (m,n,ier)
USE crunchtype
USE solver

IMPLICIT NONE

interface
SUBROUTINE ludcmp90(a,indx,d,n)
USE crunchtype
REAL(DP), DIMENSION(:,:), intent(in out)                   :: a
INTEGER(I4B), DIMENSION(:), intent(out)                    :: indx
INTEGER(I4B), INTENT(IN)                                   :: n
REAL(DP), intent(out)                                      :: d
END SUBROUTINE ludcmp90
END interface

interface
SUBROUTINE lubksb90(a,indx,b,n)
USE crunchtype
REAL(DP), DIMENSION(:,:), intent(in)                        :: a
INTEGER(I4B), DIMENSION(:), intent(in)                      :: indx
REAL(DP), DIMENSION(:), intent(in out)                      :: b
INTEGER(I4B), INTENT(IN)                                    :: n
END SUBROUTINE lubksb90
END interface

!  External variables

INTEGER(I4B), INTENT(IN)                                       :: m
INTEGER(I4B), INTENT(IN)                                       :: n
INTEGER(I4B), INTENT(OUT)                                      :: ier

!  Internal variables

INTEGER(I4B)                                                   :: nm1
INTEGER(I4B)                                                   :: nm2
INTEGER(I4B)                                                   :: k
INTEGER(I4B)                                                   :: j
INTEGER(I4B)                                                   :: i
INTEGER(I4B)                                                   :: l
INTEGER(I4B)                                                   :: km1

REAL(DP)                                                       :: det
REAL(DP)                                                       :: dpr


!-----------------------------------------------------------------------bt000040
! block-tridiagonal matrix decomposition routine.                       bt000050
! written by a. c. hindmarsh.                                           bt000060
! latest revision january 26, 1977  (ag)                                bt000070
! reference]  ucid-30150                                                bt000080
!             solution of block-tridiagonal systems of linear           bt000090
!             algebraic equations                                       bt000100
!             a.c. hindmarsh                                            bt000110
!             february 1977                                             bt000120
! the input matrix contains three blocks of elements in each block-row, bt000130
! including blocks in the (1,3) and (n,n-2) block positions.            bt000140
! decbt uses block gauss elimination and subroutines dec and sol        bt000150
! for solution of blocks.  partial pivoting is done within              bt000160
! block-rows only.                                                      bt000170
! input..                                                               bt000180
!     m = order of each block.                                          bt000190
!     n = number of blocks in each direction of the matrix.             bt000200
!         n must be 4 or more.  the complete matrix has order m*n.      bt000210
!     a = m by m by n array containing diagonal blocks.                 bt000220
!         a(i,j,k) contains the (i,j) element of the k-th block.        bt000230
!     b = m by m by n array containing the super-diagonal blocks        bt000240
!         (in b(*,*,k) for k = 1,...,n-1) and the block in the (n,n-2)  bt000250
!         block position (in b(*,*,n)).                                 bt000260
!     c = m by m by n array containing the subdiagonal blocks           bt000270
!         (in c(*,*,k) for k = 2,3,...,n) and the block in the          bt000280
!         (1,3) block position (in c(*,*,1)).                           bt000290
!    ip = integer array of length m*n for working storage.              bt000300
! output..                                                              bt000310
! a,b,c = m by m by n arrays containing the block lu decomposition      bt000320
!         of the input matrix.                                          bt000330
!    ip = m by n array of pivot information.  ip(*,k) contains          bt000340
!         information for the k-th digonal block.                       bt000350
!   ier = 0  if no trouble occurred, or                                 bt000360
!       = -1 if the input value of m or n was illegal, or               bt000370
!       = k  if a singular matrix was found in the k-th diagonal block. bt000380
! use solbt to solve the associated linear system.                      bt000390
! decbt calls subroutines  dec(m,m0,a,ip,ier)  and  sol(m,m0,a,y,ip)    bt000400
! for solution of m by m linear systems.                                bt000410
!-----------------------------------------------------------------------bt000420

!      integer nm1, nm2, km1,i,j,k,l                                     bt000430
!v    real dp                                                           bt000440

IF (m < 1 .OR. n < 4) GO TO 210
nm1 = n - 1
nm2 = n - 2
! process the first block-row. -----------------------------------------bt000480
CALL ludcmp90(aah(:,:,1),indexx(:,1),det,m)
ier = 0
k = 1
IF (ier /= 0) GO TO 200
DO j = 1,m
  CALL lubksb90(aah(:,:,1),indexx(:,1),bbh(:,j,1),m)
  CALL lubksb90(aah(:,:,1),indexx(:,1),cch(:,j,1),m)
END DO
! adjust b(*,*,2). -----------------------------------------------------bt000560
DO j = 1,m
  DO i = 1,m
    dpr = 0.
    DO l = 1,m
      dpr = dpr + cch(i,l,2)*cch(l,j,1)
    END DO
    bbh(i,j,2) = bbh(i,j,2) - dpr
  END DO
END DO
! main loop.  process block-rows 2 to n-1. -----------------------------bt000650
DO k = 2,nm1
  km1 = k - 1
  DO j = 1,m
    DO i = 1,m
      dpr = 0.
      DO l = 1,m
        dpr = dpr + cch(i,l,k)*bbh(l,j,km1)
      END DO
      aah(i,j,k) = aah(i,j,k) - dpr
    END DO
  END DO
  CALL ludcmp90(aah(:,:,k),indexx(:,k),det,m)
  IF (ier /= 0) GO TO 200
  DO j = 1,m
    CALL lubksb90(aah(:,:,k),indexx(:,k),bbh(:,j,k),m)
  END DO
END DO
! process last block-row and return. -----------------------------------bt000810
DO j = 1,m
  DO i = 1,m
    dpr = 0.
    DO l = 1,m
      dpr = dpr + bbh(i,l,n)*bbh(l,j,nm2)
    END DO
    cch(i,j,n) = cch(i,j,n) - dpr
  END DO
END DO
DO j = 1,m
  DO i = 1,m
    dpr = 0.
    DO l = 1,m
      dpr = dpr + cch(i,l,n)*bbh(l,j,nm1)
    END DO
    aah(i,j,n) = aah(i,j,n) - dpr
  END DO
END DO
CALL ludcmp90(aah(:,:,n),indexx(:,n),det,m)
k = n
IF (ier /= 0) GO TO 200
RETURN
! error returns. -------------------------------------------------------bt001020
200  ier = k
RETURN
210  ier = -1
RETURN
!-----------------------  end of subroutine decbt  ---------------------bt001070
END SUBROUTINE decbt90


