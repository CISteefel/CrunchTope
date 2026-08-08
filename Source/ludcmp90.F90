SUBROUTINE ludcmp90(a,indx,d,n)
USE crunchtype
USE CrunchFunctions

IMPLICIT NONE

! ********** INTERFACE BLOCKS ******
interface
  SUBROUTINE swap(a,b)
  USE crunchtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(INOUT)              :: a,b
  END SUBROUTINE swap
END interface

interface
  SUBROUTINE nrerror(string)
  USE crunchtype
  IMPLICIT NONE
  CHARACTER (LEN=*), intent(in)                 :: string
  END SUBROUTINE nrerror
END interface
! **********************************

!  External arrays and variables

REAL(DP), DIMENSION(:,:), intent(in out)                   :: a
INTEGER(I4B), DIMENSION(:), intent(out)                    :: indx
INTEGER(I4B), INTENT(IN)                                   :: n
REAL(DP), intent(out)                                      :: d

!  Internal arrays and variables

REAL(DP), DIMENSION(size(a,1))                             :: vv
!REAL(DP), DIMENSION(n)                                     :: vv
REAL(DP), PARAMETER                                        :: tiny=1.0E-20

INTEGER(I4B)                                               :: j
INTEGER(I4B)                                               :: imax

d=1.0
vv=maxval(ABS(a),dim=2)
IF (any(vv == 0.0)) THEN
  WRITE(*,*) '  ludcmp90: zero row in matrix (singular) -- signalling decbt90'
  d = 0.0_DP  ! sentinel: tells decbt90 to propagate error and reduce timestep
  RETURN
END IF
vv=1.0_DP/vv
DO j=1,n
  imax=(j-1)+imaxloc(vv(j:n)*ABS(a(j:n,j)))
  IF (j /= imax) THEN
    CALL swap(a(imax,:),a(j,:))
    d=-d
    vv(imax)=vv(j)
  END IF
  indx(j)=imax
  IF (a(j,j) == 0.0) a(j,j)=tiny
  a(j+1:n,j)=a(j+1:n,j)/a(j,j)
  a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
END DO

END SUBROUTINE ludcmp90



