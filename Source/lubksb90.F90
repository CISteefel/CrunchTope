SUBROUTINE lubksb90(a,indx,b,n)
USE crunchtype

IMPLICIT NONE

!  External variables and arrays

REAL(DP),  DIMENSION(:,:), INTENT(IN)                          :: a
REAL(DP),  DIMENSION(:), INTENT(IN OUT)                        :: b

INTEGER(I4B),  DIMENSION(:),INTENT(IN)                         :: indx
INTEGER(I4B), INTENT(IN)                                       :: n

!  Internal variables and arrays

INTEGER(I4B)                                                   :: i
INTEGER(I4B)                                                   :: ii
INTEGER(I4B)                                                   :: ll

REAL(DP)                                                       :: summ

ii=0
DO i=1,n
  ll=indx(i)
  summ=b(ll)
  b(ll)=b(i)
  IF (ii /= 0) THEN
    summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
  ELSE IF (summ /= 0.0) THEN
    ii=i
  END IF
  b(i)=summ
END DO
DO i=n,1,-1
  b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
END DO
END SUBROUTINE lubksb90
