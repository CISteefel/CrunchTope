SUBROUTINE mineralfind(iunit5)
USE crunchtype
USE params

IMPLICIT NONE
 
INTEGER(I4B), INTENT(IN)                            :: iunit5
CHARACTER (LEN=mls)                                 :: dummy1

100 READ(iunit5,*,END=300) dummy1

IF (dummy1 == 'End of gases') THEN
  RETURN
END IF

GO TO 100

300 RETURN
END SUBROUTINE mineralfind
