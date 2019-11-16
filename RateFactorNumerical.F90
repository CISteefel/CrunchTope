!******************        GIMRT98     ************************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:02:15
 
!************** (C) COPYRIGHT 1995,1998,1999 ******************
!*******************     C.I. Steefel      *******************
!                    All Rights Reserved

!  GIMRT98 IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR IMPLIED.
!  THE USER ASSUMES ALL RISKS OF USING GIMRT98. THERE IS NO CLAIM OF THE
!  MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO ANY OTHER
!  WORKSTATIONS
!**********************************************************************

SUBROUTINE RateFactorNumerical(ncomp,nrct,jx,jy,jz,np,k,sppTMP,RateFactor)
USE crunchtype
USE params
USE concentration, ONLY: gam,ulab
USE mineral
USE medium
USE temperature

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                        :: ncomp
INTEGER(I4B), INTENT(IN)                                        :: nrct
INTEGER(I4B), INTENT(IN)                                        :: jx
INTEGER(I4B), INTENT(IN)                                        :: jy
INTEGER(I4B), INTENT(IN)                                        :: jz
INTEGER(I4B), INTENT(IN)                                        :: np
INTEGER(I4B), INTENT(IN)                                        :: k

REAL(DP),DIMENSION(:),INTENT(IN)                                :: sppTMP

REAL(DP), INTENT(OUT)                                           :: RateFactor

!  Internal variables and arrays

REAL(DP)                                                        :: sumiap
REAL(DP)                                                        :: silogTMP
REAL(DP)                                                        :: CubicTerm
REAL(DP)                                                        :: CorrectedTerm

INTEGER(I4B)                                                    :: i

sumiap = 0.0D0
DO i = 1,ncomp
  IF (ulab(i) == 'H2O') THEN
    sumiap = sumiap + mumin(np,k,i)*(gam(i,jx,jy,jz))
  ELSE
    sumiap = sumiap + mumin(np,k,i)*(sppTMP(i)+gam(i,jx,jy,jz))
  END IF
END DO
silogTMP = (sumiap - keqmin(np,k,jx,jy,jz))/clg
  
IF (silogTMP >= 0.0d0) THEN
  RateFactor = 0.0d0
ELSE 
  CorrectedTerm = silogTMP + thresh(np,k)
  CubicTerm = CorrectedTerm*CorrectedTerm*CorrectedTerm/                 &
       ( thresh(np,k)*thresh(np,k)*thresh(np,k) )
  RateFactor = 0.5d0 - 0.75d0*CorrectedTerm/thresh(np,k) + 0.25d0*CubicTerm
END IF

   
RETURN
END SUBROUTINE RateFactorNumerical
!********************************************************
