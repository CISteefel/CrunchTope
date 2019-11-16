!!! *** Copyright Notice ***
!!! “CrunchFlow”, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory 
!!! (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.
!!! 
!!! If you have questions about your rights to use or distribute this software, please contact 
!!! Berkeley Lab's Innovation & Partnerships Office at  IPO@lbl.gov.
!!! 
!!! NOTICE.  This Software was developed under funding from the U.S. Department of Energy and the U.S. Government 
!!! consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting 
!!! on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, 
!!! prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
!!!
!!! *** License Agreement ***
!!! “CrunchFlow”, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory)
!!! subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved."
!!! 
!!! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!!! 
!!! (1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!!!
!!! (2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
!!! in the documentation and/or other materials provided with the distribution.
!!!
!!! (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy nor the names of 
!!! its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
!!!
!!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
!!! BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
!!! SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
!!! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
!!! OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
!!! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
!!! THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!!!
!!! You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or 
!!! performance of the source code ("Enhancements") to anyone; however, if you choose to make your
!!! Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without 
!!! imposing a separate written license agreement for such 
!!! Enhancements, then you hereby grant the following license: a  non-exclusive, royalty-free perpetual license to install, use, 
!!! modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or 
!!! derivative works thereof, in binary and source code form.

!!!      ****************************************

SUBROUTINE jacexchange_init(ncomp,nexchange,nexch_sec,neqn,nco)
USE crunchtype
USE params
USE concentration
USE solver

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                           :: neqn
INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nexchange
INTEGER(I4B), INTENT(IN)                           :: nexch_sec
INTEGER(I4B), INTENT(IN)                           :: nco

!  Internal variables and arrays

REAL(DP)                                           :: sum
REAL(DP)                                           :: sum1
REAL(DP)                                           :: sum2

INTEGER(I4B)                                       :: i
INTEGER(I4B)                                       :: i2
INTEGER(I4B)                                       :: nex
INTEGER(I4B)                                       :: ix
INTEGER(I4B)                                       :: ixcheck

fexch = 0.0

IF (iexc == 1 .OR. iexc == 3) THEN
  
  DO i = 1,ncomp
    IF (itype(i,nco) == 1) THEN
      DO i2 = 1,ncomp+nexchange
        sum = 0.0
        IF (equilibrate(i,nco)) THEN
          DO nex = 1,nexch_sec
            sum = sum + muexc(nex,i)*muexc(nex,i2)* spextmp10(nexchange+nex)
          END DO
        END IF
        fexch(i,i2) = sum
      END DO
    ELSE
      DO i2 = 1,ncomp+nexchange
        fexch(i,i2) = 0.0
      END DO
    END IF
  END DO
  
!  Part of "fexch" referring to total exchange concentrations (not aqueous components)
  
  DO ix = 1,nexchange
    DO i2 = 1,ncomp+nexchange
      sum2 = 0.0
      DO nex = 1,nexch_sec
        ixcheck = ixlink(nex)
        IF (ixcheck == ix) THEN
          sum2 = sum2 + muexc(nex,i2)*aexch(nex)
        END IF
        sum2 = sum2 + muexc(nex,i2)*aexch(nex)
      END DO
      fexch(ix+ncomp,i2) = sum2
    END DO
  END DO
  
ELSE    !                  Vanselow convention
  
  DO ix = 1,nexchange
    DO i2 = 1,ncomp+nexchange
      sum1 = 0.0
      sum2 = 0.0
      DO nex = 1,nexch_sec
        ixcheck = ixlink(nex)
        sum1 = sum1 + muexc(nex,ix+ncomp)*muexc(nex,i2)* aexch(nex)
        IF (ixcheck == ix) THEN
          sum2 = sum2 + muexc(nex,i2)*aexch(nex)
        END IF
      END DO
      fweight(ix+ncomp,i2) = sum1
      fexch(ix+ncomp,i2) = sum2
    END DO
  END DO
  
  
  DO i = 1,ncomp
    IF (itype(i,nco) == 1) THEN
      DO i2 = 1,ncomp+nexchange
        sum1 = 0.0
        sum2 = 0.0
        IF (equilibrate(i,nco)) THEN
          DO nex = 1,nexch_sec
            ix = ixlink(nex)
            sum1 = sum1 + muexc(nex,i)*muexc(nex,i2)* aexch(nex)*tec(ix)
            sum2 = sum2 - totexch(ix,nco)*muexc(nex,i)*aexch(nex)/  &
                (wt_aexch(ix)*wt_aexch(ix))*fweight(ix+ncomp,i2)
          END DO
        END IF
        fexch(i,i2) = sum1 + sum2
      END DO
    ELSE
      DO i2 = 1,ncomp+nexchange
        fexch(i,i2) = 0.0
      END DO
    END IF
  END DO
  
END IF

RETURN
END SUBROUTINE jacexchange_init
!***********************************************************
