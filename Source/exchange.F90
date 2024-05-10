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

SUBROUTINE exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
USE crunchtype
USE params
USE concentration
USE transport
USE medium
USE temperature

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: ncomp
INTEGER(I4B), INTENT(IN)                                   :: nexchange
INTEGER(I4B), INTENT(IN)                                   :: nexch_sec
INTEGER(I4B), INTENT(IN)                                   :: jx
INTEGER(I4B), INTENT(IN)                                   :: jy
INTEGER(I4B), INTENT(IN)                                   :: jz

!  Internal variables

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: nex
INTEGER(I4B)                                               :: i


REAL(DP)                                                   :: activity
REAL(DP)                                                   :: sum
REAL(DP)                                                   :: check
REAL(DP)                                                   :: sumtemp
REAL(DP)                                                   :: exchangetemp

IF (iexc == 1) THEN        ! Gaines-Thomas convention
  sumactivity = 0.0
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    exchangetemp = exchangesites(ix,jx,jy,jz)



      sum = 0.0
      DO i = 1,ncomp
        sum = sum + muexc(nex,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
      END DO


    sum = sum + muexc(nex,ix+ncomp)*spex(ix,jx,jy,jz)
    
    activity = EXP(-keqexc(nex) + sum + bfit(nex)*sion(jx,jy,jz) )
    aexch(nex) = activity
    spex10(nex+nexchange,jx,jy,jz) = activity*exchangetemp/(muexc(nex,ix+ncomp))  
    sumactivity(ix) = sumactivity(ix) + aexch(nex)
    spex(nex+nexchange,jx,jy,jz) = LOG(spex10(nex+nexchange,jx,jy,jz))
  END DO
ELSE IF (iexc == 2) THEN   ! Vanselow convention
  DO nex = 1,nexch_sec
    ix = ixlink(nex)



      sum = 0.0
      DO i = 1,ncomp
        sum = sum + muexc(nex,i)*(sp(i,jx,jy,jz)+gam(i,jx,jy,jz))
      END DO


    sum = sum + muexc(nex,ix+ncomp)*spex(ix,jx,jy,jz)
    aexch(nex) = EXP(-keqexc(nex) + sum + bfit(nex)*sion(jx,jy,jz))
  END DO
  wt_aexch = 0.0
  sumactivity = 0.0
  DO nex = 1,nexch_sec
    ix = ixlink(nex)
    wt_aexch(ix) = wt_aexch(ix) + muexc(nex,ix+ncomp)*aexch(nex)
    sumactivity(ix) = sumactivity(ix) + aexch(nex)
  END DO
  DO ix = 1,nexchange
    tec(ix) = exchangesites(ix,jx,jy,jz)/wt_aexch(ix)
  END DO
  DO nex=1,nexch_sec
    ix = ixlink(nex)
    spex10(nex+nexchange,jx,jy,jz) = aexch(nex)*tec(ix)
    spex(nex+nexchange,jx,jy,jz) = LOG(spex10(nex+nexchange,jx,jy,jz))
  END DO
ELSE IF (iexc == 3) THEN   ! Gapon convention--same as Gaines-Thomas, but with different reaction stoichiometry
ELSE
  CONTINUE
END IF

RETURN
END SUBROUTINE exchange
!  **************************************************************
