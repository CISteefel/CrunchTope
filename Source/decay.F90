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
    
SUBROUTINE decay(ncomp,ikin,delt,jx,jy,jz)
USE crunchtype
USE params
USE mineral
USE concentration

IMPLICIT NONE
!fp! auto_par_loops = 0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                               :: ncomp
INTEGER(I4B), INTENT(IN)                               :: ikin
REAL(DP), INTENT(IN)                                   :: delt
INTEGER(I4B), INTENT(IN)                               :: jx
INTEGER(I4B), INTENT(IN)                               :: jy
INTEGER(I4B), INTENT(IN)                               :: jz

!  Internal variables and arrays

INTEGER(I4B)                                           :: ir
INTEGER(I4B)                                           :: kd
INTEGER(I4B)                                           :: k
INTEGER(I4B)                                           :: np
INTEGER(I4B)                                           :: i

REAL(DP)                                               :: rdecay
REAL(DP)                                               :: voltemporary
REAL(DP)                                               :: frx
REAL(DP)                                               :: denom

DO ir = 1,ikin

  IF (iaqtype(ir) == 4) then
    DO kd = 1,nrad_decay(ir)
      k = kradpoint(kd,ir)
      np = npradpoint(kd,ir)
      IF (iraddecay(ir) == 0) THEN
        WRITE(*,*) ' Error in decay reactions:  parent species in mineral not found'
        READ(*,*)
        STOP
      END IF
      IF (dppt(k,jx,jy,jz) > 0.0) THEN
        voltemporary = volmol(k)*dppt(k,jx,jy,jz)*delt
        denom = voltemporary + volfx(k,jx,jy,jz)
        IF (denom == 0.0) then      !  Only happens with zero molar volume
          frx = 1.0
        ELSE
          frx = voltemporary/denom
        END IF
        DO i = 1,ncomp
          mumin_decay(np,k,i,jx,jy,jz) = frx*mumin(np,k,i) + (1.0-frx)*mumin_decay(np,k,i,jx,jy,jz)
          rdecay = mukin(ir,i)*ratek(1,ir)*mumin_decay(np,k,iraddecay(ir),jx,jy,jz)*delt
          mumin_decay(np,k,i,jx,jy,jz) = mumin_decay(np,k,i,jx,jy,jz) +  rdecay
        END DO
      ELSE
        IF (volfx(k,jx,jy,jz) > 0.0) THEN
          DO i = 1,ncomp
            rdecay = mukin(ir,i)*ratek(1,ir)*mumin_decay(np,k,iraddecay(ir),jx,jy,jz)*delt
            mumin_decay(np,k,i,jx,jy,jz) = mumin_decay(np,k,i,jx,jy,jz) + rdecay
          END DO
        END IF
      END IF
    END DO
  END IF

END DO
    

END SUBROUTINE decay
