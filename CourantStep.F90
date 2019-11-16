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

SUBROUTINE CourantStep(nx,ny,nz,dtmaxcour)

USE crunchtype
USE runtime
USE medium
USE transport
USE flow
USE modflowModule

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                 :: nx
INTEGER(I4B), INTENT(IN)                                 :: ny
INTEGER(I4B), INTENT(IN)                                 :: nz

REAL(DP), INTENT(OUT)                                    :: dtmaxcour

!  Internal variables and arrays

REAL(DP)                                                 :: dtmin,dt,cnv

INTEGER(I4B)                                             :: jx,jy,jz
INTEGER(I4B)                                             :: i

cnv = ModFlowCnv  !  Converts from MODFLOW time units to years (the CRUNCH time unit)

!   calculations for maximum time step based on courant condition

dtmin=1.e30
DO jz=1,nz
  DO jy=1,ny
    DO jx=1,nx
!fp! if_onproc( {#expr# dzz(jx,jy,jz) #});
      IF(qx(jx,jy,jz) == 0.0)THEN
        dt= 1.e30
      ELSE
        dt= ABS(courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*por(jx,jy,jz)/qx(jx,jy,jz))/cnv
!fp! end_onproc();
      END IF
      dtmin=MIN(dt,dtmin)
      IF(qy(jx,jy,jz) == 0.0)THEN
        dt=1.e30
      ELSE
        dt=ABS(courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*por(jx,jy,jz)/qy(jx,jy,jz))/cnv
      END IF
      dtmin=MIN(dt,dtmin)
    END DO
  END DO
END DO
WRITE(*,*)'qxy:',dtmin

!   wells

DO i = 1, nwells
  IF(q(i) == 0.0)THEN
    dt=1.e30
  ELSE
!fp! if_onproc( {#expr# por(jxWellLoc(i),jyWellLoc(i),jzWellLoc(i)) #});
    jx = jxWellLoc(i)
    jy = jyWellLoc(i)
    jz = jzWellLoc(i)
    dt=ABS(courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*  &
        por(jx,jy,jz)/q(i))/cnv
!fp! end_onproc();
  END IF
  dtmin=MIN(dt,dtmin)
END DO
WRITE(*,*)'well:',dtmin

!   constant head cells
DO i = 1, ncnh
  IF(qcnh(i) == 0.0)THEN
    dt=1.e30
  ELSE
!fp! if_onproc( {#expr# por(jxHeadLoc(i),jyHeadLoc(i),jzHeadLoc(i)) #});
    jx = jxHeadLoc(i)
    jy = jyHeadLoc(i)
    jz = jzHeadLoc(i)
    dt=ABS(courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*  &
        por(jx,jy,jz)/qcnh(i))/cnv
!fp! end_onproc();
  END IF
  dtmin=MIN(dt,dtmin)
END DO
WRITE(*,*)'cnh:',dtmin

!   recharge-evapotranspiration
DO jy=1,ny
  DO jx=1,nx
    IF(qrecharge(jx,jy)+qevt(jx,jy) == 0.0)THEN
      dt=1.e30
    ELSE
!fp! if_onproc( {#expr# por(jx,jy,1) #});
      dt=ABS( courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,1)*por(jx,jy,1) /  &
          (qrecharge(jx,jy) + qevt(jx,jy)) )/cnv
!fp! end_onproc();
    END IF
    dtmin=MIN(dt,dtmin)
  END DO
END DO

!   river
DO i=1,nrivers
  IF(qriver(i) == 0.0)THEN
    dt=1.e30
  ELSE
!fp! if_onproc( {#expr# por(jx,jy,1) #});
    jx = jxRiverLoc(i)
    jy = jyRiverLoc(i)
    jz = jzRiverLoc(i)
    dt=ABS(courfactor*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*  &
        por(jx,jy,jz)/qriver(i))/cnv
!fp! end_onproc();
  END IF
  dtmin=MIN(dt,dtmin)
END DO
WRITE(*,*)'river:',dtmin

dtmaxcour = dtmin

RETURN
END SUBROUTINE CourantStep
