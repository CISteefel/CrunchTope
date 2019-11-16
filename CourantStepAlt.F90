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
    
SUBROUTINE CourantStepAlt(nx,ny,nz,dtmaxcour)
USE crunchtype
USE runtime
USE medium
USE transport
USE flow
USE temperature
USE modflowModule
USE ReadFlow

IMPLICIT NONE    

INTEGER(I4B), INTENT(IN)                                 :: nx
INTEGER(I4B), INTENT(IN)                                 :: ny
INTEGER(I4B), INTENT(IN)                                 :: nz

REAL(DP), INTENT(OUT)                                    :: dtmaxcour

!  Internal variables and arrays

INTEGER(I4B)                                             :: jx
INTEGER(I4B)                                             :: jy
INTEGER(I4B)                                             :: jz
INTEGER(I4B)                                             :: i

REAL(DP)                                                 :: dtemp
REAL(DP)                                                 :: porsatro
REAL(DP)                                                 :: dtmin

IF (xflow .OR. yflow .OR. zflow) THEN

  IF (ReadNuft) THEN        !  Mass flux (kg/m**2/yr) read from NUFT
    dtmin = 1000000.0

    DO jz = 0,nz
      DO jy = 0,ny
        DO jx = 0,nx
          IF (jy /= 0 .AND. jz /= 0) THEN
            IF (qx(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
              dtemp = courfactor*dxx(jx)*porsatro/ABS(qx(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
          IF (jx /= 0 .AND. jz /= 0) THEN
            IF (qy(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
              dtemp = courfactor*dyy(jy)*porsatro/ABS(qy(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
          IF (jx /= 0 .AND. jy /= 0) THEN
            IF (qz(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
              dtemp = courfactor*dzz(jx,jy,jz)*porsatro/ABS(qz(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
        END DO
      END DO
    END DO

!   Check the recharge for the Courant criteria

    DO jy=1,ny
      DO jx=1,nx
        porsatro = por(jx,jy,1)*satliq(jx,jy,1)*ro(jx,jy,1)
        IF(qrecharge(jx,jy) == 0.0)THEN
          dtemp=1.e30
        ELSE
!fp! if_onproc( {#expr# por(jx,jy,1) #});
          dtemp=ABS( 0.5*dxx(jx)*dyy(jy)*dzz(jx,jy,1)*porsatro/qrecharge(jx,jy)  )
!fp! end_onproc();
        END IF
        dtmin=MIN(dtemp,dtmin)
      END DO
    END DO

    DO i = 1, NumSourceTerms
      jx = jxNuftSource(i)
      jy = jyNuftSource(i)
      jz = jzNuftSource(i)
      porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz)
      IF(qg(jx,jy,jz) == 0.0)THEN
        dtemp=1.e30
      ELSE
!fp! if_onproc( {#expr# por(jx,jy,jz) #});
        dtemp=ABS(0.5*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)*porsatro/qg(jx,jy,jz))
!fp! end_onproc();
      END IF
      dtmin=MIN(dtemp,dtmin)
    END DO

    dtmaxcour = dtmin

!fp! reduce( {#ident# dtmaxcour #}, min_op);

  ELSE                   !  Volumetric flux used

    dtmin = 1000000.0
    DO jz = 0,nz
      DO jy = 0,ny
        DO jx = 0,nx
          IF (jy /= 0 .AND. jz /= 0) THEN
            IF (qx(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)
              dtemp = courfactor*dxx(jx)*porsatro/ABS(qx(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
          IF (jx /= 0 .AND. jz /= 0) THEN
            IF (qy(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)
              dtemp = courfactor*dyy(jy)*porsatro/ABS(qy(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
          IF (jx /= 0 .AND. jy /= 0) THEN
            IF (qz(jx,jy,jz) /= 0.0d0) THEN
              porsatro = por(jx,jy,jz)*satliq(jx,jy,jz)
              dtemp = courfactor*dzz(jx,jy,jz)*porsatro/ABS(qz(jx,jy,jz))
              dtmin = MIN(dtemp,dtmin)
            END IF
          END IF
        END DO
      END DO
    END DO
    dtmaxcour = dtmin

!fp! reduce( {#ident# dtmaxcour #}, min_op);

  END IF
END IF 

RETURN
END SUBROUTINE CourantStepAlt
