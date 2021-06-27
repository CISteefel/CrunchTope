!!! *** Copyright Notice ***
!!! �CrunchFlow�, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory
!!! (subject to receipt of any required approvals from the U.S. Dept. of Energy).� All rights reserved.
!!!�
!!! If you have questions about your rights to use or distribute this software, please contact
!!! Berkeley Lab's Innovation & Partnerships Office at��IPO@lbl.gov.
!!!�
!!! NOTICE.� This Software was developed under funding from the U.S. Department of Energy and the U.S. Government
!!! consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting
!!! on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public,
!!! prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
!!!
!!! *** License Agreement ***
!!! �CrunchFlow�, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory)
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

SUBROUTINE watercontentRich(nx,ny,nz,dtyr)
USE crunchtype
USE params
USE medium
USE transport
USE temperature, ONLY: ro
USE flow
USE CrunchFunctions

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)                                       :: nx
INTEGER(I4B), INTENT(IN)                                       :: ny
INTEGER(I4B), INTENT(IN)                                       :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                          :: jx
INTEGER(I4B)                                                          :: jy
INTEGER(I4B)                                                          :: jz

!  ****** PARAMETERS  ****************************
REAL(DP)                                                              :: coef
REAL(DP)                                                              :: dt
REAL(DP), PARAMETER                                                   :: Ss=1.0d-6
REAL(DP), INTENT(IN)                                                  :: dtyr

CHARACTER (LEN=1)                                                     :: Coordinate

REAL(DP)                                                              :: pumpterm
INTEGER(I4B)                                                          :: npz

dt = dtyr * 365 * 86400
!   calculate darcy fluxes

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
        IF (activecellPressure(jx,jy,jz) == 1) THEN
            coef = 1.0d0 / (1.0d0 + Ss*(head(jx,jy,jz)-headOld(jx,jy,jz))/wcs(jx,jy,jz))
            wc(jx,jy,jz) = wcOld(jx,jy,jz)*coef + (dt*coef/dxx(jx))*(-qx(jx,jy,jz) + qx(jx-1,jy,jz))     &
                                        +(dt*coef/dyy(jy))*(-qy(jx,jy,jz) + qy(jx,jy-1,jz))  &
                                        + (dt*coef/dzz(jx,jy,jz))*(-qz(jx,jy,jz) + qz(jx,jy,jz-1))
            ! add source term if not along boundary
            IF (jx > 1 .AND. jx < nx) THEN
              
                IF (jy > 1 .AND. jy < ny .AND. activecellPressure(jx,jy-1,jz) == 1 .AND. qg(1,jx,jy,jz) /= 0.0) THEN
                  pumpterm = 0.0d0
                  DO npz = 1,npump(jx,jy,jz)
                    pumpterm = pumpterm + dt*qg(1,jx,jy,jz)/(secyr*dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
                  END DO
                  wc(jx,jy,jz) = wc(jx,jy,jz) + pumpterm
                END IF
                
            END IF

            ! IF (activecellPressure(jx,jy-1,jz) == 0) THEN
            !     WRITE(*,*) ' jx, wc, Ktop, qtop, Kbot, qbot = ',jx,wc(jx,jy,jz),Kfacy(jx,jy-1,jz),qy(jx,jy-1,jz),Kfacy(jx,jy,jz),qy(jx,jy,jz)
            ! END IF

            room(jx,jy,jz) = (wcs(jx,jy,jz) - wc(jx,jy,jz)) * dxx(jx) * dyy(jy) * dzz(jx,jy,jz)
            IF (room(jx,jy,jz) < 0.0) THEN
                room(jx,jy,jz) = 0.0d0
            ELSE IF (room(jx,jy,jz) > (wcs(jx,jy,jz)-wcr(jx,jy,jz)) * dxx(jx) * dyy(jy) * dzz(jx,jy,jz)) THEN
                room(jx,jy,jz) = (wcs(jx,jy,jz)-wcr(jx,jy,jz)) * dxx(jx) * dyy(jy) * dzz(jx,jy,jz)
            END IF
        END IF

    END DO
  END DO
END DO

DO jz = 1,nz
  DO jx = 1,nx
    IF (activecellPressure(jx,0,jz) == 0) THEN
        wc(jx,0,jz) = wch(jx,0,jz)
    ELSE
        wc(jx,0,jz) = wch(jx,0,jz)
    END IF
    IF (activecellPressure(jx,ny+1,jz) == 0) THEN
        wc(jx,ny+1,jz) = wch(jx,ny+1,jz)
    ELSE
        wc(jx,ny+1,jz) = wch(jx,ny+1,jz)
    END IF
  END DO
END DO

DO jz = 1,nz
  DO jy = 1,ny
    IF (activecellPressure(0,jy,jz) == 0) THEN
        wc(0,jy,jz) = wch(0,jy,jz)
    ELSE
        wc(0,jy,jz) = wch(0,jy,jz)
    END IF

    IF (activecellPressure(nx+1,jy,jz) == 0) THEN
        wc(nx+1,jy,jz) = wch(nx+1,jy,jz)
    ELSE
        wc(nx+1,jy,jz) = wch(nx+1,jy,jz)
    END IF
  END DO
END DO

DO jy = 1,ny
  DO jx = 1,nx
    IF (activecellPressure(jx,jy,0) == 0) THEN
        wc(jx,jy,0) = wch(jx,jy,0)
    ELSE
        wc(jx,jy,0) = wch(jx,jy,0)
    END IF
    IF (activecellPressure(jx,jy,nz+1) == 0) THEN
        wc(jx,jy,nz+1) = wch(jx,jy,nz+1)
    ELSE
        wc(jx,jy,nz+1) = wch(jx,jy,nz+1)
    END IF
  END DO
END DO


RETURN
END SUBROUTINE watercontentRich
