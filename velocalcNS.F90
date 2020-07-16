!!! *** Copyright Notice ***
!!! ìCrunchFlowî, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory
!!! (subject to receipt of any required approvals from the U.S. Dept. of Energy).† All rights reserved.
!!!†
!!! If you have questions about your rights to use or distribute this software, please contact
!!! Berkeley Lab's Innovation & Partnerships Office at††IPO@lbl.gov.
!!!†
!!! NOTICE.† This Software was developed under funding from the U.S. Department of Energy and the U.S. Government
!!! consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting
!!! on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public,
!!! prepare derivative works, and perform publicly and display publicly, and to permit other to do so.
!!!
!!! *** License Agreement ***
!!! ìCrunchFlowî, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory)
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

SUBROUTINE velocalcNS(nx,ny,nz,dtyr)
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

REAL(DP)                                                       :: Re
REAL(DP)                                                       :: Co
REAL(DP)                                                       :: qmax
REAL(DP)                                                       :: dt
REAL(DP), INTENT(IN)                                           :: dtyr

!  Internal variables and arrays

INTEGER(I4B)                                                          :: jx
INTEGER(I4B)                                                          :: jy
INTEGER(I4B)                                                          :: jz

!  ****** PARAMETERS  ****************************

REAL(DP), PARAMETER                                                   :: visc=0.000001d0
REAL(DP), PARAMETER                                                   :: ct=9.135E-10
REAL(DP), PARAMETER                                                   :: big=1.0d0
REAL(DP), PARAMETER                                                   :: zero=0.0d0

CHARACTER (LEN=1)                                                     :: Coordinate


dt = dtyr * 365 * 86400


!   calculate Navier Stokes flux (Only support 2D XY for now)
jz = 1
! --- Velocity on interior interfaces
DO jy = 0,ny-1
    DO jx = 0,nx
        Coordinate = 'X'
        IF (permx(jx,jy,jz) == 0.0d0) THEN
            qx(jx,jy,jz) = 0.0d0
        ELSE
            qx(jx,jy,jz) = us(jx,jy,jz) - (2.0d0*dt/(ro(jx,jy,jz)*(dxx(jx)+dxx(jx+1))))*(pres(jx+1,jy+1,jz)-pres(jx,jy+1,jz))
        END IF
    END DO
END DO
DO jy = 0,ny
    DO jx = 0,nx-1
        Coordinate = 'Y'
        IF (permy(jx,jy,jz) == 0.0d0) THEN
            qy(jx,jy,jz) = 0.0d0
        ELSE
            qy(jx,jy,jz) = vs(jx,jy,jz) - (2.0d0*dt/(ro(jx,jy,jz)*(dyy(jy)+dyy(jy+1))))*(pres(jx+1,jy+1,jz)-pres(jx+1,jy,jz))
        END IF
    END DO
END DO

!   Check Reynolds Number
! qmax = qy(60,30,1)
! Re = qmax * 20.0d0 * 2.5E-4 / visc
! Co = qmax * dt / 2.5E-4
!
! WRITE(*,*) ' Reynolds number = ', Re
! WRITE(*,*) ' Courant number = ', Co

!DO jy = 0,2!ny
!    DO jx = 0,nx-1
!        WRITE(*,*) '>(jx,jy),dt,ro,vs,dp,qy = ',jx,jy,dt,ro(jx,jy,jz),vs(jx,jy,1),(2.0d0*dt/(ro(jx,jy,jz)*(dyy(jy)+dyy(jy+1))))*(pres(jx+1,jy+1,jz)-pres(jx+1,jy,jz)),qy(jx,jy,1)
!    END DO
!END DO
!WRITE(*,*) '>>>>>>>>>>'
!DO jx = 0,nx
!    DO jy = 0,ny-1
!        WRITE(*,*) '>>> (jx,jy), dt, us, p1, p2, qx = ',jx,jy,dt,us(jx,jy,1),pres(jx+1,jy,1),pres(jx,jy,1),qx(jx,jy,1)
!    END DO
!END DO
!WRITE(*,*) '>>>>>>>>>>'



DO jy = 0,ny-1
    DO jx = 0,nx
        qx(jx,jy,jz) = qx(jx,jy,jz) * 365 * 86400
    END DO
END DO
DO jy = 0,ny
    DO jx = 0,nx-1
        qy(jx,jy,jz) = qy(jx,jy,jz) * 365 * 86400
    END DO
END DO


RETURN
END SUBROUTINE velocalcNS
