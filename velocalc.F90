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

SUBROUTINE velocalc(nx,ny,nz)
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

REAL(DP)                                                              :: vv

REAL(DP)                                                              :: RoAveLeft
REAL(DP)                                                              :: RoAveRight

!  ****** PARAMETERS  ****************************

REAL(DP), PARAMETER                                                   :: visc=0.001d0
REAL(DP), PARAMETER                                                   :: ct=9.135E-10
REAL(DP), PARAMETER                                                   :: big=1.0d0
REAL(DP), PARAMETER                                                   :: zero=0.0d0
REAL(DP), PARAMETER                                                   :: grav=9.806d0
REAL(DP)                                                              :: term1
REAL(DP)                                                              :: term2

CHARACTER (LEN=1)                                                     :: Coordinate

vv = secyr/visc  ! Convert to m/yr

!   calculate darcy fluxes

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      IF (jx /= nx) THEN
        Coordinate = 'X'  
        call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft) 
        qx(jx,jy,jz)= -2.0d0*vv*(pres(jx+1,jy,jz)-pres(jx,jy,jz))*harx(jx,jy,jz) /(dxx(jx)+dxx(jx+1))         &
                     + SignGravity*vv*harx(jx,jy,jz)*RoAveRight*grav*COSD(x_angle)
      END IF
      IF (jy /= ny) THEN
        Coordinate = 'Y'  
        call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft) 
        term1 = -2.0d0*vv*(pres(jx,jy+1,jz)-pres(jx,jy,jz))*hary(jx,jy,jz) /(dyy(jy)+dyy(jy+1))
        term2 = SignGravity*vv*hary(jx,jy,jz)*RoAveRight*grav*COSD(y_angle)
        qy(jx,jy,jz)= -2.0d0*vv*(pres(jx,jy+1,jz)-pres(jx,jy,jz))*hary(jx,jy,jz) /(dyy(jy)+dyy(jy+1))     &
                     + SignGravity*vv*hary(jx,jy,jz)*RoAveRight*grav*COSD(y_angle)
      END IF
      IF (jz /= nz) THEN
        Coordinate = 'Z'  
        call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
        qz(jx,jy,jz)= -2.0d0*vv*(pres(jx,jy,jz+1)-pres(jx,jy,jz))*harz(jx,jy,jz) /(dzz(jx,jy,jz)+dzz(jx,jy,jz+1))     &
                     + SignGravity*vv*harz(jx,jy,jz)*RoAveRight*grav*COSD(z_angle)
      END IF
    
    END DO
  END DO
END DO

DO jz = 1,nz
  DO jx = 1,nx
    IF (activecellPressure(jx,0,jz) == 0) THEN
      term1 = -2.0d0*vv*(pres(jx,1,jz)-pres(jx,0,jz))*hary(jx,0,jz) /(dyy(1))
      term2 = SignGravity*vv*hary(jx,0,jz)*ro(jx,1,jz)*grav*COSD(y_angle)
      qy(jx,0,jz)= term1 + term2    
    ELSE 
      qy(jx,0,jz) = 0.0d0
    END IF
    IF (activecellPressure(jx,ny+1,jz) == 0) THEN
      term1 = -2.0d0*vv*(pres(jx,ny+1,jz)-pres(jx,ny,jz))*hary(jx,ny,jz) /(dyy(ny))
      term2 = SignGravity*vv*hary(jx,ny,jz)*ro(jx,ny,jz)*grav*COSD(y_angle)
      qy(jx,ny,jz)= term1 + term2    
    ELSE 
      qy(jx,ny,jz) = 0.0d0
    END IF
  END DO
END DO

DO jz = 1,nz
  DO jy = 1,ny
    IF (activecellPressure(0,jy,jz) == 0) THEN
!!!      term1 = -2.0d0*vv*(pres(1,jy,jz)-pres(0,jy,jz))*harx(1,jy,jz) /(dxx(1))
      term1 = -2.0d0*vv*(pres(1,jy,jz)-pres(0,jy,jz))*harx(0,jy,jz) /(dxx(1))
      term2 = SignGravity*vv*harx(0,jy,jz)*ro(1,jy,jz)*grav*COSD(x_angle)
      qx(0,jy,jz)= term1 + term2  
    ELSE 
      qx(0,jy,jz) = 0.0d0
    END IF
      
    IF (activecellPressure(nx+1,jy,jz) == 0) THEN
      term1 = -2.0d0*vv*(pres(nx+1,jy,jz)-pres(nx,jy,jz))*harx(nx,jy,jz) /(dxx(nx))
      term2 = SignGravity*vv*harx(nx,jy,jz)*ro(nx,jy,jz)*grav*COSD(x_angle)
      qx(nx,jy,jz)= term1 + term2    
    ELSE 
      qx(nx,jy,jz) = 0.0d0
    END IF
  END DO
END DO

DO jy = 1,ny
  DO jx = 1,nx
    IF (activecellPressure(jx,jy,0) == 0) THEN
      term1 = -2.0d0*vv*(pres(jx,jy,1)-pres(jx,jy,0))*harz(jx,jy,0) /(dzz(jx,jy,1))
      term2 = SignGravity*vv*harz(jx,jy,1)*ro(jx,jy,1)*grav*COSD(z_angle)
      qz(jx,jy,0)= term1 + term2    
    ELSE 
      qz(jx,jy,0) = 0.0d0
    END IF
    IF (activecellPressure(jx,jy,nz+1) == 0) THEN
      term1 = -2.0d0*vv*(pres(jx,jy,nz+1)-pres(jx,jy,nz))*harz(jx,jy,nz) /(dzz(jx,jy,nz))
      term2 = SignGravity*vv*harz(jx,jy,nz)*ro(jx,jy,nz)*grav*COSD(z_angle)
      qz(jx,jy,nz)= term1 + term2    
    ELSE 
      qz(jx,jy,nz) = 0.0d0
    END IF
  END DO
END DO

RETURN
END SUBROUTINE velocalc
