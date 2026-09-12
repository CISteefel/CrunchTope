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

SUBROUTINE harmonic(nx,ny,nz)
USE crunchtype
USE params
USE medium
USE flow

IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                           :: nx
INTEGER(I4B), INTENT(IN)                           :: ny
INTEGER(I4B), INTENT(IN)                           :: nz

!  Internal variables

INTEGER(I4B)                                       :: jx
INTEGER(I4B)                                       :: jy
INTEGER(I4B)                                       :: jz


!   generate average (harmonic mean) conductivities at grid cell faces

DO jz = 1,nz 
  DO jy = 1,ny
    DO jx = 0,nx
      IF (permx(jx,jy,jz) == 0.0 .AND. permx(jx+1,jy,jz) == 0.0) THEN
        harx(jx,jy,jz) = 0.0
      ELSE
        harx(jx,jy,jz) = permx(jx,jy,jz)*permx(jx+1,jy,jz) *(dxx(jx)+dxx(jx+1))/(dxx(jx)*permx(jx+1,jy,jz) + dxx(jx+1)*permx(jx,jy,jz))
        continue
      END IF
    END DO
  END DO
END DO

DO jz = 1,nz 
  DO jy = 0,ny
    DO jx = 1,nx
      IF (permy(jx,jy,jz) == 0.0 .AND. permy(jx,jy+1,jz) == 0.0) THEN
        hary(jx,jy,jz) = 0.0
      ELSE
        hary(jx,jy,jz) = permy(jx,jy,jz)*permy(jx,jy+1,jz) *(dyy(jy)+dyy(jy+1))/(dyy(jy)*permy(jx,jy+1,jz) + dyy(jy+1)*permy(jx,jy,jz))
      END IF
    END DO
  END DO
END DO

DO jz = 0,nz 
  DO jy = 1,ny
    DO jx = 1,nx
      IF (permz(jx,jy,jz) == 0.0 .AND. permz(jx,jy,jz+1) == 0.0) THEN
        harz(jx,jy,jz) = 0.0
      ELSE
        harz(jx,jy,jz) = permz(jx,jy,jz)*permz(jx,jy,jz+1) *(dzz(jx,jy,jz)+dzz(jx,jy,jz+1))  &
          /(dzz(jx,jy,jz)*permz(jx,jy,jz+1) + dzz(jx,jy,jz+1)*permz(jx,jy,jz))
      END IF
    END DO
  END DO
END DO

RETURN
END SUBROUTINE harmonic


