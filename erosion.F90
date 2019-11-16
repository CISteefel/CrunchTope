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

SUBROUTINE erosion(nx,ny,nz)
USE crunchtype
USE params
USE medium
USE transport
USE temperature
 
IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                                   :: nx
INTEGER(I4B), INTENT(IN)                                   :: ny
INTEGER(I4B), INTENT(IN)                                   :: nz

!  Internal variables

INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

REAL(DP)                                                   :: zero
REAL(DP)                                      :: porp
REAL(DP)                                      :: pore
REAL(DP)                                      :: porw
REAL(DP)                                      :: dxe
REAL(DP)                                      :: dxw
REAL(DP)                                      :: dys
REAL(DP)                                      :: dyn
REAL(DP)                                      :: porharm
REAL(DP)                                      :: porsum
REAL(DP)                                      :: avgro
REAL(DP)                                      :: fe
REAL(DP)                                      :: ae
REAL(DP)                                      :: fw
REAL(DP)                                      :: aw
REAL(DP)                                      :: fn
REAL(DP)                                      :: an
REAL(DP)                                      :: fs
REAL(DP)                                      :: as
REAL(DP)                                      :: porn
REAL(DP)                                      :: pors
REAL(DP)                                      :: apx
REAL(DP)                                      :: apy

zero = 0.0D0

jz = 1
DO jy = 1,ny
  DO jx = 1,nx
    
    IF (jx == 1) THEN
      
      fe = dyy(jy)*SolidBuryX(jx)                 !! Now working with sites per bulk porous medium
      ae = DMAX1(-fe,zero)
      
      fw = dyy(jy)*SolidBuryX(jx-1)
      aw = DMAX1(fw,zero)
      
      apx = DMAX1(-fw,zero) + DMAX1(fe,zero)
      
    ELSE IF (jx == nx) THEN
      
      fe = dyy(jy)*SolidBuryX(jx)                 !! Now working with sites per bulk porous medium
      ae = DMAX1(-fe,zero)
      
      fw = dyy(jy)*SolidBuryX(jx-1)
      aw = DMAX1(fw,zero)
      
      apx = DMAX1(-fw,zero) + DMAX1(fe,zero)
      
    ELSE
      
      fe = dyy(jy)*SolidBuryX(jx)                 !! Now working with sites per bulk porous medium
      ae = DMAX1(-fe,zero)
      
      fw = dyy(jy)*SolidBuryX(jx-1)
      aw = DMAX1(fw,zero)
      
      apx = DMAX1(-fw,zero) + DMAX1(fe,zero)
      
    END IF
    
    300     CONTINUE
    IF (ny == 1) GO TO 400
    
    IF (jy == 1) THEN
      
      fn = dxx(jx)*SolidBuryY(jy)
      an = DMAX1(-fn,zero)
      
      fs = dxx(jx)*SolidBuryY(jy)
      as = DMAX1(fs,zero)
      
      apy = DMAX1(-fs,zero) + DMAX1(fn,zero)
      
    ELSE IF (jy == ny) THEN
      
      fn = dxx(jx)*SolidBuryY(jy)
      an = DMAX1(-fn,zero)
      
      fs = dxx(jx)*SolidBuryY(jy)
      as = DMAX1(fs,zero)
      
      apy = DMAX1(-fs,zero) + DMAX1(fn,zero)
      
    ELSE
      
      fn = dxx(jx)*SolidBuryY(jy)
      an = DMAX1(-fn,zero)
      
      fs = dxx(jx)*SolidBuryY(jy)
      as = DMAX1(fs,zero)
      
      apy = DMAX1(-fs,zero) + DMAX1(fn,zero)
      
    END IF
    
    
    400     CONTINUE
    IF (nx == 1) THEN
      abu(jx,jy,jz) = zero
      bbu(jx,jy,jz) = zero
      cbu(jx,jy,jz) = zero
    ELSE
      
      abu(jx,jy,jz) = -aw
      cbu(jx,jy,jz) = -ae
      bbu(jx,jy,jz) = apx
      
    END IF
    
    IF (ny == 1) THEN
      dbu(jx,jy,jz) = zero
      fbu(jx,jy,jz) = zero
      ebu(jx,jy,jz) = zero
    ELSE
      
      dbu(jx,jy,jz) = -an
      fbu(jx,jy,jz) = -as
      ebu(jx,jy,jz) = apy
      
    END IF
    
  END DO
END DO


RETURN
END SUBROUTINE erosion
