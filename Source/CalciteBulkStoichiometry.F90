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
    
subroutine CalciteBulkStoichiometry(nx,ny,nz,dt)
USE crunchtype
USE mineral
USE concentration

IMPLICIT NONE

!!  External arrays and variables

INTEGER(I4B), INTENT(IN)                                   :: nx
INTEGER(I4B), INTENT(IN)                                   :: ny
INTEGER(I4B), INTENT(IN)                                   :: nz

REAL(DP), INTENT(IN)                                       :: dt

!!  Internal arrays and variables

REAL(DP)                                                   :: volcheck
REAL(DP)                                                   :: MineralSum
REAL(DP)                                                   :: rdecay
REAL(DP)                                                   :: tot

INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      volcheck = volmol(kUCalcite)*dppt(kUCalcite,jx,jy,jz)*dt
      MineralSum = volcheck + volfx(kUCalcite,jx,jy,jz)

      muUranium234Bulk(jx,jy,jz) = muUranium234Bulk(jx,jy,jz)*volfx(kUCalcite,jx,jy,jz)/MineralSum + volcheck*muUranium234(jx,jy,jz)/MineralSum
      muUranium238Bulk(jx,jy,jz) = muUranium238Bulk(jx,jy,jz)*volfx(kUCalcite,jx,jy,jz)/MineralSum + volcheck*muUranium238(jx,jy,jz)/MineralSum
      muCalciumBulk(jx,jy,jz)    = muCalciumBulk(jx,jy,jz)*volfx(kUCalcite,jx,jy,jz)/MineralSum    + volcheck*muCalcium(jx,jy,jz)/MineralSum

      rdecay = -2.88e-06*muUranium234Bulk(jx,jy,jz)*dt  + 1.55e-10*muUranium238Bulk(jx,jy,jz)*dt   !!  234U
      muUranium234Bulk(jx,jy,jz) = muUranium234Bulk(jx,jy,jz) + rdecay

      rdecay = -1.55e-10*muUranium238Bulk(jx,jy,jz)*dt                                                !!  238U
      muUranium238Bulk(jx,jy,jz) = muUranium238Bulk(jx,jy,jz) + rdecay

    END DO
  END DO
END DO


END SUBROUTINE CalciteBulkStoichiometry
