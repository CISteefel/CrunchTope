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

SUBROUTINE coeffDiffuse(nx,ny,nz)
USE crunchtype
USE medium
USE transport
USE concentration
USE temperature
USE params
USE CrunchFunctions

IMPLICIT NONE

!  External arrays and variables

INTEGER(I4B), INTENT(IN)                      :: nx
INTEGER(I4B), INTENT(IN)                      :: ny
INTEGER(I4B), INTENT(IN)                      :: nz

! Internal arrays and variables

REAL(DP)                                      :: sum
REAL(DP)                                      :: porp
REAL(DP)                                      :: pore
REAL(DP)                                      :: porw
REAL(DP)                                      :: satw
REAL(DP)                                      :: satp
REAL(DP)                                      :: sate
REAL(DP)                                      :: satn
REAL(DP)                                      :: sats
REAL(DP)                                      :: dumw
REAL(DP)                                      :: dume
REAL(DP)                                      :: dumpx
REAL(DP)                                      :: dumn
REAL(DP)                                      :: dums
REAL(DP)                                      :: dumpy
REAL(DP)                                      :: dxe
REAL(DP)                                      :: dxw
REAL(DP)                                      :: dys
REAL(DP)                                      :: dyn

REAL(DP)                                      :: satup
REAL(DP)                                      :: satdn
REAL(DP)                                      :: dumup
REAL(DP)                                      :: dumdn
REAL(DP)                                      :: dumpz
REAL(DP)                                      :: dzdn
REAL(DP)                                      :: dzup

REAL(DP)                                      :: dumsum
REAL(DP)                                      :: avgro
REAL(DP)                                      :: ae
REAL(DP)                                      :: aw
REAL(DP)                                      :: an
REAL(DP)                                      :: as
REAL(DP)                                      :: aup
REAL(DP)                                      :: adn
REAL(DP)                                      :: tk
REAL(DP)                                      :: quirk
REAL(DP)                                      :: porn
REAL(DP)                                      :: pors
REAL(DP)                                      :: porup
REAL(DP)                                      :: pordn
REAL(DP)                                      :: dharm
REAL(DP)                                      :: dspe
REAL(DP)                                      :: dspw
REAL(DP)                                      :: de
REAL(DP)                                      :: dw
REAL(DP)                                      :: ds
REAL(DP)                                      :: dn
REAL(DP)                                      :: dsps
REAL(DP)                                      :: dspn
REAL(DP)                                      :: apx
REAL(DP)                                      :: apy

REAL(DP)                                      :: ddn
REAL(DP)                                      :: dup
REAL(DP)                                      :: dspdn
REAL(DP)                                      :: dspup
REAL(DP)                                      :: apz
REAL(DP)                                      :: tort

REAL(DP)                                      :: PorPow
REAL(DP)                                      :: SatPow

INTEGER(I4B)                                  :: jx
INTEGER(I4B)                                  :: jy
INTEGER(I4B)                                  :: jz

SatPow = 10.0d0/3.0d0
PorPow = 4.0d0/3.0d0

IF (idiffus == 0) THEN
  d_25 = dzero
ELSE
  d_25 = dcoeff
END IF

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      tk = 273.15 + t(jx,jy,jz)
      IF (idiffus == 0) THEN
        dstar(jx,jy,jz) = dzero*EXP((activation/rgas)*(tk25 - 1.0/tk))/formation
      ELSE
        dstar(jx,jy,jz) = dcoeff/formation
      END IF
    END DO
  END DO
END DO

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
    
      porp = por(jx,jy,jz)
      satp = satliq(jx,jy,jz)
    
      IF (nx == 1) GO TO 100
    
    IF (jx == 1) THEN
      
      dxe = 0.5d0*(dxx(jx)+dxx(jx+1))
      dxw = 0.5d0*dxx(1)
      pore = por(jx+1,jy,jz)
      porw = por(jx,jy,jz)
      sate = satliq(jx+1,jy,jz)
      satw = satliq(jx,jy,jz)
      IF (UseThresholdPorosity) THEN
        IF (pore > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dume = ro(jx+1,jy,jz)*dstar(jx+1,jy,jz)*sate*pore*tort
        IF (porp > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpx = ro(jx,jy,jz)*dstar(jx,jy,jz)*sate*porp*tort
        dumw = dumpx
      ELSE IF (MillingtonQuirk) THEN
        dume = ro(jx+1,jy,jz)*(sate)**(SatPow)*(pore)**(PorPow)*dstar(jx+1,jy,jz)
        dumpx = ro(jx,jy,jz)*(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)
        dumw = dumpx
      ELSE
        dume = ro(jx+1,jy,jz)*sate*pore*dstar(jx+1,jy,jz)*tortuosity(jx+1,jy,jz)
        dumpx = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*tortuosity(jx,jy,jz)
        dumw = ro(jx-1,jy,jz)*sate*pore*dstar(jx-1,jy,jz)*tortuosity(jx-1,jy,jz)
!!!        dumw = dumpx
      END IF
    ELSE IF (jx == nx) THEN
      dxw = 0.5d0*(dxx(jx)+dxx(jx-1))
      dxe = 0.5d0*dxx(nx)
      pore = por(jx,jy,jz)
      porw = por(jx-1,jy,jz)
      sate = satliq(jx,jy,jz)
      satw = satliq(jx-1,jy,jz)
      IF (UseThresholdPorosity) THEN
        IF (porw > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumw = ro(jx-1,jy,jz)*dstar(jx-1,jy,jz)*satw*porw*tort
        IF (porp > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpx = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort  
        dume = dumpx 
      ELSE IF (MillingtonQuirk) THEN
        dumw = ro(jx-1,jy,jz)*(satw)**(SatPow)*(porw)**(PorPow)*dstar(jx-1,jy,jz)
        dumpx = ro(jx,jy,jz)*(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)
        dume = dumpx
      ELSE
        dumw = ro(jx-1,jy,jz)*satw*porw*dstar(jx-1,jy,jz)*tortuosity(jx-1,jy,jz)
        dumpx = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*tortuosity(jx,jy,jz)
        dume = ro(jx+1,jy,jz)*satw*porw*dstar(jx+1,jy,jz)*tortuosity(jx+1,jy,jz)
!!!        dume = dumpx
      END IF
    ELSE
      dxe = 0.5d0*(dxx(jx)+dxx(jx+1))
      dxw = 0.5d0*(dxx(jx)+dxx(jx-1))
      pore = por(jx+1,jy,jz)
      porw = por(jx-1,jy,jz)
      sate = satliq(jx+1,jy,jz)
      satw = satliq(jx-1,jy,jz)
      IF (UseThresholdPorosity) THEN
        IF (pore > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dume = ro(jx+1,jy,jz)*dstar(jx+1,jy,jz)*sate*pore*tort
        IF (porw > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumw = ro(jx-1,jy,jz)*dstar(jx-1,jy,jz)*sate*porw*tort
        IF (porp > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpx = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort
      ELSE IF (MillingtonQuirk) THEN
        dume = ro(jx+1,jy,jz)*(sate)**(SatPow)*(pore)**(PorPow)*dstar(jx+1,jy,jz)
        dumpx = ro(jx,jy,jz)*(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)
        dumw = ro(jx-1,jy,jz)*(satw)**(SatPow)*(porw)**(PorPow)*dstar(jx-1,jy,jz)
      ELSE
        dume = ro(jx+1,jy,jz)*sate*pore*dstar(jx+1,jy,jz)*tortuosity(jx+1,jy,jz)
        dumpx = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*tortuosity(jx,jy,jz)
        dumw = ro(jx-1,jy,jz)*satw*porw*dstar(jx-1,jy,jz)*tortuosity(jx-1,jy,jz)
      END IF
    END IF
    
    100     CONTINUE
    IF (ny == 1) GO TO 200
    
    IF (jy == 1) THEN
      dyn = 0.5d0*(dyy(jy)+dyy(jy+1))
      dys = 0.5d0*dyy(1)
      porn = por(jx,jy+1,jz)
      pors = por(jx,jy,jz)
      satn = satliq(jx,jy+1,jz)
      sats = satliq(jx,jy,jz)
      IF (UseThresholdPorosity) THEN
        IF (porn > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumn = ro(jx,jy+1,jz)*dstar(jx,jy+1,jz)*satn*porn*tort*anisotropyY
        IF (por(jx,jy,jz) > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpy = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort*anisotropyY
        dums = dumpy
      ELSE IF (MillingtonQuirk) THEN
        dumn = ro(jx,jy+1,jz)*(satn)**(SatPow)*(por(jx,jy+1,jz))**(PorPow)*dstar(jx,jy+1,jz)*anisotropyY
        dumpy = ro(jx,jy,jz)*(satp)**(SatPow)*(por(jx,jy,jz))**(PorPow)*dstar(jx,jy,jz)*anisotropyY
        dums = dumpy
      ELSE
        dumn = ro(jx,jy+1,jz)*satn*porn*dstar(jx,jy+1,jz)*anisotropyY*tortuosity(jx,jy+1,jz)
        dumpy = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*anisotropyY*tortuosity(jx,jy,jz)
        dums = ro(jx,jy-1,jz)*satn*porn*dstar(jx,jy-1,jz)*anisotropyY*tortuosity(jx,jy-1,jz)
!!!        dums = dumpy
      END IF
    ELSE IF (jy == ny) THEN
      dys = 0.5d0*(dyy(jy)+dyy(jy-1))
      dyn = 0.5d0*dyy(ny)
      porn = por(jx,jy,jz)
      pors = por(jx,jy-1,jz)
      satn = satliq(jx,jy,jz)
      sats = satliq(jx,jy-1,jz)
      IF (UseThresholdPorosity) THEN
        IF (pors > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dums = ro(jx,jy-1,jz)*dstar(jx,jy-1,jz)*sats*pors*tort*anisotropyY
        IF (por(jx,jy,jz) > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpy = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort*anisotropyY
        dumn = dumpy
      ELSE IF (MillingtonQuirk) THEN
        dums = ro(jx,jy-1,jz)*(sats)**(SatPow)*(pors)**(PorPow)*dstar(jx,jy-1,jz)*anisotropyY
        dumpy = ro(jx,jy,jz)*(satp)**(SatPow)*(por(jx,jy,jz))**(PorPow)*dstar(jx,jy,jz)*anisotropyY
        dumn = dumpy
      ELSE
        dums = ro(jx,jy-1,jz)*sats*pors*dstar(jx,jy-1,jz)*anisotropyY*tortuosity(jx,jy-1,jz)
        dumpy = ro(jx,jy,jz)*satp*porP*dstar(jx,jy,jz)*anisotropyY*tortuosity(jx,jy,jz)
        dumn = ro(jx,jy+1,jz)*sats*pors*dstar(jx,jy+1,jz)*anisotropyY*tortuosity(jx,jy+1,jz)
!!!        dumn = dumpy
      END IF
    ELSE
      dyn = 0.5d0*(dyy(jy)+dyy(jy+1))
      dys = 0.5d0*(dyy(jy)+dyy(jy-1))
      porn = por(jx,jy+1,jz)
      pors = por(jx,jy-1,jz)
      satn = satliq(jx,jy+1,jz)
      sats = satliq(jx,jy-1,jz)

      IF (UseThresholdPorosity) THEN
        IF (porn > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumn = ro(jx,jy+1,jz)*dstar(jx,jy+1,jz)*satn*porn*tort*anisotropyY
        IF (pors > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dums = ro(jx,jy-1,jz)*dstar(jx,jy-1,jz)*sats*pors*tort*anisotropyY
        IF (por(jx,jy,jz) > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpy = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort*anisotropyY
      ELSE IF (MillingtonQuirk) THEN
        dumn = ro(jx,jy+1,jz)*(satn)**(SatPow)*(porn)**(PorPow)*dstar(jx,jy+1,jz)*anisotropyY
        dums = ro(jx,jy-1,jz)*(sats)**(SatPow)*(pors)**(PorPow)*dstar(jx,jy-1,jz)*anisotropyY
        dumpy = ro(jx,jy,jz)*(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)*anisotropyY
      ELSE
        dumn = ro(jx,jy+1,jz)*satn*porn*dstar(jx,jy+1,jz)*anisotropyY*tortuosity(jx,jy+1,jz)
        dums = ro(jx,jy-1,jz)*sats*pors*dstar(jx,jy-1,jz)*anisotropyY*tortuosity(jx,jy-1,jz)
        dumpy = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*anisotropyY*tortuosity(jx,jy,jz)
      END IF
    END IF
    
    200     CONTINUE

    IF (nz == 1) GO TO 250
    
    IF (jz == 1) THEN
      
      dzup= 0.5*(dzz(jx,jy,jz)+dzz(jx,jy,jz+1))
      dzdn= 0.5*dzz(jx,jy,jz)
      porup = por(jx,jy,jz+1)
      pordn = por(jx,jy,jz)
      satup = satliq(jx,jy,jz+1)
      satdn = satliq(jx,jy,jz)

      IF (UseThresholdPorosity) THEN
        IF (porup > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumup = ro(jx,jy,jz+1)*dstar(jx,jy,jz+1)*satup*porup*tort*anisotropyZ(jx,jy,jz+1)
        IF (por(jx,jy,jz) > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpz = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort*anisotropyZ(jx,jy,jz)
      ELSE IF (MillingtonQuirk) THEN
        dumup = ro(jx,jy,jz+1)*(satup)**(SatPow)*(porup)**(PorPow)*dstar(jx,jy,jz+1)*anisotropyZ(jx,jy,jz+1)
        dumpz = ro(jx,jy,jz)*(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)*anisotropyZ(jx,jy,jz)
      ELSE
        dumup = ro(jx,jy,jz+1)*satup*porup*dstar(jx,jy,jz+1)*anisotropyZ(jx,jy,jz+1)*tortuosity(jx,jy,jz+1)
        dumpz = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*anisotropyZ(jx,jy,jz)*tortuosity(jx,jy,jz)
        dumdn = ro(jx,jy,jz-1)*satup*porup*dstar(jx,jy,jz-1)*anisotropyZ(jx,jy,jz-1)*tortuosity(jx,jy,jz-1)
      END IF

!!!      dumdn = dumpz
    ELSE IF (jz == nz) THEN
      dzdn = 0.5*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))
      dzup = 0.5*dzz(jx,jy,jz)
      porup = por(jx,jy,jz)
      pordn = por(jx,jy,jz-1)
      satup = satliq(jx,jy,jz)
      satdn = satliq(jx,jy,jz-1)

      IF (UseThresholdPorosity) THEN
        IF (pordn > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumdn = ro(jx,jy,jz-1)*dstar(jx,jy,jz-1)*satdn*pordn*tort*anisotropyZ(jx,jy,jz-1)
        IF (por(jx,jy,jz) > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpz = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort*anisotropyZ(jx,jy,jz)
      ELSE IF (MillingtonQuirk) THEN
        dumdn = ro(jx,jy,jz-1)*(satdn)**(SatPow)*(pordn)**(PorPow)*dstar(jx,jy,jz-1)*anisotropyZ(jx,jy,jz-1)
        dumpz = ro(jx,jy,jz)*(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)*anisotropyZ(jx,jy,jz)
      ELSE
        dumdn = ro(jx,jy,jz-1)*satdn*pordn*dstar(jx,jy,jz-1)*anisotropyZ(jx,jy,jz-1)*tortuosity(jx,jy,jz-1)
        dumpz = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*anisotropyZ(jx,jy,jz)*tortuosity(jx,jy,jz)
        dumup = ro(jx,jy,jz+1)*satdn*pordn*dstar(jx,jy,jz+1)*anisotropyZ(jx,jy,jz+1)*tortuosity(jx,jy,jz+1)
      END IF

!!!      dumup = dumpz
    ELSE
      dzup = 0.5d0*(dzz(jx,jy,jz)+dzz(jx,jy,jz+1))
      dzdn = 0.5d0*(dzz(jx,jy,jz)+dzz(jx,jy,jz-1))
      porup = por(jx,jy,jz+1)
      pordn = por(jx,jy,jz-1)
      satup = satliq(jx,jy,jz+1)
      satdn = satliq(jx,jy,jz-1)

      IF (UseThresholdPorosity) THEN
        IF (porup > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumup = ro(jx,jy,jz+1)*dstar(jx,jy,jz+1)*satup*porup*tort*anisotropyZ(jx,jy,jz+1)
        IF (pordn > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumdn = ro(jx,jy,jz-1)*dstar(jx,jy,jz-1)*satdn*pordn*tort*anisotropyZ(jx,jy,jz-1)
        IF (por(jx,jy,jz) > ThresholdPorosity) THEN
          tort = TortuosityAboveThreshold
        ELSE
          tort = TortuosityBelowThreshold
        END IF
        dumpz = ro(jx,jy,jz)*dstar(jx,jy,jz)*satp*porp*tort*anisotropyZ(jx,jy,jz)
      ELSE IF (MillingtonQuirk) THEN
        dumup = ro(jx,jy,jz+1)*(satup)**(SatPow)*(porup)**(PorPow)*dstar(jx,jy,jz+1)*anisotropyZ(jx,jy,jz+1)
        dumdn = ro(jx,jy,jz-1)*(satdn)**(SatPow)*(pordn)**(PorPow)*dstar(jx,jy,jz-1)*anisotropyZ(jx,jy,jz-1)
        dumpz = ro(jx,jy,jz)*(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)*anisotropyZ(jx,jy,jz)
      ELSE
        dumup = ro(jx,jy,jz+1)*satup*porup*dstar(jx,jy,jz+1)*anisotropyZ(jx,jy,jz+1)*tortuosity(jx,jy,jz+1)
        dumdn = ro(jx,jy,jz-1)*satdn*pordn*dstar(jx,jy,jz-1)*anisotropyZ(jx,jy,jz-1)*tortuosity(jx,jy,jz-1)
        dumpz = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*anisotropyZ(jx,jy,jz)*tortuosity(jx,jy,jz)
      END IF
    END IF
    
      250     CONTINUE
      
      IF (nx == 1) GO TO 300
    
      IF (jx == 1) THEN
      
        avgro = 0.5*( ro(jx+1,jy,jz) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dume,dumpx)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dume,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dume,dumpx)
      ELSE
        dharm = GeometricMean(dume,dumpx)
      END IF
        dspe = avgro*dspx(jx,jy,jz) + dharm
        de = dspe*dyy(jy)*dzz(jx,jy,jz)/dxe

        ae = de
        aw = 0.0                     !  Assume no diffusion across boundary
        apx = de                     !  Assume no diffusion across boundary
      
      ELSE IF (jx == nx) THEN
      
        avgro = 0.5*( ro(jx-1,jy,jz) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dumw,dumpx)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dume,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumw,dumpx)
      ELSE
        dharm = GeometricMean(dumw,dumpx)
      END IF
        dspw = avgro*dspx(jx-1,jy,jz) + dharm
        dw = dspw*dyy(jy)*dzz(jx,jy,jz)/dxw

        aw = dw
        ae = 0.0                     !  Assume no diffusion across boundary
        apx = dw                     !  Assume no diffusion across boundary
      
      ELSE
        
        avgro = 0.5*( ro(jx+1,jy,jz) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dume,dumpx)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dume,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dume,dumpx)
      ELSE
        dharm = GeometricMean(dume,dumpx)
      END IF
        dspe = avgro*dspx(jx,jy,jz) + dharm
        de = dspe*dyy(jy)*dzz(jx,jy,jz)/dxe
      
        avgro = 0.5*( ro(jx-1,jy,jz) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dumw,dumpx)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumw,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumw,dumpx)
      ELSE
        dharm = GeometricMean(dumw,dumpx)
      END IF
        dspw = avgro*dspx(jx-1,jy,jz) + dharm
        dw = dspw*dyy(jy)*dzz(jx,jy,jz)/dxw

        aw = dw
        ae = de
        apx = dw + de 
      
      END IF
    
      300     CONTINUE
      IF (ny == 1) GO TO 400
    
      IF (jy == 1) THEN
      
        avgro = 0.5*( ro(jx,jy+1,jz) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dumn,dumpy)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumn,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumn,dumpy)
      ELSE
        dharm = GeometricMean(dumn,dumpy)
      END IF
        dspn = avgro*dspy(jx,jy,jz) + dharm
        dn = dspn*dxx(jx)*dzz(jx,jy,jz)/dyn

        an = dn
        as = 0.0
        apy = dn

      ELSE IF (jy == ny) THEN
         
        avgro = 0.5*( ro(jx,jy-1,jz) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dums,dumpy)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dums,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dums,dumpy)
      ELSE
        dharm = GeometricMean(dums,dumpy)
      END IF
        dsps = avgro*dspy(jx,jy-1,jz) + dharm
        ds = dsps*dxx(jx)*dzz(jx,jy,jz)/dys

        as = ds   
        an = 0.00
        apy = ds 
      
      ELSE
      
        avgro = 0.5*( ro(jx,jy+1,jz) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dumn,dumpy)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumn,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumn,dumpy)
      ELSE
        dharm = GeometricMean(dumn,dumpy)
      END IF
        dspn = avgro*dspy(jx,jy,jz) + dharm
        dn = dspn*dxx(jx)*dzz(jx,jy,jz)/dyn
      
        avgro = 0.5*( ro(jx,jy-1,jz) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dums,dumpy)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dums,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dums,dumpy)
      ELSE
        dharm = GeometricMean(dums,dumpy)
      END IF
        dsps = avgro*dspy(jx,jy-1,jz) + dharm
        ds = dsps*dxx(jx)*dzz(jx,jy,jz)/dys

        an = dn
        as = ds
        apy = ds + dn 
      
      END IF
    
      400     CONTINUE

      IF (nz == 1) GO TO 500
    
      IF (jz == 1) THEN
      
        avgro = 0.5*( ro(jx,jy,jz+1) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dumup,dumpz)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumup,dumpz)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumup,dumpz)
      ELSE
        dharm = GeometricMean(dumup,dumpz)
      END IF
        dspup = avgro*dspz(jx,jy,jz) + dharm
        dup = dspup*dxx(jx)*dyy(jy)/dzup
        aup = dup
      
        adn = 0.0
      
        apz = dup

      
      ELSE IF (jz == nz) THEN
      
        aup= 0.00
      
        avgro = 0.5*( ro(jx,jy,jz-1) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dumdn,dumpz)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumdn,dumpz)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumdn,dumpz)
      ELSE
        dharm = GeometricMean(dumdn,dumpz)
      END IF
        dspdn = avgro*dspz(jx,jy,jz-1) + dharm
        ddn = dspdn*dxx(jx)*dyy(jy)/dzdn
        adn = ddn
      
        apz = ddn 

      ELSE
      
        avgro = 0.5*( ro(jx,jy,jz+1) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dumup,dumpz)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumup,dumpz)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumup,dumpz)
      ELSE
        dharm = GeometricMean(dumup,dumpz)
      END IF
        dspup = avgro*dspz(jx,jy,jz) + dharm
        dup = dspup*dxx(jx)*dyy(jy)/dzup
        aup = dup
      
        avgro = 0.5*( ro(jx,jy,jz-1) + ro(jx,jy,jz) )
!!        dharm = ArithmeticMean(dumdn,dumpz)
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumdn,dumpz)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumdn,dumpz)
      ELSE
        dharm = GeometricMean(dumdn,dumpz)
      END IF
        dspdn = avgro*dspz(jx,jy,jz-1) + dharm
        ddn = dspdn*dxx(jx)*dyy(jy)/dzdn
        adn = ddn
      
        apz = ddn + dup 
      
      END IF
    
    
      500     CONTINUE

      IF (nx == 1) THEN
        CONTINUE
      ELSE      
        aDD(jx,jy,jz) = -aw
        cDD(jx,jy,jz) = -ae
        bDD(jx,jy,jz) = apx     
      END IF
    
      IF (ny == 1) THEN
        CONTINUE
      ELSE    
        dDD(jx,jy,jz) = -an
        fDD(jx,jy,jz) = -as
        eDD(jx,jy,jz) = apy   
      END IF

      IF (nz == 1) THEN
        CONTINUE
      ELSE    
        gDD(jx,jy,jz) = -aup
        iDD(jx,jy,jz) = -adn
        hDD(jx,jy,jz) = apz   
      END IF

      dxy(jx,jy,jz) = dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
    
    END DO
  END DO
END DO


RETURN
END SUBROUTINE coeffDiffuse
