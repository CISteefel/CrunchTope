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

SUBROUTINE coeff(nx,ny,nz)
USE crunchtype
USE medium
USE transport
USE concentration
USE temperature
USE params
USE CrunchFunctions
USE flow

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
REAL(DP)                                      :: dumsum
REAL(DP)                                      :: avgro
REAL(DP)                                      :: fe
REAL(DP)                                      :: ae
REAL(DP)                                      :: fw
REAL(DP)                                      :: aw
REAL(DP)                                      :: fn
REAL(DP)                                      :: an
REAL(DP)                                      :: fs
REAL(DP)                                      :: as
REAL(DP)                                      :: tk
REAL(DP)                                      :: Uliaq
REAL(DP)                                      :: Quirkaq
REAL(DP)                                      :: porn
REAL(DP)                                      :: pors
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
REAL(DP)                                      :: tort
REAL(DP)                                      :: AreaE
REAL(DP)                                      :: AreaW
REAL(DP)                                      :: AreaS
REAL(DP)                                      :: AreaN
REAL(DP)                                      :: mu_water_20

REAL(DP)                                      :: PorPow
REAL(DP)                                      :: SatPow

INTEGER(I4B)                                  :: jx
INTEGER(I4B)                                  :: jy
INTEGER(I4B)                                  :: jz
INTEGER(I4B)                                  :: j

jz = 1

SatPow = 10.0d0/3.0d0
PorPow = 4.0d0/3.0d0

IF (idiffus == 0) THEN
  d_25 = dzero
ELSE
  d_25 = dcoeff
END IF

!!jz = 1
!!DO jy = 1,ny
!!  DO jx = 1,nx
!!    j = (jy-1)*nx+jx
!!    tk = 273.15d0 + t(jx,jy,jz)
!!    IF (idiffus == 0) THEN
!!      dstar(jx,jy,jz) = dzero*DEXP((activation/rgas)*(tk25 - 1.0d0/tk))/formation
!!    ELSE
!!      dstar(jx,jy,jz) = dcoeff/formation
!!    END IF
!!  END DO
!!END DO

jz = 1
DO jy = 1,ny
  DO jx = 1,nx
    j = (jy-1)*nx+jx
    tk = 273.15 + t(jx,jy,jz)
    IF (idiffus == 0) THEN
      dstar(jx,jy,jz) = dzero*EXP((activation/rgas)*(tk25 - 1.0/tk))/formation
    ELSE
      dstar(jx,jy,jz) = (dcoeff/formation)
    END IF
  END DO
  dstar(0,jy,jz) = dstar(1,jy,jz)
  dstar(nx+1,jy,jz) = dstar(nx,jy,jz)
END DO
DO jx = 1,nx
  dstar(jx,0,jz) = dstar(jx,1,jz)
  dstar(jx,ny+1,jz) = dstar(jx,ny,jz)
END DO

jz = 1
DO jy = 1,ny
  DO jx = 1,nx

    j = (jy-1)*nx+jx
    
    porp = por(jx,jy,jz)
!!!    satp = satliq(jx,jy,jz)
    satp = 0.5*( satliq(jx,jy,jz)+satliqold(jx,jy,jz) )
    
    IF (nx == 1) GO TO 100
    
    IF (jx == 1) THEN
      
      dxe = 0.5d0*(dxx(jx)+dxx(jx+1))
      dxw = 0.5d0*dxx(1)
      pore = por(jx+1,jy,jz)
      porw = por(jx,jy,jz)
      sate = 0.5*( satliq(jx+1,jy,jz)+satliqold(jx+1,jy,jz) )
      satw = 0.5*( satliq(jx,jy,jz)+satliqold(jx,jy,jz) )
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
        dumpx = ro(jx,jy,jz)*(satp)**(SatPow) * (porp)**(PorPow)*dstar(jx,jy,jz)

        dumw = dumpx
      ELSE
        dume = ro(jx+1,jy,jz)*sate*pore*dstar(jx+1,jy,jz)*tortuosity(jx+1,jy,jz)
        dumpx = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*tortuosity(jx,jy,jz)
        dumw = dumpx

      END IF
      
    ELSE IF (jx == nx) THEN
      dxw = 0.5d0*(dxx(jx)+dxx(jx-1))
      dxe = 0.5d0*dxx(nx)
      pore = por(jx,jy,jz)
      porw = por(jx-1,jy,jz)
      sate = 0.5*( satliq(jx,jy,jz)+satliqold(jx,jy,jz) )
      satw = 0.5*( satliq(jx-1,jy,jz)+satliqold(jx-1,jy,jz) )
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
        dume = dumpx
      END IF
      
!!! Hardwired by Lucien
!!!      IF (east_river .and. ny == 1 .and. nz == 1) THEN
!!!        dume = 0
!!!      ENDIF
      
    ELSE
      dxe = 0.5d0*(dxx(jx)+dxx(jx+1))
      dxw = 0.5d0*(dxx(jx)+dxx(jx-1))
      pore = por(jx+1,jy,jz)
      porw = por(jx-1,jy,jz)
      sate = 0.5*( satliq(jx+1,jy,jz)+satliqold(jx+1,jy,jz) )
      satw = 0.5*( satliq(jx-1,jy,jz)+satliqold(jx-1,jy,jz) )
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
!!!      satn = satliq(jx,jy+1,jz)
!!!      sats = satliq(jx,jy,jz)
      satn = 0.5*(satliq(jx,jy+1,jz) + satliqold(jx,jy+1,jz) )
      sats = 0.5*(satliq(jx,jy,jz) + satliqold(jx,jy,jz) )
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
        dums = dumpy
      END IF
    ELSE IF (jy == ny) THEN
      dys = 0.5d0*(dyy(jy)+dyy(jy-1))
      dyn = 0.5d0*dyy(ny)
      porn = por(jx,jy,jz)
      pors = por(jx,jy-1,jz)
!!!      satn = satliq(jx,jy,jz)
!!!      sats = satliq(jx,jy-1,jz)
      satn = 0.5*(satliq(jx,jy,jz) + satliqold(jx,jy,jz) )
      sats = 0.5*(satliq(jx,jy-1,jz) + satliqold(jx,jy-1,jz) )
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
        dumpy = ro(jx,jy,jz) *(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)*anisotropyY
        dumn = dumpy
      ELSE
        dums = ro(jx,jy-1,jz)*sats*pors*dstar(jx,jy-1,jz)*anisotropyY*tortuosity(jx,jy-1,jz)
        dumpy = ro(jx,jy,jz) *satp*porP*dstar(jx,jy,jz)  *anisotropyY*tortuosity(jx,jy,jz)
        dumn = dumpy
      END IF
    ELSE
      
      dyn = 0.5d0*(dyy(jy)+dyy(jy+1))
      dys = 0.5d0*(dyy(jy)+dyy(jy-1))
      porn = por(jx,jy+1,jz)
      pors = por(jx,jy-1,jz)
!!!      satn = satliq(jx,jy+1,jz)
!!!      sats = satliq(jx,jy-1,jz)
      satn = 0.5*(satliq(jx,jy+1,jz) + satliqold(jx,jy+1,jz) )
      sats = 0.5*(satliq(jx,jy-1,jz) + satliqold(jx,jy-1,jz) )

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
        dumn = ro(jx,jy+1,jz) *(satn)**(SatPow)* (porn)**(PorPow) *dstar(jx,jy+1,jz)*anisotropyY
        dums = ro(jx,jy-1,jz) *(sats)**(SatPow)* (pors)**(PorPow) *dstar(jx,jy-1,jz)*anisotropyY
        dumpy = ro(jx,jy,jz)  *(satp)**(SatPow)* (porp)**(PorPow) *dstar(jx,jy,jz)  *anisotropyY
      ELSE
        dumn = ro(jx,jy+1,jz)*satn*porn*dstar(jx,jy+1,jz)*anisotropyY*tortuosity(jx,jy+1,jz)
        dums = ro(jx,jy-1,jz)*sats*pors*dstar(jx,jy-1,jz)*anisotropyY*tortuosity(jx,jy-1,jz)
        dumpy = ro(jx,jy,jz)*satp*porp*dstar(jx,jy,jz)*anisotropyY*tortuosity(jx,jy,jz)
      END IF
    END IF
    
    200     CONTINUE
    IF (nx == 1) GO TO 300
    
    IF (jx == 1) THEN
      !!!  EAST
      avgro = 0.5d0*( ro(jx+1,jy,jz) + ro(jx,jy,jz) )

      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dume,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dume,dumpx)
      ELSE
        dharm = GeometricMean(dume,dumpx)
      END IF

      AreaE = dyy(jy)*dzz(jx,jy,jz)
      
      dspe = avgro*dspx(jx,jy,jz) + dharm
      de = AreaE*dspe/dxe
      fe = AreaE*avgro*(qx(jx,jy,jz) + FluidBuryX(jx))
      ae = DMAX1(-fe,0.0D0) + de
      netflowX(1,jy,jz) = fe
      netDiffuseX(1,jy,jz) = de
      
      !!!  WEST
      
      avgro = ro(jx,jy,jz)
      dharm = dumpx
      AreaW = dyy(jy)*dzz(jx,jy,jz)
      dspw = avgro*dspx(jx-1,jy,jz) + dharm
      dw = AreaW*dspw/dxw
      fw = AreaW*avgro*(qx(jx-1,jy,jz) + FluidBuryX(jx-1))
      netflowX(0,jy,jz) = fw
      
      IF (jc(1) == 2 .or. JcByGrid(jx-1,jy,jz) == 2) THEN  
        aw = DMAX1(fw,0.0D0)       !  Pure advective boundary
        netDiffuseX(jx-1,jy,jz) = 0.0d0
      ELSE
        aw = DMAX1(fw,0.0D0) + dw
        netDiffuseX(jx-1,jy,jz) = dw
      END IF
      
      IF (jc(1) == 2 .or. JcByGrid(jx-1,jy,jz) == 2) THEN  
        apx = de +      DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)  !  Pure advective boundary
      ELSE
        apx = dw + de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
      END IF
      
    ELSE IF (jx == nx) THEN
         
   !!!  WEST
      
      avgro = 0.5d0*( ro(jx-1,jy,jz) + ro(jx,jy,jz) )

      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumw,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumw,dumpx)
      ELSE
        dharm = GeometricMean(dumw,dumpx)
      END IF

      AreaW = dyy(jy)*dzz(jx,jy,jz)
      
      dspw = avgro*dspx(nx-1,jy,jz) + dharm
      dw = AreaW*dspw/dxw
      fw = AreaW*avgro*(qx(nx-1,jy,jz) + FluidBuryX(jx-1))
      aw = DMAX1(fw,0.0D0) + dw
      netflowX(nx-1,jy,jz) = fw
      netDiffuseX(nx-1,jy,jz) = dw
      
  !!!  EAST
      
      avgro = ro(jx,jy,jz)
      dharm = dumpx
      AreaE = dyy(jy)*dzz(jx,jy,jz)
      dspe = avgro*dspx(jx+1,jy,jz) + dharm
      de = AreaE*dspe/dxe
      fe = AreaE*avgro*(qx(jx,jy,jz) + FluidBuryX(jx))
      netflowX(nx,jy,jz) = fe
      

      IF (jc(2) == 2 .or. JcByGrid(jx+1,jy,jz) == 2) THEN  

        ae = DMAX1(-fe,0.0D0)       !  Pure advective boundary
        netDiffuseX(jx,jy,jz) = 0.0d0
      ELSE
        ae = DMAX1(-fe,0.0D0) + de
        netDiffuseX(jx,jy,jz) = de
      END IF
      
      IF (jc(2) == 2 .or. JcByGrid(jx+1,jy,jz) == 2) THEN  
        apx = dw +      DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)       !  Pure advective boundary
      ELSE
        apx = dw + de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
      END IF
      
    ELSE     !!! Not jx /= 1 .or. jx /= nx
      
  !!!  EAST
      
      avgro = 0.5d0*( ro(jx+1,jy,jz) + ro(jx,jy,jz) ) 
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dume,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dume,dumpx)
      ELSE
        dharm = GeometricMean(dume,dumpx)
      END IF
      
      AreaE = dyy(jy)*dzz(jx,jy,jz)
      dspe = avgro*dspx(jx,jy,jz) + dharm
      de = AreaE*dspe/dxe
      fe = AreaE*avgro*(qx(jx,jy,jz) + FluidBuryX(jx))
      ae = DMAX1(-fe,0.0D0) + de
      netflowX(jx,jy,jz) = fe
      netDiffuseX(jx,jy,jz) = de
      
  !!!  WEST
      
      avgro = 0.5d0*( ro(jx-1,jy,jz) + ro(jx,jy,jz) )  
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumw,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumw,dumpx)
      ELSE
        dharm = GeometricMean(dumw,dumpx)
      END IF
      
      AreaW = dyy(jy)*dzz(jx,jy,jz)
      dspw = avgro*dspx(jx-1,jy,jz) + dharm
      dw = AreaW*dspw/dxw
      fw = AreaW*avgro*(qx(jx-1,jy,jz) + FluidBuryX(jx-1))
      aw = DMAX1(fw,0.0D0) + dw
      
      apx = dw + de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
      
    END IF
    
!!!   ********* Now JY coordinate direction *****************************
    
    300     CONTINUE
    IF (ny == 1) GO TO 400
    
    IF (jy == 1) THEN
      
      avgro = 0.5d0*( ro(jx,jy+1,jz) + ro(jx,jy,jz) )
      
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumn,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumn,dumpy)
      ELSE
        dharm = GeometricMean(dumn,dumpy)
      END IF
      
      AreaN = dxx(jx)*dzz(jx,jy,jz)
      dspn = avgro*dspy(jx,jy,jz) + dharm
      dn = AreaN*dspn/dyn
      fn = AreaN*avgro*(qy(jx,jy,jz) + FluidBuryY(jy))
      an = DMAX1(-fn,0.0D0) + dn
      
      dharm = dumpy
      avgro = ro(jx,jy,jz)
      AreaS = dxx(jx)*dzz(jx,jy,jz)
      dsps = avgro*dspy(jx,jy,jz) + dharm
      ds = AreaS*dsps/dys
      fs = AreaS*avgro*(qy(jx,jy-1,jz) + FluidBuryY(jy-1))
      
      IF (jc(3) == 2 .or. JcByGrid(jx,jy-1,jz) == 2) THEN  
        as = DMAX1(fs,0.0D0)      ! Pure advective boundary
      ELSE
        as = DMAX1(fs,0.0D0) + ds
      END IF
      
      IF (jc(3) == 2 .or. JcByGrid(jx,jy-1,jz) == 2) THEN  
        apy = dn + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0) ! Pure advective boundary
      ELSE
        apy = ds + dn + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0)
      END IF
      
    ELSE IF (jy == ny) THEN
      
      avgro = ro(jx,jy,jz)
      dharm = dumpy
      AreaN = dxx(jx)*dzz(jx,jy,jz)
      dspn  = avgro*dspy(jx,jy-1,jz) + dharm
      dn    = AreaN*dspn/dyn
      fn    = AreaN*avgro*(qy(jx,jy,jz) + FluidBuryY(jy))
      
     IF (jc(4) == 2 .or. JcByGrid(jx,jy+1,jz) == 2) THEN  
        an = DMAX1(-fn,0.0D0)     ! Pure advective boundary
      ELSE
        an = DMAX1(-fn,0.0D0) + dn
      END IF
      
      avgro = 0.5d0*( ro(jx,jy-1,jz) + ro(jx,jy,jz) )

      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dums,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dums,dumpy)
      ELSE
        dharm = GeometricMean(dums,dumpy)
      END IF
      
      AreaS = dxx(jx)*dzz(jx,jy,jz)
      dsps  = avgro*dspy(jx,jy-1,jz) + dharm
      ds    = AreaS*dsps/dys
      fs    = AreaS*avgro*(qy(jx,jy-1,jz) + FluidBuryY(jy-1))
      as    = DMAX1(fs,0.0D0) + ds
      
    IF (jc(4) == 2 .or. JcByGrid(jx,jy+1,jz) == 2) THEN  
        apy = ds + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0) ! Pure advective boundary
      ELSE
        apy = ds + dn + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0)
      END IF
      
    ELSE     !!! JY /= 1 or NY
      
      avgro = 0.5d0*( ro(jx,jy+1,jz) + ro(jx,jy,jz) )
      
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumn,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumn,dumpy)
      ELSE
        dharm = GeometricMean(dumn,dumpy)
      END IF
      
      AreaN = dxx(jx)*dzz(jx,jy,jz)
      dspn = avgro*dspy(jx,jy,jz) + dharm
      dn = AreaN*dspn/dyn
      fn = AreaN*avgro*(qy(jx,jy,jz) + FluidBuryY(jy))
      an = DMAX1(-fn,0.0D0) + dn
      
      avgro = 0.5d0*( ro(jx,jy-1,jz) + ro(jx,jy,jz) )
      
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dums,dumpy)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dums,dumpy)
      ELSE
        dharm = GeometricMean(dums,dumpy)
      END IF
      
      AreaS = dxx(jx)*dzz(jx,jy,jz)
      dsps = avgro*dspy(jx,jy-1,jz) + dharm
      ds = AreaS*dsps/dys
      fs = AreaS*avgro*(qy(jx,jy-1,jz) + FluidBuryY(jy-1))
      as = DMAX1(fs,0.0D0) + ds
      
      apy = ds + dn + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0)
      
    END IF
    
    400     CONTINUE

    IF (nx == 1) THEN
      a(jx,jy,jz) = 0.0d0
      b(jx,jy,jz) = 0.0d0
      c(jx,jy,jz) = 0.0d0
    ELSE
      
      a(jx,jy,jz) = -aw
      c(jx,jy,jz) = -ae
      b(jx,jy,jz) = apx
      
    END IF
    
    IF (ny == 1) THEN
      d(jx,jy,jz) = 0.0d0
      f(jx,jy,jz) = 0.0d0
      e(jx,jy,jz) = 0.0d0
    ELSE
      
      d(jx,jy,jz) = -an
      f(jx,jy,jz) = -as
      e(jx,jy,jz) = apy
      
    END IF

    dxy(jx,jy,jz) = dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
    
  END DO
END DO

RETURN
END SUBROUTINE coeff
