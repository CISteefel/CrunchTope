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

SUBROUTINE coeffSphericalNew_d(nx,ny,nz,ncomp,nspec)
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
INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nspec

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
REAL(DP)                                      :: porharm
REAL(DP)                                      :: porsum
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
REAL(DP)                                      :: quirk
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
REAL(DP)                                      :: OneOver
REAL(DP)                                      :: pi
REAL(DP)                                      :: AreaE
REAL(DP)                                      :: AreaW
REAL(DP)                                      :: AreaS
REAL(DP)                                      :: AreaN
REAL(DP)                                      :: tort

REAL(DP)                                      :: PorPow
REAL(DP)                                      :: SatPow

INTEGER(I4B)                                  :: jx
INTEGER(I4B)                                  :: jy
INTEGER(I4B)                                  :: jz
INTEGER(I4B)                                  :: j
INTEGER(I4B)                                  :: ik

jz = 1
jy = 1     !!  Assumes a 1D system with X as radial direction 

pi = DACOS(-1.0d0)

SatPow = 10.0d0/3.0d0
PorPow = 4.0d0/3.0d0

IF (idiffus == 0) THEN
  d_25 = dzero
ELSE
  d_25 = dcoeff
END IF

DO ik = 1,ncomp+nspec
  d_correct(ik) = d_sp(ik)/d_25
END DO

jz = 1
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

dstar(0,jy,jz) = dstar(1,jy,jz)
dstar(nx+1,jy,jz) = dstar(nx,jy,jz)

!!DO jy = 1,ny
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
        dumw = dumpx
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
        dume = dumpx
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
    
    AreaE = 4.0*pi*(x(jx) + dxx(jx)/2.0)**2
    AreaW = 4.0*pi*(x(jx) - dxx(jx)/2.0)**2

    IF (nx == 1) GO TO 300
    
    IF (jx == 1) THEN

      porsum = pore*sate+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pore*sate*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = 0.5*( ro(jx+1,jy,jz) + ro(jx,jy,jz) )
      dumsum = dume+dumpx
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dume,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dume,dumpx)
      ELSE
        dharm = GeometricMean(dume,dumpx)
      END IF
      dspe = dharm
      de = AreaE*dspe/dxe
      ae = de
      
      porsum = porw*satw+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porw*satw*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = ro(jx,jy,jz)
      dspw = dumpx
      dw = AreaW*dspw/dxw
      IF (jc(1) == 2) THEN
        aw = 0.0
      ELSE
        aw = dw
      END IF
      
      IF (jc(1) == 2) THEN
        apx = de               !  No diffusion flux across boundary
      ELSE
        apx = dw + de 
      END IF
      
    ELSE IF (jx == nx) THEN

      porsum = pore*sate+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pore*sate*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = ro(jx,jy,jz)
      dspe = dumpx
      de = AreaE*dspe/dxe
      IF (jc(2) == 2) THEN
        ae = 0.0
      ELSE
        ae = de
      END IF
      
      porsum = porw*satw+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porw*satw*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = 0.5*( ro(jx-1,jy,jz) + ro(jx,jy,jz) )
      dumsum = dumw+dumpx
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumw,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumw,dumpx)
      ELSE
        dharm = GeometricMean(dumw,dumpx)
      END IF
      dspw = dharm
      dw = AreaW*dspw/dxw
      aw = dw
      
      IF (jc(2) == 2) THEN
        apx = dw         !  No diffusion flux across boundary
      ELSE
        apx = dw + de 
      END IF
      
    ELSE

      porsum = pore*sate+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pore*sate*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = 0.5*( ro(jx+1,jy,jz) + ro(jx,jy,jz) )
      dumsum = dume+dumpx
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dume,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dume,dumpx)
      ELSE
        dharm = GeometricMean(dume,dumpx)
      END IF
      dspe = dharm
      de = AreaE*dspe/dxe
      ae = de
      
      porsum = porw*satw+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porw*satw*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = 0.5*( ro(jx-1,jy,jz) + ro(jx,jy,jz) )
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumw,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumw,dumpx)
      ELSE
        dharm = GeometricMean(dumw,dumpx)
      END IF
      dspw = dharm
      dw = AreaW*dspw/dxw
      aw = dw
      
      apx = dw + de 
      
    END IF
    
    300     CONTINUE

    IF (nx == 1) THEN
      a(jx,jy,jz) = 0.0
      b(jx,jy,jz) = 0.0
      c(jx,jy,jz) = 0.0

      a_d(jx,jy,jz) = 0.0
      b_d(jx,jy,jz) = 0.0
      c_d(jx,jy,jz) = 0.0

    ELSE
      
      a(jx,jy,jz) = 0.0               !!  Assumes for the moment no flow
      b(jx,jy,jz) = 0.0
      c(jx,jy,jz) = 0.0

      a_d(jx,jy,jz) = -aw
      b_d(jx,jy,jz) = apx
      c_d(jx,jy,jz) = -ae
      
    END IF
    
    d(jx,jy,jz) = 0.0
    f(jx,jy,jz) = 0.0
    e(jx,jy,jz) = 0.0

    d_d(jx,jy,jz) = 0.0
    f_d(jx,jy,jz) = 0.0
    e_d(jx,jy,jz) = 0.0

    dxy(jx,jy,jz) = (4.0/3.0)*pi*( (x(jx)+0.5*dxx(jx))**3 - (x(jx)-0.5*dxx(jx))**3  )
    
  END DO
!!END DO


RETURN
END SUBROUTINE coeffSphericalNew_d
