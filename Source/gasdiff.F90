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
 
SUBROUTINE gasdiff(nx,ny,nz)
  USE crunchtype
  USE params
  USE concentration
  USE medium
  USE transport
  USE flow
  USE temperature
  USE CrunchFunctions
  
  IMPLICIT NONE
  
  !  External arrays and variables
  
  INTEGER(I4B), INTENT(IN)                      :: nx
  INTEGER(I4B), INTENT(IN)                      :: ny
  INTEGER(I4B), INTENT(IN)                      :: nz
  
  ! Internal arrays and variables
  
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
  REAL(DP)                                      :: fe
  REAL(DP)                                      :: ae
  REAL(DP)                                      :: fw
  REAL(DP)                                      :: aw
  REAL(DP)                                      :: fn
  REAL(DP)                                      :: an
  REAL(DP)                                      :: fs
  REAL(DP)                                      :: as
  REAL(DP)                                      :: tk
  REAL(DP)                                      :: porn
  REAL(DP)                                      :: pors
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
  REAL(DP)                                      :: zero
  REAL(DP)                                      :: gasd
  REAL(DP)                                      :: tempe
  REAL(DP)                                      :: tempw
  REAL(DP)                                      :: temps
  REAL(DP)                                      :: tempn
  REAL(DP)                                      :: dharm
  REAL(DP)                                      :: AreaE
  REAL(DP)                                      :: AreaW
  REAL(DP)                                      :: AreaS
  REAL(DP)                                      :: AreaN
  
  REAL(DP)                                      :: PorPow
  REAL(DP)                                      :: SatPow
  
  INTEGER(I4B)                                  :: jx
  INTEGER(I4B)                                  :: jy
  INTEGER(I4B)                                  :: jz
  INTEGER(I4B)                                  :: j
  LOGICAL(LGT)                                  :: gad_diff_T


  gad_diff_T = .TRUE.

  zero = 0.0d0 
  jz = 1
  
  SatPow = 7.0d0/3.0d0
  PorPow = 1.0d0/3.0d0
  
  !!! Added by Steefel July 17
  
  DO jy = 1,ny
    DO jx = 1,nx
      IF (satliq(jx,jy,jz) > 1.0d0) THEN
        satliq(jx,jy,jz) = 1.0d0
      END IF
    END DO
  END DO
  
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
  
      j = (jy-1)*nx+jx
      tk = 273.15 + t(jx,jy,jz)
      
      porp = por(jx,jy,jz)
    
      satp = 1.0-satliq(jx,jy,jz)    
!!!      satp = 0.5*( 1.0-satliq(jx,jy,jz) + 1.0-satliqold(jx,jy,jz) )
      
      IF (nx == 1) GO TO 100
      
      IF (jx == 1) THEN
        
        dxe = 0.5*(dxx(jx)+dxx(jx+1))
        dxw = 0.5*dxx(1)
        pore = por(jx+1,jy,jz)
        sate = 1.0-satliq(jx+1,jy,jz)
!!!        sate = 0.5*( 1.0-satliq(jx+1,jy,jz) + 1.0-satliqold(jx+1,jy,jz) )
        porw = por(jx,jy,jz)
        satw = 1.0-satliq(jx,jy,jz)
!!!        satw = 0.5*( 1.0-satliq(jx,jy,jz) + 1.0-satliqold(jx,jy,jz) )
        if (east_river .and. gad_diff_T) then
          gasd = (pore)**(PorPow)*(sate)**(SatPow)*dgas*(((t(jx+1,jy,jz)+273.15)/273.15)**1.81)
          dume = pore*sate*gasd
          gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas*(((t(jx,jy,jz)+273.15)/273.15)**1.81)
          dumpx = porp*satp*gasd
        else
          gasd = (pore)**(PorPow)*(sate)**(SatPow)*dgas
          dume = pore*sate*gasd
          gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas
          dumpx = porp*satp*gasd
        endif
        dumw = dumpx
        
      ELSE IF (jx == nx .and. gad_diff_T) THEN
        
        dxw = 0.5*(dxx(jx)+dxx(jx-1))
        dxe = 0.5*dxx(nx)
        pore = por(jx,jy,jz)
        sate = 1.0-satliq(jx,jy,jz)
!!!        sate = 0.5*( 1.0-satliq(jx,jy,jz) + 1.0-satliqold(jx,jy,jz) )
        porw = por(jx-1,jy,jz)
        satw = 1.0-satliq(jx-1,jy,jz)
!!!        satw = 0.5*( 1.0-satliq(jx-1,jy,jz) + 1.0-satliqold(jx-1,jy,jz) )
        if (east_river .and. gad_diff_T) then
          gasd = (porw)**(PorPow)*(satw)**(SatPow)*dgas*(((t(jx-1,jy,jz)+273.15)/273.15)**1.81)
          dumw = porw*satw*gasd
          gasd = porp**(PorPow)*(satp)**(SatPow)*dgas*(((t(jx,jy,jz)+273.15)/273.15)**1.81)
          dumpx = porp*satp*gasd
        else
          gasd = (porw)**(PorPow)*(satw)**(SatPow)*dgas
          dumw = porw*satw*gasd
          gasd = porp**(PorPow)*(satp)**(SatPow)*dgas
          dumpx = porp*satp*gasd
        endif
        dume = dumpx
        
      ELSE
        
        dxe = 0.5*(dxx(jx)+dxx(jx+1))
        dxw = 0.5*(dxx(jx)+dxx(jx-1))
        pore = por(jx+1,jy,jz)
        sate = 1.0-satliq(jx+1,jy,jz)
!!!        sate = 0.5*( 1.0-satliq(jx+1,jy,jz) + 1.0-satliqold(jx+1,jy,jz) )
        porw = por(jx-1,jy,jz)
        satw = 1.0-satliq(jx-1,jy,jz)
!!!        satw = 0.5*( 1.0-satliq(jx-1,jy,jz) + 1.0-satliqold(jx-1,jy,jz) )
        if (east_river .and. gad_diff_T) then
          gasd = (pore)**(PorPow)*(sate)**(SatPow)*dgas*(((t(jx+1,jy,jz)+273.15)/273.15)**1.81)
          dume = pore*sate*gasd
          gasd = (porw)**(PorPow)*(satw)**(SatPow)*dgas*(((t(jx-1,jy,jz)+273.15)/273.15)**1.81)
          dumw = porw*satw*gasd
          gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas*(((t(jx,jy,jz)+273.15)/273.15)**1.81)
          dumpx = porp*satp*gasd
        else
        gasd = (pore)**(PorPow)*(sate)**(SatPow)*dgas
        dume = pore*sate*gasd
        gasd = (porw)**(PorPow)*(satw)**(SatPow)*dgas
        dumw = porw*satw*gasd
        gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas
        dumpx = porp*satp*gasd
        endif
      END IF
      
      100     CONTINUE
      IF (ny == 1) GO TO 200
      
      IF (jy == 1) THEN
        dyn = 0.5*(dyy(jy)+dyy(jy+1))
        dys = 0.5*dyy(1)
        porn = por(jx,jy+1,jz)
        !satn = 1.0-satliq(jx,jy+1,jz)
        satn = 0.5*( 1.0-satliq(jx,jy+1,jz) + 1.0-satliqold(jx,jy+1,jz) )
        pors = por(jx,jy,jz)
        !sats = 1.0-satliq(jx,jy,jz)
        sats = 0.5*( 1.0-satliq(jx,jy,jz) + 1.0-satliqold(jx,jy,jz) )
        if (east_river .and. gad_diff_T) then
        gasd = (porn)**(PorPow)*(satn)**(SatPow)*dgas*(((t(jx,jy+1,jz)+273.15)/273.15)**1.81)
        dumn = porn*satn*gasd
        gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas*(((t(jx,jy,jz)+273.15)/273.15)**1.81)
        dumpy = porp*satp*gasd
        else
        gasd = (porn)**(PorPow)*(satn)**(SatPow)*dgas
        dumn = porn*satn*gasd
        gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas
        dumpy = porp*satp*gasd
        endif
        dums = dumpy
      ELSE IF (jy == ny) THEN
        dys = 0.5*(dyy(jy)+dyy(jy-1))
        dyn = 0.5*dyy(ny)
        porn = por(jx,jy,jz)
        !satn = 1.0-satliq(jx,jy,jz)
        satn = 0.5*( 1.0-satliq(jx,jy,jz) + 1.0-satliqold(jx,jy,jz) )
        pors = por(jx,jy-1,jz)
        !sats = 1.0-satliq(jx,jy-1,jz)
        sats = 0.5*( 1.0-satliq(jx,jy-1,jz) + 1.0-satliqold(jx,jy-1,jz) )
        if (east_river .and. gad_diff_T) then
        gasd = (pors)**(PorPow)*(sats)**(SatPow)*dgas*(((t(jx,jy-1,jz)+273.15)/273.15)**1.81)
        dums = pors*sats*gasd
        gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas*(((t(jx,jy,jz)+273.15)/273.15)**1.81)
        dumpy = porp*satp*gasd
        else
        gasd = (pors)**(PorPow)*(sats)**(SatPow)*dgas
        dums = pors*sats*gasd
        gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas
        dumpy = porp*satp*gasd
        endif
        dumn = dumpy
      ELSE
        dyn = 0.5*(dyy(jy)+dyy(jy+1))
        dys = 0.5*(dyy(jy)+dyy(jy-1))
        porn = por(jx,jy+1,jz)
        !satn = 1.0-satliq(jx,jy+1,jz)
        satn = 0.5*( 1.0-satliq(jx,jy+1,jz) + 1.0-satliqold(jx,jy+1,jz) )
        pors = por(jx,jy-1,jz)
        !sats = 1.0-satliq(jx,jy-1,jz)
        sats = 0.5*( 1.0-satliq(jx,jy-1,jz) + 1.0-satliqold(jx,jy-1,jz) )
        if (east_river .and. gad_diff_T) then
        gasd = (pors)**(PorPow)*(sats)**(SatPow)*dgas*(((t(jx,jy-1,jz)+273.15)/273.15)**1.81)
        dums = pors*sats*gasd
        gasd = (porn)**(PorPow)*(satn)**(SatPow)*dgas*(((t(jx,jy+1,jz)+273.15)/273.15)**1.81)
        dumn = porn*satn*gasd
        gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas*(((t(jx,jy,jz)+273.15)/273.15)**1.81)
        dumpy = porp*satp*gasd
        else
        gasd = (pors)**(PorPow)*(sats)**(SatPow)*dgas
        dums = pors*sats*gasd
        gasd = (porn)**(PorPow)*(satn)**(SatPow)*dgas
        dumn = porn*satn*gasd
        gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas
        dumpy = porp*satp*gasd
        endif
      END IF
      
      200     CONTINUE
      IF (nx == 1) GO TO 300
      
!!! NOTE: Below assumes gas is always a Dirichlet boundary condition.  
!!!         This agrees with what is is in FxTopeGlobal
      
                !!! Steefel Checked
        
    !!! East
        
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dume,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dume,dumpx)
      ELSE
        dharm = GeometricMean(dume,dumpx)
      END IF

      if (jx==nx) then
        continue
      end if

      AreaE = dyy(jy)*dzz(jx,jy,jz) 
      dspe = dharm
      de = AreaE*dspe/dxe
      fe = AreaE*qxgas(jx,jy,jz)
      ae = DMAX1(-fe,0.0D0) + de
      
   !!!  WEST
      
      IF (MeanDiffusion == 1) THEN
        dharm = ArithmeticMean(dumw,dumpx)
      ELSE IF (MeanDiffusion == 2) THEN
        dharm = HarmonicMean(dumw,dumpx)
      ELSE
        dharm = GeometricMean(dumw,dumpx)
      END IF

      AreaW = dyy(jy)*dzz(jx,jy,jz) 
      dspw = dharm
      dw = AreaW*dspw/dxw
      fw = AreaW*qxgas(jx-1,jy,jz)
      aw = DMAX1(fw,0.0D0) + dw

      apx = dw + de + DMAX1(-fw,zero) + DMAX1(fe,zero)
      
  !!  **************************************************************
      
      300     CONTINUE
      IF (ny == 1) GO TO 400
      
      tempn= dumn + dumpy
      temps = dums + dumpy
      IF (tempn /= zero) THEN
        dspn = 2.0D0*dumn*dumpy/(tempn)
      ELSE
        dspn = zero
      END IF
      IF (temps /= zero) THEN
        dsps = 2.0D0*dums*dumpy/(temps)
      ELSE
        dsps = zero
      END IF
      dn = dspn*dxx(jx)/dyn
      ds = dsps*dxx(jx)/dys
  
      IF (jy == 1 .AND. jc(3) == 2) THEN
        ds = 0.00
      END IF
  
      IF (jy == Ny .AND. jc(4) == 2) THEN
        dn = 0.00
      END IF
      
  !!  The following forces Dirichlet conditions for gases at all times (unless above is uncommented)
  
      fn = dxx(jx)*qygas(jx,jy,jz)
      an = DMAX1(-fn,zero) + dn
      
      fs = dxx(jx)*qygas(jx,jy-1,jz)
      as = DMAX1(fs,zero) + ds
      
      apy = ds + dn + DMAX1(-fs,zero) + DMAX1(fn,zero)
      
      400     CONTINUE
      
      IF (nx == 1) THEN
        ag(jx,jy,jz) = zero
        cg(jx,jy,jz) = zero
        bg(jx,jy,jz) = zero
      ELSE
        ag(jx,jy,jz) = -aw
        cg(jx,jy,jz) = -ae
        bg(jx,jy,jz) = apx
      END IF
      
      IF (ny == 1) THEN
        dg(jx,jy,jz) = zero
        fg(jx,jy,jz) = zero
        eg(jx,jy,jz) = zero
      ELSE
        dg(jx,jy,jz) = -an
        fg(jx,jy,jz) = -as
        eg(jx,jy,jz) = apy
      END IF
      
    END DO
  END DO
  
  RETURN
  END SUBROUTINE gasdiff
  