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
SUBROUTINE coeffNuft(nx,ny,nz)
USE crunchtype
USE medium
USE transport
USE concentration
USE temperature
USE params

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

REAL(DP)                                      :: PorPow
REAL(DP)                                      :: SatPow

INTEGER(I4B)                                  :: jx
INTEGER(I4B)                                  :: jy
INTEGER(I4B)                                  :: jz
INTEGER(I4B)                                  :: j

REAL(DP)                                      :: erodex
REAL(DP)                                      :: erodey

!!  Set erosion to zero here, add appropriate erosion/burial arrays later

erodex = 0.0
erodey = 0.0

jz = 1

SatPow = 10.0d0/3.0d0
PorPow = 4.0d0/3.0d0

IF (idiffus == 0) THEN
  d_25 = dzero
ELSE
  d_25 = dcoeff
END IF

jz = 1
DO jy = 1,ny
  DO jx = 1,nx
!fp! set_index({#ident# jz #});
    tk = 273.15 + t(jx,jy,jz)
    IF (idiffus == 0) THEN
      dstar(jx,jy,jz) = dzero*EXP((activation/rgas)*(tk25 - 1.0/tk))/formation
    ELSE
      dstar(jx,jy,jz) = dcoeff/formation
    END IF
  END DO
END DO

jz = 1
DO jy = 1,ny
  DO jx = 1,nx
!fp! set_index({#ident# jz #});

    porp = por(jx,jy,jz)
    satp = satliq(jx,jy,jz)
    
    IF (nx == 1) GO TO 100
    
    IF (jx == 1) THEN
      
      dxe = 0.5*(dxx(jx)+dxx(jx+1))
      dxw = 0.5*dxx(1)
      pore = por(jx+1,jy,jz)
      porw = por(jx,jy,jz)
      sate = satliq(jx+1,jy,jz)
      satw = satliq(jx,jy,jz)
      dume = ro(jx+1,jy,jz)*(sate)**(SatPow)*(pore)**(PorPow)*dstar(jx+1,jy,jz)
      dumpx = ro(jx,jy,jz)*(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)
      dumw = dumpx
    ELSE IF (jx == nx) THEN
      dxw = 0.5*(dxx(jx)+dxx(jx-1))
      dxe = 0.5*dxx(nx)
      pore = por(jx,jy,jz)
      porw = por(jx-1,jy,jz)
      sate = satliq(jx,jy,jz)
      satw = satliq(jx-1,jy,jz)
      dumw = ro(jx-1,jy,jz)*(satw)**(SatPow)*(porw)**(PorPow)*dstar(jx-1,jy,jz)
      dumpx = ro(jx,jy,jz)*(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)
      dume = dumpx
    ELSE
      dxe = 0.5*(dxx(jx)+dxx(jx+1))
      dxw = 0.5*(dxx(jx)+dxx(jx-1))
      pore = por(jx+1,jy,jz)
      porw = por(jx-1,jy,jz)
      sate = satliq(jx+1,jy,jz)
      satw = satliq(jx-1,jy,jz)
      dume = ro(jx+1,jy,jz)*(sate)**(SatPow)*(pore)**(PorPow)*dstar(jx+1,jy,jz)
      dumpx = ro(jx,jy,jz)*(satp)**(SatPow)*(porp)**(PorPow)*dstar(jx,jy,jz)
      dumw = ro(jx-1,jy,jz)*(satw)**(SatPow)*(porw)**(PorPow)*dstar(jx-1,jy,jz)
    END IF
    
    100     CONTINUE
    IF (ny == 1) GO TO 200
    
    IF (jy == 1) THEN
      dyn = 0.5*(dyy(jy)+dyy(jy+1))
      dys = 0.5*dyy(1)
      porn = por(jx,jy+1,jz)
      pors = por(jx,jy,jz)
      satn = satliq(jx,jy+1,jz)
      sats = satliq(jx,jy,jz)
      dumn = ro(jx,jy+1,jz)*(satn)**(SatPow)*(por(jx,jy+1,jz))**(PorPow)*dstar(jx,jy+1,jz)
      dumpy = ro(jx,jy,jz)*(satp)**(SatPow)*(por(jx,jy,jz))**(PorPow)*dstar(jx,jy,jz)
      dums = dumpy
    ELSE IF (jy == ny) THEN
      dys = 0.5*(dyy(jy)+dyy(jy-1))
      dyn = 0.5*dyy(ny)
      porn = por(jx,jy,jz)
      pors = por(jx,jy-1,jz)
      satn = satliq(jx,jy,jz)
      sats = satliq(jx,jy-1,jz)
      dums = ro(jx,jy-1,jz)*(sats)**(SatPow)*(por(jx,jy-1,jz))**(PorPow)*dstar(jx,jy-1,jz)
      dumpy = ro(jx,jy,jz)*(satp)**(SatPow)*(por(jx,jy,jz))**(PorPow)*dstar(jx,jy,jz)
      dumn = dumpy
    ELSE
      dyn = 0.5*(dyy(jy)+dyy(jy+1))
      dys = 0.5*(dyy(jy)+dyy(jy-1))
      porn = por(jx,jy+1,jz)
      pors = por(jx,jy-1,jz)
      satn = satliq(jx,jy+1,jz)
      sats = satliq(jx,jy-1,jz)
      dumn = ro(jx,jy+1,jz)*(satn)**(SatPow)*(por(jx,jy+1,jz))**(PorPow)*dstar(jx,jy+1,jz)
      dums = ro(jx,jy-1,jz)*(sats)**(SatPow)*(por(jx,jy-1,jz))**(PorPow)*dstar(jx,jy-1,jz)
      dumpy = ro(jx,jy,jz)*(satp)**(SatPow)*(por(jx,jy,jz))**(PorPow)*dstar(jx,jy,jz)
    END IF
    
    200     CONTINUE
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
      IF (dumsum /= 0) THEN
        dharm = 2.0*dume*dumpx/(dume + dumpx)
      ELSE
        dharm = 0.0
      END IF
      dspe = dspx(jx,jy,jz) + dharm
      de = dspe*dyy(jy)/dxe
      fe = dyy(jy)*(qx(jx,jy,jz)+avgro*porharm*erodex)
      ae = DMAX1(-fe,0.0D0) + de
      netflowx(1,jy,jz) = qx(jx,jy,jz) + avgro*porharm*erodex
      
      porsum = porw*satw+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porw*satw*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = ro(jx,jy,jz)
      dspw = dspx(jx,jy,jz) + dumpx
      dw = dspw*dyy(jy)/dxw
      fw = dyy(jy)*(qx(jx-1,jy,jz)+avgro*porharm*erodex)
      IF (jc(1) == 2) THEN
        aw = DMAX1(fw,0.0D0)       !  Pure advective boundary
      ELSE
        aw = DMAX1(fw,0.0D0) + dw
      END IF
      
      netflowx(0,jy,jz) = qx(jx-1,jy,jz) + avgro*porharm*erodex
      
      IF (jc(1) == 2) THEN
        apx = de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)  !  Pure advective boundary
      ELSE
        apx = dw + de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
      END IF
      
    ELSE IF (jx == nx) THEN
      
      porsum = pore*sate+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pore*sate*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = ro(jx,jy,jz)
      dspe = dspx(jx,jy,jz) + dumpx
      de = dspe*dyy(jy)/dxe
      fe = dyy(jy)*(qx(jx,jy,jz)+avgro*porharm*erodex)
      IF (jc(2) == 2) THEN
        ae = DMAX1(-fe,0.0D0)      !  Pure advective boundary
      ELSE
        ae = DMAX1(-fe,0.0D0) + de
      END IF
      netflowx(jx,jy,jz) = qx(jx,jy,jz) + avgro*porharm*erodex
      
      porsum = porw*satw+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porw*satw*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = 0.5*( ro(jx-1,jy,jz) + ro(jx,jy,jz) )
      dumsum = dumw+dumpx
      IF (dumsum /= 0) THEN
        dharm = 2.0*dumw*dumpx/(dumw + dumpx)
      ELSE
        dharm = 0.0
      END IF
      dspw = dspx(jx-1,jy,jz) + dharm
      dw = dspw*dyy(jy)/dxw
      fw = dyy(jy)*(qx(jx-1,jy,jz)+avgro*porharm*erodex)
      aw = DMAX1(fw,0.0D0) + dw
      netflowx(jx-1,jy,jz) = qx(jx-1,jy,jz) + avgro*porharm*erodex
      
      IF (jc(2) == 2) THEN
        apx = dw + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)      !  Pure advective boundary
      ELSE
        apx = dw + de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
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
      IF (dumsum /= 0) THEN
        dharm = 2.0*dume*dumpx/(dumsum)
      ELSE
        dharm = 0.0
      END IF
      dspe = dspx(jx,jy,jz) + dharm
      de = dspe*dyy(jy)/dxe
      fe = dyy(jy)*(qx(jx,jy,jz)+avgro*porharm*erodex)
      ae = DMAX1(-fe,0.0D0) + de
      netflowx(jx,jy,jz) = qx(jx,jy,jz) + porharm*erodex
      
      porsum = porw*satw+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porw*satw*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = 0.5*( ro(jx-1,jy,jz) + ro(jx,jy,jz) )
      dumsum = dumw+dumpx
      IF (dumsum /= 0) THEN
        dharm = 2.0*dumw*dumpx/(dumsum)
      ELSE
        dharm = 0.0
      END IF
      dspw = dspx(jx-1,jy,jz) + dharm
      dw = dspw*dyy(jy)/dxw
      fw = dyy(jy)*(qx(jx-1,jy,jz)+avgro*porharm*erodex)
      aw = DMAX1(fw,0.0D0) + dw
      netflowx(jx-1,jy,jz) = qx(jx-1,jy,jz) + porharm*erodex
      
      apx = dw + de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
      
    END IF
    
    300     CONTINUE
    IF (ny == 1) GO TO 400
    
    IF (jy == 1) THEN
      
      porsum = porn*satn+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porn*satn*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = 0.5*( ro(jx,jy+1,jz) + ro(jx,jy,jz) )
      dumsum = dumn+dumpy
      IF (dumsum /= 0) THEN
        dharm = 2.0*dumn*dumpy/(dumn + dumpy)
      ELSE
        dharm = 0.0
      END IF
      dspn = dspy(jx,jy,jz) + dharm
      dn = dspn*dxx(jx)/dyn
      fn = dxx(jx)*(qy(jx,jy,jz)+avgro*porharm*erodey)
      an = DMAX1(-fn,0.0D0) + dn
      
      porsum = pors*sats+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pors*sats*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = ro(jx,jy,jz)
      dsps = dspy(jx,jy,jz) + dumpy
      ds = dsps*dxx(jx)/dys
      fs = dxx(jx)*(qy(jx,jy-1,jz)+avgro*porharm*erodey)
      IF (jc(3) == 2) THEN
        as = DMAX1(fs,0.0D0)      ! Pure advective boundary
      ELSE
        as = DMAX1(fs,0.0D0) + ds
      END IF
      
      IF (jc(3) == 2) THEN
        apy = dn + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0) ! Pure advective boundary
      ELSE
        apy = ds + dn + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0)
      END IF
      
    ELSE IF (jy == ny) THEN
      
      porsum = porn*satn+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porn*satn*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = ro(jx,jy,jz)
      dspn = dspy(jx,jy-1,jz) + dumpy
      dn = dspn*dxx(jx)/dyn
      fn = dxx(jx)*(qy(jx,jy,jz)+avgro*porharm*erodey)
      IF (jc(4) == 2) THEN
        an = DMAX1(-fn,0.0D0)     ! Pure advective boundary
      ELSE
        an = DMAX1(-fn,0.0D0) + dn
      END IF
      
      porsum = pors*sats+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pors*sats*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = 0.5*( ro(jx,jy-1,jz) + ro(jx,jy,jz) )
      dumsum = dums+dumpy
      IF (dumsum /= 0) THEN
        dharm = 2.0*dums*dumpy/(dums + dumpy)
      ELSE
        dharm = 0.0
      END IF
      dsps = dspy(jx,jy-1,jz) + dharm
      ds = dsps*dxx(jx)/dys
      fs = dxx(jx)*(qy(jx,jy-1,jz)+avgro*porharm*erodey)
      as = DMAX1(fs,0.0D0) + ds
      
      IF (jc(4) == 2) THEN
        apy = ds + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0) ! Pure advective boundary
      ELSE
        apy = ds + dn + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0)
      END IF
      
    ELSE
      
      porsum = porn*satn+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porn*satn*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = 0.5*( ro(jx,jy+1,jz) + ro(jx,jy,jz) )
      dumsum = dumn+dumpy
      IF (dumsum /= 0) THEN
        dharm = 2.0*dumn*dumpy/(dumn + dumpy)
      ELSE
        dharm = 0.0
      END IF
      dspn = dspy(jx,jy,jz) + dharm
      dn = dspn*dxx(jx)/dyn
      fn = dxx(jx)*(qy(jx,jy,jz)+avgro*porharm*erodey)
      an = DMAX1(-fn,0.0D0) + dn
      
      porsum = pors*sats+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pors*sats*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      avgro = 0.5*( ro(jx,jy-1,jz) + ro(jx,jy,jz) )
      dumsum = dums+dumpy
      IF (dumsum /= 0) THEN
        dharm = 2.0*dums*dumpy/(dums + dumpy)
      ELSE
        dharm = 0.0
      END IF
      dsps = dspy(jx,jy-1,jz) + dharm
      ds = dsps*dxx(jx)/dys
      fs = dxx(jx)*(qy(jx,jy-1,jz)+avgro*porharm*erodey)
      as = DMAX1(fs,0.0D0) + ds
      
      apy = ds + dn + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0)
      
    END IF
    
    
    400     CONTINUE

    IF (nx == 1) THEN
      a(jx,jy,jz) = 0.0
      b(jx,jy,jz) = 0.0
      c(jx,jy,jz) = 0.0
    ELSE
      
      a(jx,jy,jz) = -aw
      c(jx,jy,jz) = -ae
      b(jx,jy,jz) = apx
      
    END IF
    
    IF (ny == 1) THEN
      d(jx,jy,jz) = 0.0
      f(jx,jy,jz) = 0.0
      e(jx,jy,jz) = 0.0
    ELSE
      
      d(jx,jy,jz) = -an
      f(jx,jy,jz) = -as
      e(jx,jy,jz) = apy
      
    END IF

    dxy(jx,jy,jz) = dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
    
  END DO
END DO


RETURN
END SUBROUTINE coeffNuft
