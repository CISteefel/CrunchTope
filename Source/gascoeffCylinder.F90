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

SUBROUTINE gascoeffCylinder(nx,ny,nz)
USE crunchtype
USE medium
USE transport
USE concentration
USE temperature
USE flow
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
REAL(DP)                                      :: OneOver
REAL(DP)                                      :: pi
REAL(DP)                                      :: AreaE
REAL(DP)                                      :: AreaW
REAL(DP)                                      :: AreaS
REAL(DP)                                      :: AreaN
REAL(DP)                                      :: zero
REAL(DP)                                      :: gasd

REAL(DP)                                      :: PorPow
REAL(DP)                                      :: SatPow


INTEGER(I4B)                                  :: jx
INTEGER(I4B)                                  :: jy
INTEGER(I4B)                                  :: jz
INTEGER(I4B)                                  :: j

jz = 1

SatPow = 7.0/3.0
PorPow = 1.0/3.0

zero = 0.0
pi = DACOS(-1.0d0)

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
      sate = 1.0-satliq(jx+1,jy,jz)
      porw = por(jx,jy,jz)
      satw = 1.0-satliq(jx,jy,jz)
      gasd = (pore)**(PorPow)*(sate)**(SatPow)*dgas
      dume = pore*sate*gasd
      gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas
      dumpx = porp*satp*gasd
      dumw = dumpx
    ELSE IF (jx == nx) THEN
      dxw = 0.5*(dxx(jx)+dxx(jx-1))
      dxe = 0.5*dxx(nx)
      pore = por(jx,jy,jz)
      sate = 1.0-satliq(jx,jy,jz)
      porw = por(jx-1,jy,jz)
      satw = 1.0-satliq(jx-1,jy,jz)
      gasd = (porw)**(PorPow)*(satw)**(SatPow)*dgas
      dumw = porw*satw*gasd
      gasd = porp**(PorPow)*(satp)**(SatPow)*dgas
      dumpx = porp*satp*gasd
      dume = dumpx
    ELSE
      dxe = 0.5*(dxx(jx)+dxx(jx+1))
      dxw = 0.5*(dxx(jx)+dxx(jx-1))
      pore = por(jx+1,jy,jz)
      sate = 1.0-satliq(jx+1,jy,jz)
      porw = por(jx-1,jy,jz)
      satw = 1.0-satliq(jx-1,jy,jz)
      gasd = (pore)**(PorPow)*(sate)**(SatPow)*dgas
      dume = pore*sate*gasd
      gasd = (porw)**(PorPow)*(satw)**(SatPow)*dgas
      dumw = porw*satw*gasd
      gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas
      dumpx = porp*satp*gasd
    END IF
    
    100     CONTINUE
    IF (ny == 1) GO TO 200
    
    IF (jy == 1) THEN
      dyn = 0.5*(dyy(jy)+dyy(jy+1))
      dys = 0.5*dyy(1)
      porn = por(jx,jy+1,jz)
      satn = 1.0-satliq(jx,jy+1,jz)
      pors = por(jx,jy,jz)
      sats = 1.0-satliq(jx,jy,jz)
      gasd = (porn)**(PorPow)*(satn)**(SatPow)*dgas
      dumn = porn*satn*gasd
      gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas
      dumpy = porp*satp*gasd
      dums = dumpy
    ELSE IF (jy == ny) THEN
      dys = 0.5*(dyy(jy)+dyy(jy-1))
      dyn = 0.5*dyy(ny)
      porn = por(jx,jy,jz)
      satn = 1.0-satliq(jx,jy,jz)
      pors = por(jx,jy-1,jz)
      sats = 1.0-satliq(jx,jy-1,jz)
      gasd = (pors)**(PorPow)*(sats)**(SatPow)*dgas
      dums = pors*sats*gasd
      gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas
      dumpy = porp*satp*gasd
      dumn = dumpy
    ELSE
      dyn = 0.5*(dyy(jy)+dyy(jy+1))
      dys = 0.5*(dyy(jy)+dyy(jy-1))
      porn = por(jx,jy+1,jz)
      satn = 1.0-satliq(jx,jy+1,jz)
      pors = por(jx,jy-1,jz)
      sats = 1.0-satliq(jx,jy-1,jz)
      gasd = (pors)**(PorPow)*(sats)**(SatPow)*dgas
      dums = pors*sats*gasd
      gasd = (porn)**(PorPow)*(satn)**(SatPow)*dgas
      dumn = porn*satn*gasd
      gasd = (porp)**(PorPow)*(satp)**(SatPow)*dgas
      dumpy = porp*satp*gasd
    END IF
    
    200     CONTINUE
    IF (nx == 1) GO TO 300
    
    IF (jx == 1) THEN

      AreaE = dyy(jy)*2.0*pi*(x(jx) + dxx(jx)/2.0)
      AreaW = dyy(jy)*2.0*pi*(x(jx) - dxx(jx)/2.0)
      porsum = pore*sate+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pore*sate*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      dumsum = dume+dumpx
      IF (dumsum /= 0) THEN
        dharm = 2.0*dume*dumpx/(dume + dumpx)
      ELSE
        dharm = 0.0
      END IF

      dspe = dharm
      de = AreaE*dspe/dxe
      fe = AreaE*qxgas(jx,jy,jz)
      ae = DMAX1(-fe,0.0D0) + de
      
      porsum = porw*satw+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porw*satw*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      dspw = dumpx
      dw = AreaW*dspw/dxw
      fw = AreaW*qxgas(jx-1,jy,jz)
      IF (jc(1) == 2) THEN
        aw = DMAX1(fw,0.0D0)       !  Pure advective boundary
      ELSE
        aw = DMAX1(fw,0.0D0) + dw
      END IF
      
      IF (jc(1) == 2) THEN
        apx = de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)  !  Pure advective boundary
      ELSE
        apx = dw + de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
      END IF
      
    ELSE IF (jx == nx) THEN

      AreaE = dyy(jy)*2.0*pi*(x(jx) + dxx(jx)/2.0)
      AreaW = dyy(jy)*2.0*pi*(x(jx) - dxx(jx)/2.0)

      porsum = pore*sate+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pore*sate*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      dspe = dumpx
      de = AreaE*dspe/dxe
      fe = AreaE*qxgas(jx,jy,jz)
      IF (jc(2) == 2) THEN
        ae = DMAX1(-fe,0.0D0)      !  Pure advective boundary
      ELSE
        ae = DMAX1(-fe,0.0D0) + de
      END IF
      
      porsum = porw*satw+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porw*satw*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      dumsum = dumw+dumpx
      IF (dumsum /= 0) THEN
        dharm = 2.0*dumw*dumpx/(dumw + dumpx)
      ELSE
        dharm = 0.0
      END IF
      dspw = dharm
      dw = AreaW*dspw/dxw
      fw = AreaW*qxgas(jx-1,jy,jz)
      aw = DMAX1(fw,0.0D0) + dw
      
      IF (jc(2) == 2) THEN
        apx = dw + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)      !  Pure advective boundary
      ELSE
        apx = dw + de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
      END IF
      
    ELSE

      AreaE = dyy(jy)*2.0*pi*(x(jx) + dxx(jx)/2.0)
      AreaW = dyy(jy)*2.0*pi*(x(jx) - dxx(jx)/2.0)

      porsum = pore*sate+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pore*sate*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      dumsum = dume+dumpx
      IF (dumsum /= 0) THEN
        dharm = 2.0*dume*dumpx/(dumsum)
      ELSE
        dharm = 0.0
      END IF

      dspe = dharm
      de = AreaE*dspe/dxe
      fe = AreaE*qxgas(jx,jy,jz)
      ae = DMAX1(-fe,0.0D0) + de
      
      porsum = porw*satw+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porw*satw*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      dumsum = dumw+dumpx
      IF (dumsum /= 0) THEN
        dharm = 2.0*dumw*dumpx/(dumsum)
      ELSE
        dharm = 0.0
      END IF

      dspw = dharm
      dw = AreaW*dspw/dxw
      fw = AreaW*qxgas(jx-1,jy,jz)
      aw = DMAX1(fw,0.0D0) + dw
      
      apx = dw + de + DMAX1(-fw,0.0D0) + DMAX1(fe,0.0D0)
      
    END IF
    
    300     CONTINUE
    IF (ny == 1) GO TO 400
    
    IF (jy == 1) THEN
      
      AreaN = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )
      AreaS = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )

      porsum = porn*satn+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porn*satn*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      dumsum = dumn+dumpy
      IF (dumsum /= 0) THEN
        dharm = 2.0*dumn*dumpy/(dumn + dumpy)
      ELSE
        dharm = 0.0
      END IF

      dspn = dharm
      dn = AreaN*dspn/dyn
      fn = AreaN*qygas(jx,jy,jz)
      an = DMAX1(-fn,0.0D0) + dn
      porsum = pors*sats+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pors*sats*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF

      dsps = dumpy
      ds = AreaS*dsps/dys
      fs = AreaS*qygas(jx,jy-1,jz)
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
      
      AreaN = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )
      AreaS = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )

      porsum = porn*satn+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porn*satn*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      dspn = dumpy
      dn = AreaN*dspn/dyn
      fn = AreaN*qygas(jx,jy,jz)
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
      dumsum = dums+dumpy
      IF (dumsum /= 0) THEN
        dharm = 2.0*dums*dumpy/(dums + dumpy)
      ELSE
        dharm = 0.0
      END IF
      dsps = dharm
      ds = AreaS*dsps/dys
      fs = AreaS*qygas(jx,jy-1,jz)
      as = DMAX1(fs,0.0D0) + ds
      
      IF (jc(4) == 2) THEN
        apy = ds + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0) ! Pure advective boundary
      ELSE
        apy = ds + dn + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0)
      END IF
      
    ELSE
      
      AreaN = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )
      AreaS = pi*( (x(jx)+0.5*dxx(jx) )**2 - ( x(jx)-0.5*dxx(jx) )**2  )

      porsum = porn*satn+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*porn*satn*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      dumsum = dumn+dumpy
      IF (dumsum /= 0) THEN
        dharm = 2.0*dumn*dumpy/(dumn + dumpy)
      ELSE
        dharm = 0.0
      END IF
      dspn = dharm
      dn = AreaN*dspn/dyn
      fn = AreaN*qygas(jx,jy,jz)
      an = DMAX1(-fn,0.0D0) + dn
      
      porsum = pors*sats+porp*satp
      IF (porsum /= 0) THEN
        porharm= 2.0*pors*sats*porp*satp/porsum
      ELSE
        porharm = 0.0
      END IF
      dumsum = dums+dumpy
      IF (dumsum /= 0) THEN
        dharm = 2.0*dums*dumpy/(dums + dumpy)
      ELSE
        dharm = 0.0
      END IF
      dsps = dharm
      ds = AreaS*dsps/dys
      fs = AreaS*qygas(jx,jy-1,jz)
      as = DMAX1(fs,0.0D0) + ds
      
      apy = ds + dn + DMAX1(-fs,0.0D0) + DMAX1(fn,0.0D0)
      
    END IF
    
    
    400     CONTINUE

    IF (nx == 1) THEN
      ag(jx,jy,jz) = 0.0
      bg(jx,jy,jz) = 0.0
      cg(jx,jy,jz) = 0.0
    ELSE
      
      ag(jx,jy,jz) = -aw
      cg(jx,jy,jz) = -ae
      bg(jx,jy,jz) = apx
      
    END IF
    
    IF (ny == 1) THEN
      dg(jx,jy,jz) = 0.0
      fg(jx,jy,jz) = 0.0
      eg(jx,jy,jz) = 0.0
    ELSE
      
      dg(jx,jy,jz) = -an
      fg(jx,jy,jz) = -as
      eg(jx,jy,jz) = apy
      
    END IF

    dxy(jx,jy,jz) = pi*dyy(jy)*( (x(jx)+0.5*dxx(jx))**2 - (x(jx)-0.5*dxx(jx))**2  )
    
  END DO
END DO


RETURN
END SUBROUTINE gascoeffCylinder
