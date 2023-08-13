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

SUBROUTINE FxTopeLocal(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,nrct,  &
    nspec,ngas,neqn,dt,jx,jy,jz,nx,ny,nz,ne,time)
USE crunchtype
USE params
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE temperature
USE RunTime

IMPLICIT NONE
!fp! auto_par_loops=0;

!  External variables and arrays

integer nx,ny,nz,ne
INTEGER(I4B), INTENT(IN)                                  :: ncomp
INTEGER(I4B), INTENT(IN)                                  :: nexchange
INTEGER(I4B), INTENT(IN)                                  :: nexch_sec
INTEGER(I4B), INTENT(IN)                                  :: nsurf
INTEGER(I4B), INTENT(IN)                                  :: nsurf_sec
INTEGER(I4B), INTENT(IN)                                  :: nrct
INTEGER(I4B), INTENT(IN)                                  :: nspec
INTEGER(I4B), INTENT(IN)                                  :: ngas
INTEGER(I4B), INTENT(IN)                                  :: neqn
INTEGER(I4B), INTENT(IN)                                  :: jx
INTEGER(I4B), INTENT(IN)                                  :: jy
INTEGER(I4B), INTENT(IN)                                  :: jz

REAL(DP), INTENT(IN)                                      :: dt
REAL(DP), INTENT(IN)                                      :: time

!  Internal variables and arrays

REAL(DP), DIMENSION(ncomp)                                :: scorr
REAL(DP)                                                  :: r
REAL(DP)                                                  :: satgasnew
REAL(DP)                                                  :: satgasold
REAL(DP)                                                  :: df
REAL(DP)                                                  :: source
REAL(DP)                                                  :: aq_accum
REAL(DP)                                                  :: gas_accum
REAL(DP)                                                  :: ex_accum
REAL(DP)                                                  :: surf_accum
REAL(DP)                                                  :: satl
REAL(DP)                                                  :: satlold
REAL(DP)                                                  :: Retardation

INTEGER(I4B)                                              :: i
INTEGER(I4B)                                              :: is
INTEGER(I4B)                                              :: ix
INTEGER(I4B)                                              :: ind

REAL(DP)                                                  :: CellVolume
REAL(DP)                                                  :: MultiplyCell
REAL(DP)                                                  :: pi

!  This routine calculates (and assembles) the function residuals.

pi = DACOS(-1.0d0)

r = 1.0D0/dt

IF (cylindrical) THEN
  CellVolume = dyy(jy)*pi*( (x(jx)+dxx(jx)/2.0d0 )**2.0d0 - ( x(jx)-dxx(jx)/2.0d0 )**2.0d0  )
  df = 1.0d0
  MultiplyCell = CellVolume
  IF (CylindricalDivideVolume) THEN
      df = 1.0/CellVolume
      MultiplyCell = 1.0
  END IF

ELSE
#ifndef ALQUIMIA
   df = 1.0/(dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
#endif
END IF

!!   CellVolume = dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
!!   MultiplyCell = CellVolume

Retardation = 0.001d0*SolidDensitygrid(jx,jy,jz)*(1.0-por(jx,jy,jz))/por(jx,jy,jz)

fxmax = 0.0
satl = satliq(jx,jy,jz)
satgasnew = 1.0 - satl
satlOld = satl
!!satlold = satliqold(jx,jy,jz)
satgasold = 1.0 - satlold

DO i = 1,ncomp  

  ind = i

  source = 0.0
  
!!      Assumes Kd is dimensionless  

!!  aq_accum = xgram(jx,jy,jz)*r*por(jx,jy,jz)*ro(jx,jy,jz)*   &
!!           (satl*s(i,jx,jy,jz)-satlold*sn(i,jx,jy,jz)) *(1.0 + Retardation*distrib(i) ) 

  if (i /= ikh2o) THEN
    aq_accum = r*por(jx,jy,jz)*ro(jx,jy,jz)*   &
           ( H2Oreacted(jx,jy,jz) * xgram(jx,jy,jz) * satl * s(i,jx,jy,jz) * ( 1.0 + Retardation*distrib(i) ) - &
             xgramOld(jx,jy,jz) * satlold * sn(i,jx,jy,jz) ) - &
            r * skdold(i,jx,jy,jz)
  ELSE
    aq_accum = r*por(jx,jy,jz)*ro(jx,jy,jz)*   &
           ( xgram(jx,jy,jz) * satl * s(i,jx,jy,jz) * ( 1.0 + Retardation*distrib(i) ) - &
             xgramOld(jx,jy,jz) * satlold * sn(i,jx,jy,jz) ) - &
           r * skdold(i,jx,jy,jz) 
  END IF

  IF (isaturate == 1) THEN
    gas_accum = por(jx,jy,jz)*r*(satgasnew*sgas(i,jx,jy,jz)-satgasold*sgasn(i,jx,jy,jz))
  ELSE
    gas_accum = 0.0
  END IF

  ex_accum = r*(sNCexch_local(i) - sexold(i,jx,jy,jz))      &
            +    r*(sNCsurf_local(i) - ssurfold(i,jx,jy,jz))

  fxx(ind) = aq_accum + gas_accum + ex_accum - source 
  
  IF (ABS(fxx(ind)) > fxmax(i)) THEN
    fxmax(i) = ABS(fxx(ind))
  END IF
  
END DO

DO is = 1,nsurf
  ind = is+ncomp+nexchange
  surf_accum = r*(ssurf_local(is) - ssurfn(is,jx,jy,jz))
  fxx(ind) = surf_accum
END DO

RETURN
END SUBROUTINE FxTopeLocal
!***************************************************************************
