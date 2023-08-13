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

SUBROUTINE AssembleLocal(ncomp,nspec,nkin,nrct,ngas,ikin,  &
    nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,delt,jx,jy,jz,nx,ny,nz,ne,time)
USE crunchtype
USE params
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE RunTime
!USE watsolv_mod
!USE watsolv_main
USE temperature

IMPLICIT NONE
!fp! auto_par_loops=0;
integer nx,ny,nz,ne

INTEGER(I4B), INTENT(IN)                      :: ncomp
INTEGER(I4B), INTENT(IN)                      :: nspec
INTEGER(I4B), INTENT(IN)                      :: nkin
INTEGER(I4B), INTENT(IN)                      :: nrct
INTEGER(I4B), INTENT(IN)                      :: ngas
INTEGER(I4B), INTENT(IN)                      :: ikin
INTEGER(I4B), INTENT(IN)                      :: nexchange
INTEGER(I4B), INTENT(IN)                      :: nexch_sec
INTEGER(I4B), INTENT(IN)                      :: nsurf
INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
INTEGER(I4B), INTENT(IN)                      :: npot
INTEGER(I4B), INTENT(IN)                      :: ndecay
INTEGER(I4B), INTENT(IN)                      :: jx
INTEGER(I4B), INTENT(IN)                      :: jy
INTEGER(I4B), INTENT(IN)                      :: jz


REAL(DP), INTENT(IN)                          :: delt
REAL(DP), INTENT(IN)                          :: time

!  Internal arrays

!  Internal real variables

REAL(DP)                                      :: r
REAL(DP)                                      :: satgas
REAL(DP)                                      :: df
REAL(DP)                                      :: source
REAL(DP)                                      :: sumrct
REAL(DP)                                      :: sumkin
REAL(DP)                                      :: rxnmin
REAL(DP)                                      :: rxnaq
REAL(DP)                                      :: aq_accum
REAL(DP)                                      :: source_jac
REAL(DP)                                      :: ex_accum
REAL(DP)                                      :: surf_accum
REAL(DP)                                      :: sum
REAL(DP)                                      :: sqrt_sion
REAL(DP)                                      :: portemp
REAL(DP)                                      :: rotemp
REAL(DP)                                      :: xgtemp
REAL(DP)                                      :: satl
REAL(DP)                                      :: faraday
REAL(DP)                                      :: delta_z
REAL(DP)                                      :: correct
REAL(DP)                                      :: Retardation

!  Internal integers

INTEGER(I4B)                                  :: neqn
INTEGER(I4B)                                  :: i
INTEGER(I4B)                                  :: jdum
INTEGER(I4B)                                  :: i2
INTEGER(I4B)                                  :: ix
INTEGER(I4B)                                  :: ix2
INTEGER(I4B)                                  :: is
INTEGER(I4B)                                  :: is2
INTEGER(I4B)                                  :: k
INTEGER(I4B)                                  :: ir
INTEGER(I4B)                                  :: np
INTEGER(I4B)                                  :: ind
INTEGER(I4B)                                  :: ind2
INTEGER(I4B)                                  :: ik
INTEGER(I4B)                                  :: npt2
INTEGER(I4B)                                  :: npt
INTEGER(I4B)                                  :: jpoint
INTEGER(I4B)                                  :: nrow
INTEGER(I4B)                                  :: ncol
INTEGER(I4B)                                  :: ns

REAL(DP)                                                  :: CellVolume
REAL(DP)                                                  :: MultiplyCell
REAL(DP)                                                  :: pi

pi = DACOS(-1.0d0)

aaa = 0.0 
source = 0.0

!! CellVolume = dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
!! MultiplyCell = CellVolume

Retardation = 0.001d0*SolidDensitygrid(jx,jy,jz)*(1.0-por(jx,jy,jz))/por(jx,jy,jz)

r = 1.0/delt      
neqn = ncomp + nsurf + nexchange + npot
IF (cylindrical) THEN
  CellVolume = dyy(jy)*pi*( (x(jx)+dxx(jx)/2.0d0 )**2.0d0 - ( x(jx)-dxx(jx)/2.0d0 )**2.0d0  )
  df = 1.0d0
  MultiplyCell = CellVolume
  IF (CylindricalDivideVolume) THEN
      df = 1.0/CellVolume
      MultiplyCell = 1.0
  END IF
!!  df = 1.0/(3.1416*dxx(jx)*dxx(jx)*dzz(jx,jy,jz))
ELSE IF (spherical) THEN
  df = 1.0
ELSE
#ifndef ALQUIMIA
   df = 1.0/(dxx(jx)*dyy(jy)*dzz(jx,jy,jz))
#endif
END IF
satl = satliq(jx,jy,jz)
satgas = 1.0 - satl
portemp = por(jx,jy,jz)
rotemp = ro(jx,jy,jz)
xgtemp = xgram(jx,jy,jz)
jpoint = jinit(jx,jy,jz)
faraday = 96485.0

IF (ActiveCell(jx,jy,jz) == 0) GO TO 800

!  Surface charge calculation

CALL SurfaceCharge(ncomp,nspec,nsurf,nsurf_sec,npot,jx,jy,jz,time)

sum = 0.0D0
DO ik = 1,ncomp+nspec
  sum = sum + sp10(ik,jx,jy,jz)*chg(ik)*chg(ik)
END DO
sqrt_sion = SQRT(0.50D0*sum)

DO npt = 1,npot
  ind = npt+ncomp+nexchange+nsurf
  fxx(ind) = 0.1174*sqrt_sion*SINH( LogPotential(npt,jx,jy,jz) ) - surfcharge(ksurf(ispot(npt)))
END DO

DO i = 1,ncomp
  ind = i

  IF (nradmax > 0) THEN
    sumrct = 0.0
    DO k = 1,nkin
      DO np = 1,nreactmin(k)
        IF (rmin(np,k) >= 0.0) then
          sumrct = sumrct + decay_correct(i,k)*mumin(np,k,i)*rmin(np,k)
        ELSE
          sumrct = sumrct + decay_correct(i,k)*mumin_decay(np,k,i,jx,jy,jz)*rmin(np,k)
        END IF
      END DO
    END DO
  ELSE
    sumrct = 0.0
    DO k = 1,nkin
      DO np = 1,nreactmin(k)   
        sumrct = sumrct + decay_correct(i,k)*mumin(np,k,i)*rmin(np,k)
      END DO
    END DO
  ENDIF
      
  sumkin = 0.0
  DO ir = 1,ikin
    sumkin = sumkin - mukin(ir,i)*raq_tot(ir,jx,jy,jz)
  END DO
      
!  Update the residual, adding reaction terms and exchange terms
      
  fxx(ind) = fxx(ind) + sumrct  &
    + satl*xgtemp*portemp*rotemp*sumkin
      
  sumrd = 0.0
  sumjackin = 0.0
      
  IF (nradmax > 0) THEN
    DO k = 1,nkin
      DO np = 1,nreactmin(k)
        IF (mumin(np,k,i) /= 0.0) THEN
          IF (rmin(np,k) >= 0.0) THEN
            DO i2 = 1,ncomp+nexchange+nsurf
              sumrd(i2) = sumrd(i2) + decay_correct(i,k)*mumin(np,k,i)*jac_rmin(i2,np,k)
            END DO
          ELSE
            DO i2 = 1,ncomp+nexchange+nsurf
              sumrd(i2) = sumrd(i2) + decay_correct(i,k)*mumin_decay(np,k,i,jx,jy,jz)*jac_rmin(i2,np,k)
            END DO
          END IF 
        END IF
      END DO
    END DO
  ELSE
    DO k = 1,nkin
      DO np = 1,nreactmin(k)
        IF (mumin(np,k,i) /= 0.0) THEN
          DO i2 = 1,ncomp+nexchange+nsurf
            sumrd(i2) = sumrd(i2) + decay_correct(i,k)*mumin(np,k,i)*jac_rmin(i2,np,k)
          END DO
        END IF
      END DO
    END DO
  END IF

  DO ir = 1,ikin
    IF (mukin(ir,i) /= 0.0) THEN
      DO i2 = 1,ncomp
        sumjackin(i2) = sumjackin(i2) - mukin(ir,i)*rdkin(ir,i2)
      END DO
    END IF
  END DO
      
  DO i2 = 1,ncomp             
    rxnmin = sumrd(i2)
    rxnaq = satl*xgtemp*portemp*rotemp*sumjackin(i2)

    IF (i /= ikh2o) THEN
      aq_accum = H2Oreacted(jx,jy,jz)*satl*xgtemp*r*portemp*rotemp*   &
          fjac_loc(i2,i) *(1.0 + Retardation*distrib(i) ) 
    ELSE
      aq_accum = satl*xgtemp*r*portemp*rotemp*   &
          fjac_loc(i2,i) *(1.0 + Retardation*distrib(i) ) 
    END IF

    source_jac = source*fjac_loc(i2,i)     
    ex_accum = r*fch_local(i,i2)  
    aaa(i,i2) = rxnmin + rxnaq + aq_accum - source_jac + ex_accum 
  END DO   ! end of I2 loop
      
  DO ix2 = 1,nexchange
    ind2 = ix2 + ncomp
    rxnmin = sumrd(ix2+ncomp)
    ex_accum = r*fch_local(i,ix2+ncomp)     
    aaa(i,ind2) = ex_accum + rxnmin
  END DO
      
  DO is2 = 1,nsurf      
    ind2 = is2+ncomp+nexchange
    rxnmin = sumrd(is2+ncomp+nexchange)
    surf_accum = r*fch_local(i,is2+ncomp+nexchange)
    aaa(i,ind2) = surf_accum + rxnmin
  END DO

  DO npt2 = 1,npot
    ind2 = npt2+ncomp+nexchange+nsurf
!    rxnmin = sumrd(is2+ncomp+nexchange)
    aaa(i,ind2) = r*fjpotncomp_local(npt2,i)
  END DO
      
  IF (isaturate == 1) THEN
    DO i2 = 1,ncomp         
      aaa(i,i2) = aaa(i,i2) + satgas*portemp*r*fgas_local(i2,i) 
    END DO   ! end of I2 loop
  END IF

END DO     !   End of I loop
    
DO ix = 1,nexchange
  ind = ix+ncomp
  IF (iexc == 1 .OR. iexc == 3) THEN
!!    fxx(ind) = (por(jx,jy,jz)*satliq(jx,jy,jz)*ro(jx,jy,jz))*totex(ix) - exchangesites(ix,jx,jy,jz)
    fxx(ind) = (sumactivity(ix) - 1.0)
  ELSE
    fxx(ind) = sumactivity(ix) - 1.0
  END IF
  DO i2 = 1,ncomp+nexchange
    ind2 = i2
    aaa(ix+ncomp,ind2) = fexch(ix+ncomp,i2)
  END DO

END DO     !  End of IX loop
    
DO is = 1,nsurf
      
  DO i2 = 1,ncomp
    ind2 = i2
    surf_accum = r*fsurf_local(is,i2)  
    aaa(is+ncomp+nexchange,ind2) = surf_accum
  END DO
      
  DO is2 = 1,nsurf
    ind2 = is2+ncomp+nexchange
    surf_accum = r*fsurf_local(is,is2+ncomp)        
    aaa(is+ncomp+nexchange,ind2) = surf_accum
  END DO
      
  DO npt2 = 1,npot
    ind2 = npt2+ncomp+nexchange+nsurf
    aaa(is+ncomp+nexchange,ind2) =r*fjpotnsurf_local(npt2,is)
  END DO

END DO     !   End of IS loop

DO npt = 1,npot
  nrow = npt + ncomp + nexchange + nsurf
  is = ispot(npt)
  k = ksurf(is)
  correct = wtmin(k)*specificByGrid(k,jx,jy,jz)*volfx(k,jx,jy,jz)/volmol(k) 
 
  DO i2 = 1,ncomp
    ind2 = i2
    sum = 0.0
    DO ns = 1,nsurf_sec
      IF (ksurf(islink(ns)) == ksurf(is)) THEN
!!      IF (islink(ns) == is) THEN
        sum = sum - musurf(ns,i2)*zsurf(ns+nsurf)*spsurf10(ns+nsurf,jx,jy,jz)*faraday/correct
      END IF
    END DO
    aaa(nrow,ind2) = sum
  END DO

  DO is2 = 1,nsurf
    ind2 = is2+ncomp+nexchange
    sum = 0.0
    DO ns = 1,nsurf_sec
      IF (ksurf(islink(ns)) == ksurf(is)) THEN
!!        IF (islink(ns) == is) THEN
        sum = sum - musurf(ns,is2+ncomp)*zsurf(ns+nsurf)*spsurf10(ns+nsurf,jx,jy,jz)*faraday/correct
      END IF
    END DO
    aaa(nrow,ind2) = sum
    IF (ksurf(is2) == ksurf(is)) THEN
       aaa(nrow,ind2) =  aaa(nrow,ind2)  &
             - zsurf(is2)*spsurf10(is2,jx,jy,jz)*faraday/correct
        END IF
  END DO

END DO

DO npt = 1,npot
  nrow = npt + ncomp + nexchange + nsurf
  is = ispot(npt)
  k = ksurf(is)
  correct = wtmin(k)*specificByGrid(k,jx,jy,jz)*volfx(k,jx,jy,jz)/volmol(k) 

  DO npt2 = 1,npot
    ncol = npt2 + ncomp + nexchange + nsurf
    IF (ksurf(ispot(npt)) == ksurf(ispot(npt2))) THEN
      sum = 0.0
      DO ns = 1,nsurf_sec
        delta_z = zsurf(ns+nsurf) - zsurf(islink(ns))
        IF (islink(ns) == ispot(npt2)) THEN
          sum = sum - zsurf(ns+nsurf)*spsurf10(ns+nsurf,jx,jy,jz)*delta_z*2.0
        END IF
      END DO
    END IF

    IF (nrow == ncol) THEN
      aaa(nrow,ncol) = 0.1174*sqrt_sion*COSH(LogPotential(npt,jx,jy,jz)) -     & 
             sum*faraday/correct
    ELSE
      aaa(nrow,ncol) =  -sum*faraday/correct
    END IF

  END DO

END DO


800 IF (activecell(jx,jy,jz) == 0) THEN
  fxx = 0.0
  aaa = 0.0
  DO i = 1,ncomp
    fxx(i) = sp10(i,jx,jy,jz) - spcond10(i,jinit(jx,jy,jz))
  END DO
  DO ix = 1,nexchange
    fxx(ix+ncomp) = spex(ix,jx,jy,jz) - spcondex(ix,jinit(jx,jy,jz))*por(jx,jy,jz)*ro(jx,jy,jz)*satliq(jx,jy,jz) 
  END DO
  DO ix = 1,nexch_sec
    spex10(ix+nexchange,jx,jy,jz) = spcondex10(ix+nexchange,jinit(jx,jy,jz))*por(jx,jy,jz)*ro(jx,jy,jz)*satliq(jx,jy,jz)
  END DO
  DO is = 1,nsurf
    fxx(is+ncomp+nexchange) = spsurf10(is,jx,jy,jz) - spcondsurf10(is,jinit(jx,jy,jz))*por(jx,jy,jz)*ro(jx,jy,jz)*satliq(jx,jy,jz)
  END DO
  DO i = 1,ncomp
    aaa(I,i) = sp10(i,jx,jy,jz)
  END DO
  DO ix = 1,nexchange
    aaa(ix+ncomp,ix+ncomp) = 1.0
  END DO
END IF
     
RETURN
END SUBROUTINE AssembleLocal
