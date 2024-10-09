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

    
SUBROUTINE os3d_newton(ncomp,nspec,nkin,nrct,ngas,ikin,       &
    nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,neqn,igamma,   & 
    delt,corrmax,jx,jy,jz,iterat,icvg,nx,ny,nz,time,AqueousToBulk)
USE crunchtype
USE runtime, ONLY: H2Opresent, Duan, Duan2006
USE params
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE temperature

IMPLICIT NONE


!  **********************  INTERFACE BLOCKS  **************************
INTERFACE
  SUBROUTINE ludcmp90(a,indx,d,n)
  USE crunchtype
  REAL(DP), DIMENSION(:,:), intent(in out)                   :: a
  INTEGER(I4B), DIMENSION(:), intent(out)                    :: indx
  INTEGER(I4B), INTENT(IN)                                   :: n
  REAL(DP), intent(out)                                      :: d
  END SUBROUTINE ludcmp90
END INTERFACE

INTERFACE
  SUBROUTINE lubksb90(a,indx,b,n)
  USE crunchtype
  IMPLICIT NONE
  REAL(DP),  DIMENSION(:,:), INTENT(IN)                          :: a
  REAL(DP),  DIMENSION(:), INTENT(IN OUT)                        :: b
  INTEGER(I4B),  DIMENSION(:),INTENT(IN)                         :: indx
  INTEGER(I4B), INTENT(IN)                                       :: n
  END SUBROUTINE lubksb90
END INTERFACE
!  ********************************************************************
integer nx,ny,nz

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
INTEGER(I4B), INTENT(IN)                      :: neqn
INTEGER(I4B), INTENT(IN)                      :: igamma
INTEGER(I4B), INTENT(IN)                      :: jx
INTEGER(I4B), INTENT(IN)                      :: jy
INTEGER(I4B), INTENT(IN)                      :: jz
INTEGER(I4B), INTENT(OUT)                     :: iterat
INTEGER(I4B), INTENT(OUT)                     :: icvg

REAL(DP), INTENT(IN)                          :: delt
REAL(DP), INTENT(IN)                          :: corrmax
REAL(DP), INTENT(IN)                          :: time
REAL(DP), INTENT(IN)                          :: AqueousToBulk

!  Internal arrays
#if defined(ALQUIMIA)
INTEGER(I4B), PARAMETER                       :: newton=50
#else
INTEGER(I4B), PARAMETER                       :: newton=50
#endif
INTEGER(I4B), PARAMETER                       :: ione=1
INTEGER(I4B)                                  :: ind
INTEGER(I4B)                                  :: ne
INTEGER(I4B)                                  :: i
INTEGER(I4B)                                  :: i2
INTEGER(I4B)                                  :: ix
INTEGER(I4B)                                  :: is
INTEGER(I4B)                                  :: npt
INTEGER(I4B)                                  :: info
INTEGER(I4B)                                  :: round

REAL(DP), PARAMETER                           :: atol=1.e-09
REAL(DP), PARAMETER                           :: rtol=1.e-06
REAL(DP)                                      :: errmax
REAL(DP)                                      :: tolmax
REAL(DP)                                      :: det
REAL(DP)                                      :: aascale
REAL(DP)                                      :: check

CHARACTER (LEN=1)                             :: trans

trans = 'N'
icvg = 0
iterat = 0

round = 1

DO ne = 1,newton
  iterat = iterat + 1

  CALL xmassNodeByNode(jx,jy,jz,ncomp,nspec)
  IF (H2Opresent) THEN
    CALL WaterReactedNodeByNode(jx,jy,jz)
  END IF

  IF (igamma == 2) THEN
!!!    IF (Duan .OR. Duan2006) THEN
!!!      CALL gamma_co2(ncomp,nspec,ngas,jx,jy,jz)
!!!    ELSE
!!!      CALL gamma(ncomp,nspec,jx,jy,jz)
!!!    END IF
    
    CALL gammaUpdated(ncomp,nspec,jx,jy,jz)
  END IF

  CALL SpeciesLocal(ncomp,nspec,jx,jy,jz)
  CALL SurfLocal(ncomp,nsurf,nsurf_sec,jx,jy,jz,AqueousToBulk)
  CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
  IF (isaturate == 1) THEN
    CALL gases(ncomp,ngas,jx,jy,jz)
  END IF

  CALL totconc(ncomp,nspec,jx,jy,jz)

  CALL totexchange_local(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)
  CALL totsurf_local(ncomp,nsurf,nsurf_sec,jx,jy,jz)
  CALL tot_ex(ncomp,nexchange,nexch_sec,jx,jy,jz)
  IF (isaturate == 1) THEN
    CALL totgas(ncomp,nspec,ngas,jx,jy,jz)
  END IF

  CALL jac_local(ncomp,nspec,jx,jy,jz)
  IF (isaturate == 1) THEN
    CALL jacgas_local(ncomp,ngas,jx,jy,jz)
  END IF
  CALL jac_exchange_local(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,neqn,jx,jy,jz)
  CALL jacsurf_local(ncomp,nsurf,nsurf_sec,jx,jy,jz)
  CALL JacPotentialLocal(ncomp,nsurf,nsurf_sec,npot,jx,jy,jz)

  CALL reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt,time)
  CALL jacmin(ncomp,nspec,nexchange,nsurf,nkin,nrct,jx,jy,jz,time,round)
    
  CALL reactkin(ncomp,nspec,nrct,ikin,jx,jy,jz,AqueousToBulk,time)
  CALL jacrkin(ncomp,nspec,nrct,ikin,jx,jy,jz,AqueousToBulk)

  CALL FxTopeLocal(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,nrct,  &
    nspec,ngas,neqn,delt,jx,jy,jz,nx,ny,nz,ne,time)

  CALL AssembleLocal(ncomp,nspec,nkin,nrct,ngas,ikin,  &
    nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,delt,jx,jy,jz,nx,ny,nz,ne,time)

!  Solve the sucker

  DO i = 1,neqn
    bb(i) = -fxx(i)
  END DO

!!  CALL ludcmp90(aaa,indd,det,neqn)
!!  CALL lubksb90(aaa,indd,bb,neqn)


  CALL dgetrf(neqn,neqn,aaa,neqn,indd,info)
  CALL dgetrs(trans,neqn,ione,aaa,neqn,indd,bb,neqn,info)

  errmax = 0.0D0
  DO i = 1,ncomp
    IF (DABS(bb(i)) > errmax) THEN
      errmax = DABS(bb(i))
    END IF
    IF (DABS(bb(i)) > corrmax) THEN
      bb(i) = SIGN(corrmax,bb(i))
    ELSE
      CONTINUE
    END IF
    sp(i,jx,jy,jz) = sp(i,jx,jy,jz) + bb(i)
    sp10(i,jx,jy,jz) = DEXP(sp(i,jx,jy,jz))
  END DO

  DO ix = 1,nexchange
    ind = ix + ncomp
    IF (DABS(bb(ind)) > errmax) THEN
      errmax = ABS(bb(ind))
    END IF
    IF (DABS(bb(ind)) > 1.0d0) THEN
      bb(ind) = SIGN(1.0d0,bb(ind))
    ELSE
      CONTINUE
    END IF
    spex(ix,jx,jy,jz) = spex(ix,jx,jy,jz) + bb(ind)
!!    spex10(ix,jx,jy,jz) = EXP(spex(ix,jx,jy,jz))
  END DO
  DO is = 1,nsurf
    ind =is+ncomp+nexchange
    IF (DABS(bb(ind)) > errmax) THEN
      errmax = DABS(bb(ind))
    END IF
    IF (DABS(bb(ind)) > corrmax) THEN
      bb(ind) = SIGN(corrmax,bb(ind))
    ELSE
      CONTINUE
    END IF
    spsurf(is,jx,jy,jz) = spsurf(is,jx,jy,jz) + bb(ind)
    spsurf10(is,jx,jy,jz) = DEXP(spsurf(is,jx,jy,jz))
  END DO
  DO npt = 1,npot
    ind = npt+ncomp+nexchange+nsurf
    IF (DABS(bb(ind)) > errmax) THEN
      errmax = DABS(bb(ind))
    END IF
    IF (DABS(bb(ind)) > 0.1d0) THEN
      bb(ind) = SIGN(0.1d0,bb(ind))
    ELSE
      CONTINUE
    END IF
    LogPotential(npt,jx,jy,jz) = LogPotential(npt,jx,jy,jz) + bb(ind)
  END DO

  IF (DABS(errmax) < 1.e-15) THEN
    icvg = 0
    EXIT
  END IF

  icvg = 0

  DO i = 1,ncomp
    ind = i
!!    tolmax = atol + rtol*sp10(i,jx,jy,jz)
    tolmax = atol
    IF (DABS(delt*fxx(ind)) > tolmax) THEN
      icvg = 1
      EXIT
    END IF 
  END DO

  DO ix = 1,nexchange
    tolmax = 1.D-06
    ind = ix+ncomp
    IF (DABS(fxx(ind)) > tolmax) THEN
      icvg = 1
      EXIT
    END IF
  END DO

  DO is = 1,nsurf
    tolmax = atol
    ind = is+ncomp+nexchange
    IF (DABS(delt*fxx(ind)) > tolmax) THEN
      icvg = 1
      EXIT
    END IF
  END DO

  DO npt = 1,npot
    tolmax = 1.e-07
     ind = npt+ncomp+nexchange+nsurf
     IF (DABS(fxx(ind)) > tolmax) THEN
       icvg = 1
       EXIT
     END IF
   END DO

  IF (icvg == 0) THEN
    EXIT
  END IF

END DO

IF (icvg == 1) THEN
  WRITE(*,*) ' No convergence at: ',jx,jy,jz
END IF

RETURN

END SUBROUTINE os3d_newton
