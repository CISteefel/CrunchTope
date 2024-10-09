SUBROUTINE PointByPointSolver(ncomp,nspec,nkin,nrct,ngas,ikin,       &
    nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,neqn,igamma,   & 
    delt,corrmax,jx,jy,jz,iterat,icvg,nx,ny,nz,time,AqueousToBulk)
USE crunchtype
USE runtime, ONLY: H2Opresent, Duan
USE params
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE temperature

IMPLICIT NONE

!!EXTERNAL dgetrf
!!EXTERNAL dgetrs


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

INTEGER(I4B), PARAMETER                       :: newton=50
INTEGER(I4B), PARAMETER                       :: ione=1
INTEGER(I4B)                                  :: ind
INTEGER(I4B)                                  :: ne
INTEGER(I4B)                                  :: i
INTEGER(I4B)                                  :: i2
INTEGER(I4B)                                  :: ix
INTEGER(I4B)                                  :: is
INTEGER(I4B)                                  :: npt
INTEGER(I4B)                                  :: info

REAL(DP), PARAMETER                           :: atol=1.e-08
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
DO ne = 1,newton
    
  iterat = iterat + 1

!! IF (igamma == 2) THEN
  IF (UpdateGammaInNewton) THEN
    IF (Duan) THEN
!!!      CALL gamma_co2(ncomp,nspec,ngas,jx,jy,jz)
    ELSE
      CALL gamma(ncomp,nspec,jx,jy,jz)
    END IF
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
  CALL jacmin(ncomp,nspec,nexchange,nsurf,nkin,nrct,jx,jy,jz,time)
    
  CALL reactkin(ncomp,nspec,nrct,ikin,jx,jy,jz,AqueousToBulk,time)
  CALL jacrkin(ncomp,nspec,nrct,ikin,jx,jy,jz,AqueousToBulk)

  CALL fx_local(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,nrct,  &
    nspec,ngas,neqn,delt,jx,jy,jz,nx,ny,nz,ne,time)

  CALL assemble_local(ncomp,nspec,nkin,nrct,ngas,ikin,  &
    nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,delt,jx,jy,jz,nx,ny,nz,ne,time)

!  Solve the sucker

  DO i = 1,neqn
!!    aascale = aaa(i,i)
    aascale = 1.0
    DO i2 = 1,neqn
      aaa(i,i2) = aaa(i,i2)/aascale
    END DO
    bb(i) = -fxx(i)/aascale
  END DO

!!  CALL ludcmp90(aaa,indd,det,neqn)
!!  CALL lubksb90(aaa,indd,bb,neqn)

  check = -fxx(1)/aaa(1,1)

  CALL dgetrf(neqn,neqn,aaa,neqn,indd,info)
  CALL dgetrs(trans,neqn,ione,aaa,neqn,indd,bb,neqn,info)

  errmax = 0.0D0
  DO i = 1,ncomp
    IF (ABS(bb(i)) > errmax) THEN
      errmax = ABS(bb(i))
    END IF
    IF (ABS(bb(i)) > corrmax) THEN
      bb(i) = SIGN(corrmax,bb(i))
    ELSE
      CONTINUE
    END IF
    sp(i,jx,jy,jz) = sp(i,jx,jy,jz) + bb(i)
    sp10(i,jx,jy,jz) = EXP(sp(i,jx,jy,jz))
  END DO
!!fp! compare_elem("sp end+ newton",{#expr# sp(1:ncomp,jx(1:nx),+jy(1:ny),jz(1:nz))#},0,{#expr#ne==1#},0,real8,20);

  DO ix = 1,nexchange
    ind = ix + ncomp
    IF (ABS(bb(ind)) > errmax) THEN
      errmax = ABS(bb(ind))
    END IF
    IF (ABS(bb(ind)) > 1.0d0) THEN
      bb(ind) = SIGN(1.0d0,bb(ind))
    ELSE
      CONTINUE
    END IF
    spex(ix,jx,jy,jz) = spex(ix,jx,jy,jz) + bb(ind)
!!    spex10(ix,jx,jy,jz) = EXP(spex(ix,jx,jy,jz))
  END DO
  DO is = 1,nsurf
    ind =is+ncomp+nexchange
    IF (ABS(bb(ind)) > errmax) THEN
      errmax = ABS(bb(ind))
    END IF
    IF (ABS(bb(ind)) > corrmax) THEN
      bb(ind) = SIGN(corrmax,bb(ind))
    ELSE
      CONTINUE
    END IF
    spsurf(is,jx,jy,jz) = spsurf(is,jx,jy,jz) + bb(ind)
    spsurf10(is,jx,jy,jz) = EXP(spsurf(is,jx,jy,jz))
  END DO
  DO npt = 1,npot
    ind = npt+ncomp+nexchange+nsurf
    IF (ABS(bb(ind)) > errmax) THEN
      errmax = ABS(bb(ind))
    END IF
    IF (ABS(bb(ind)) > 0.1) THEN
      bb(ind) = SIGN(0.1d0,bb(ind))
    ELSE
      CONTINUE
    END IF
    LogPotential(npt,jx,jy,jz) = LogPotential(npt,jx,jy,jz) + bb(ind)
  END DO

  IF (ABS(errmax) < 1.e-15) THEN
    icvg = 0
    EXIT
  END IF

  icvg = 0

  DO i = 1,ncomp
    ind = i
!!    tolmax = atol + rtol*sp10(i,jx,jy,jz)
    tolmax = atol
    IF (ABS(delt*fxx(ind)) > tolmax) THEN
      icvg = 1
      EXIT
    END IF 
  END DO

  DO ix = 1,nexchange
    tolmax = 1.e-06
    ind = ix+ncomp
    IF (ABS(fxx(ind)) > tolmax) THEN
      icvg = 1
      EXIT
    END IF
  END DO

  DO is = 1,nsurf
    tolmax = atol
    ind = is+ncomp+nexchange
    IF (ABS(delt*fxx(ind)) > tolmax) THEN
      icvg = 1
      EXIT
    END IF
  END DO

  DO npt = 1,npot
    tolmax = 1.e-07
     ind = npt+ncomp+nexchange+nsurf
     IF (ABS(fxx(ind)) > tolmax) THEN
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
