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

SUBROUTINE FxTopeGlobal(nx,ny,nz,ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,nrct,  &
    nspec,ngas,neqn,dt,jx,jy,jz,nBoundaryConditionZone)
USE crunchtype
USE params
USE runtime
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE temperature
USE flow
USE ReadFlow

IMPLICIT NONE

! ********** INTERFACE BLOCKS ******
INTERFACE
  SUBROUTINE GasInjection(ncomp,nspec,ngas,scorr,jx,jy,jz)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                           :: ncomp
  INTEGER(I4B), INTENT(IN)                           :: nspec
  INTEGER(I4B), INTENT(IN)                           :: ngas
  INTEGER(I4B), INTENT(IN)                           :: jz
  INTEGER(I4B), INTENT(IN)                           :: jx
  INTEGER(I4B), INTENT(IN)                           :: jy
  REAL(DP), DIMENSION(:), INTENT(OUT)                :: scorr
  END SUBROUTINE GasInjection
END INTERFACE
INTERFACE
  SUBROUTINE bdexchange(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,nbnd,sexb)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                      :: ncomp
  INTEGER(I4B), INTENT(IN)                      :: nspec
  INTEGER(I4B), INTENT(IN)                      :: nexchange
  INTEGER(I4B), INTENT(IN)                      :: nexch_sec
  INTEGER(I4B), INTENT(IN)                      :: nsurf
  INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
  INTEGER(I4B), INTENT(IN OUT)                  :: nbnd
  REAL(DP), DIMENSION(:), INTENT(OUT)           :: sexb
  END SUBROUTINE bdexchange
END INTERFACE
INTERFACE
  SUBROUTINE bdgas(ncomp,nspec,nrct,ngas,nbnd,scorr)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                      :: ncomp
  INTEGER(I4B), INTENT(IN)                      :: nspec
  INTEGER(I4B), INTENT(IN)                      :: nrct
  INTEGER(I4B), INTENT(IN)                      :: ngas
  INTEGER(I4B), INTENT(IN)                      :: nbnd
  REAL(DP), DIMENSION(:), INTENT(OUT)           :: scorr
  END SUBROUTINE bdgas
END INTERFACE
INTERFACE
  SUBROUTINE bdgas_by_grid(ncomp,nspec,nrct,ngas,jx,jy,jz,scorr)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                      :: ncomp
  INTEGER(I4B), INTENT(IN)                      :: nspec
  INTEGER(I4B), INTENT(IN)                      :: nrct
  INTEGER(I4B), INTENT(IN)                      :: ngas
  INTEGER(I4B), INTENT(IN)                      :: jx
  INTEGER(I4B), INTENT(IN)                      :: jy
  INTEGER(I4B), INTENT(IN)                      :: jz
  REAL(DP), DIMENSION(:), INTENT(OUT)           :: scorr
  END SUBROUTINE bdgas_by_grid
END INTERFACE
INTERFACE
  SUBROUTINE bdrecalc(ncomp,nspec,nbnd,scorr)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                      :: ncomp
  INTEGER(I4B), INTENT(IN)                      :: nspec
  INTEGER(I4B), INTENT(IN)                      :: nbnd
  REAL(DP), DIMENSION(:), INTENT(OUT)           :: scorr
  END SUBROUTINE bdrecalc
END INTERFACE
INTERFACE
SUBROUTINE bdrecalc_by_grid(ncomp,nspec,jx,jy,jz,scorr)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                      :: ncomp
  INTEGER(I4B), INTENT(IN)                      :: nspec
  INTEGER(I4B), INTENT(IN)                      :: jx
  INTEGER(I4B), INTENT(IN)                      :: jy
  INTEGER(I4B), INTENT(IN)                      :: jz
  REAL(DP), DIMENSION(:), INTENT(OUT)           :: scorr
  END SUBROUTINE bdrecalc_by_grid
END INTERFACE
INTERFACE
  SUBROUTINE bdsurf(ncomp,nsurf,nsurf_sec,nbnd,scorr)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                      :: ncomp
  INTEGER(I4B), INTENT(IN)                      :: nsurf
  INTEGER(I4B), INTENT(IN)                      :: nsurf_sec
  INTEGER(I4B), INTENT(IN)                      :: nbnd
  REAL(DP), DIMENSION(:), INTENT(OUT)           :: scorr
  END SUBROUTINE bdsurf
END INTERFACE
!  *******  END INTERFACE BLOCKS  *******

!fp! auto_par_loops=0;

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                  :: nx
INTEGER(I4B), INTENT(IN)                                  :: ny
INTEGER(I4B), INTENT(IN)                                  :: nz
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
INTEGER(I4B), INTENT(IN)                                  :: nBoundaryConditionZone

REAL(DP), INTENT(IN)                                      :: dt
!  Internal variables and arrays

REAL(DP)                                                  :: r
REAL(DP)                                                  :: satgas
REAL(DP)                                                  :: satgasOld
REAL(DP)                                                  :: df
REAL(DP)                                                  :: sumchgbd
REAL(DP)                                                  :: xvectors
REAL(DP)                                                  :: xvecgas
REAL(DP)                                                  :: xvec_ex
REAL(DP)                                                  :: xspecdiffw
REAL(DP)                                                  :: xspecdiffe
REAL(DP)                                                  :: xspecdiffs
REAL(DP)                                                  :: xspecdiffn
REAL(DP)                                                  :: xspecdiffw_tmp
REAL(DP)                                                  :: xspecdiffe_tmp
REAL(DP)                                                  :: yvectors
REAL(DP)                                                  :: yvecgas
REAL(DP)                                                  :: yvec_ex
REAL(DP)                                                  :: xbdflux
REAL(DP)                                                  :: ybdflux
REAL(DP)                                                  :: source
REAL(DP)                                                  :: aq_accum
REAL(DP)                                                  :: gas_accum
REAL(DP)                                                  :: gas_transport
REAL(DP)                                                  :: ex_accum
REAL(DP)                                                  :: ex_transport
REAL(DP)                                                  :: surf_accum
REAL(DP)                                                  :: xvec_surf
REAL(DP)                                                  :: yvec_surf
REAL(DP)                                                  :: surf_transport
REAL(DP)                                                  :: portemp
REAL(DP)                                                  :: rotemp
REAL(DP)                                                  :: rotempOld
REAL(DP)                                                  :: satl
REAL(DP)                                                  :: satlOld
REAL(DP)                                                  :: CellVolume
REAL(DP)                                                  :: recharge
REAL(DP)                                                  :: pi
REAL(DP)                                                  :: MultiplyCell

INTEGER(I4B)                                              :: npz
INTEGER(I4B)                                              :: nbnd
INTEGER(I4B)                                              :: jdum
INTEGER(I4B)                                              :: jdum2
INTEGER(I4B)                                              :: jz
INTEGER(I4B)                                              :: j
INTEGER(I4B)                                              :: i
INTEGER(I4B)                                              :: is
INTEGER(I4B)                                              :: ix
INTEGER(I4B)                                              :: nb
INTEGER(I4B)                                              :: ind

REAL(DP)                                                  :: GasSource
REAL(DP)                                                  :: Retardation
REAL(DP)                                                  :: check
REAL(DP)                                                  :: A_transpi

  !!! Added July 17 by Carl (hopefully not stomped on)

!  This routine calculates (and assembles) the function residuals.

surf_accum = 0.0d0
aq_accum = 0.0d0

jz = 1
j = (jz-1)*nx*ny + (jy-1)*nx + jx

Retardation = 0.001d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))/por(jx,jy,jz)

pi = DACOS(-1.0d0)
portemp = por(jx,jy,jz)
rotemp = ro(jx,jy,jz)
rotempOld = roOld(jx,jy,jz)
satl = satliq(jx,jy,jz)
satlOld = satliqOld(jx,jy,jz)
satgas = 1.0 - satl
satgasOld = 1.0 - satlOld

r = 1.0D0/dt
IF (cylindrical) THEN
  CellVolume = dyy(jy)*pi*( (x(jx)+dxx(jx)/2.0d0 )**2.0d0 - ( x(jx)-dxx(jx)/2.0d0 )**2.0d0  )
  df = 1.0d0
  MultiplyCell = CellVolume
  IF (CylindricalDivideVolume) THEN
      df = 1.0/CellVolume
      MultiplyCell = 1.0
  END IF
ELSE IF (spherical) THEN
  CellVolume = (4.0d0/3.0d0)*pi*( (x(jx)+0.5d0*dxx(jx))**3.0d0 - (x(jx)-0.5d0*dxx(jx))**3.0d0  )
  df = 1.0d0
  MultiplyCell = CellVolume
  IF (CylindricalDivideVolume) THEN
      df = 1.0/CellVolume
      MultiplyCell = 1.0
  END IF
!!  df = 1.0/CellVolume
!!  MultiplyCell = 1.0
ELSE
      CellVolume = dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
       df = 1.0d0
       MultiplyCell = CellVolume
!!      df = 1.0/CellVolume
!!      MultiplyCell = 1.0
END IF

fxmax = 0.0d0 

IF (species_diffusion .and. nBoundaryConditionZone == 0) THEN

  nbnd = 1
  CALL bd_diffuse(ncomp,nspec,nbnd,sumchgbd)
  DO i = 1,ncomp
    s_dsp(i,0,jy,jz) = sdsp(i)
    s_chg(i,0,jy,jz) = schg(i)
  END DO
  sumwtchg(0,jy,jz) = sumchgbd
  
  nbnd = 2
  CALL bd_diffuse(ncomp,nspec,nbnd,sumchgbd)
  DO i = 1,ncomp
    s_dsp(i,nx+1,jy,jz) = sdsp(i)
    s_chg(i,nx+1,jy,jz) = schg(i)
  END DO
  sumwtchg(nx+1,jy,jz) = sumchgbd

  nbnd = 3
  CALL bd_diffuse(ncomp,nspec,nbnd,sumchgbd)
  DO i = 1,ncomp
    s_dsp(i,jx,0,jz) = sdsp(i)
    s_chg(i,jx,0,jz) = schg(i)
  END DO
  sumwtchg(jx,0,jz) = sumchgbd

  nbnd = 4
  CALL bd_diffuse(ncomp,nspec,nbnd,sumchgbd)
  DO i = 1,ncomp
    s_dsp(i,jx,ny+1,jz) = sdsp(i)
    s_chg(i,jx,ny+1,jz) = schg(i)
  END DO
  sumwtchg(jx,ny+1,jz) = sumchgbd

END IF

IF (ierode /= 1) THEN
  IF (nsurf > 0 .OR. nexchange > 0) THEN
    CALL totexchange_local(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)
  END IF
END IF
IF (ierode /= 1 .AND. nsurf > 0) THEN
  CALL totsurf_local(ncomp,nsurf,nsurf_sec,jx,jy,jz)
END IF

IF (activecell(jx,jy,jz) == 0) GO TO 800

IF (nx == 1) GO TO 100

IF (jx == 1) THEN

  jdum = jx+1
  DO i = 1,ncomp
    sce(i) = s(i,jdum,jy,jz)*xgram(jdum,jy,jz)
  END DO

  IF (nBoundaryConditionZone > 0) THEN   !! Boundary cells by grid

!!!  BC by grid cells

    jdum2 = 0
    IF (isaturate == 1) THEN
      DO i = 1,ncomp
        sge(i) = sgas(i,jdum,jy,jz)
      END DO
      CALL bdgas_by_grid(ncomp,nspec,nrct,ngas,jdum2,jy,jz,sgw)
    END IF

    IF (JcByGrid(jx-1,jy,jz) == 1 .OR. netflowx(jdum2,jy,jz) > 0.0) THEN
      CALL bdrecalc_by_grid(ncomp,nspec,jdum2,jy,jz,scw)
    END IF

  ELSE    !!  Conventional face treatment of BC

    nbnd = 1
    IF (isaturate == 1) THEN
      DO i = 1,ncomp
        sge(i) = sgas(i,jdum,jy,jz)
      END DO
      CALL bdgas(ncomp,nspec,nrct,ngas,nbnd,sgw)
    END IF
    IF (jc(1) == 1 .OR. netflowx(0,jy,jz) > 0.0) THEN
      CALL bdrecalc(ncomp,nspec,nbnd,scw)
    END IF

  END IF
 
  IF (ierode == 1) THEN
    DO i = 1,ncomp
      sex_east(i) = sch(i,jdum,jy,jz)
    END DO
    IF (SolidBuryX(0) > 0.0) THEN
      CALL bdexchange(ncomp,nspec,nexchange,nexch_sec,  &
          nsurf,nsurf_sec,nbnd,sex_west)
    ELSE
      DO i = 1,ncomp
        sex_west(i) = 0.0
      END DO
    END IF
    DO is = 1,nsurf
      surf_east(is) = ssurf(is,jdum,jy,jz)
    END DO
    IF (SolidBuryX(0) > 0.0) THEN
      CALL bdsurf(ncomp,nsurf,nsurf_sec,nbnd,surf_west)
    ELSE
      DO is = 1,nsurf
        surf_west(is) = 0.0
      END DO
    END IF
  END IF    !  Block only used if ierode = 1 (erosion or burial)

ELSE IF (jx == nx) THEN

  jdum = jx-1
  DO i = 1,ncomp
    scw(i) = s(i,jdum,jy,jz)*xgram(jdum,jy,jz)
  END DO

  IF (nBoundaryConditionZone > 0) THEN   !! Boundary cells by grid

!!!  BC by grid cells

    jdum2 = nx+1
    IF (isaturate == 1) THEN
      DO i = 1,ncomp
        sgw(i) = sgas(i,jdum,jy,jz)
      END DO
      CALL bdgas_by_grid(ncomp,nspec,nrct,ngas,jdum2,jy,jz,sge)
    END IF

    IF (JcByGrid(jx+1,jy,jz) == 1 .OR. netflowx(nx,jy,jz) < 0.0) THEN
      CALL bdrecalc_by_grid(ncomp,nspec,jdum2,jy,jz,sce)
    END IF

  ELSE    !!  Conventional face treatment of BC

    nbnd = 2
    IF (isaturate == 1) THEN
      DO i = 1,ncomp
        sgw(i) = sgas(i,jdum,jy,jz)
      END DO
      CALL bdgas(ncomp,nspec,nrct,ngas,nbnd,sge)
    END IF
    IF (jc(2) == 1 .OR. netflowx(nx,jy,jz) < 0.0) THEN
      CALL bdrecalc(ncomp,nspec,nbnd,sce)
            continue
    END IF

  END IF

  IF (ierode == 1) THEN
    DO i = 1,ncomp
      sex_west(i) = sch(i,jdum,jy,jz)
    END DO
    IF (SolidBuryX(nx) < 0.0) THEN
      CALL bdexchange(ncomp,nspec,nexchange,nexch_sec,  &
          nsurf,nsurf_sec,nbnd,sex_east)
    ELSE
      DO i = 1,ncomp
        sex_east(i) = 0.0
      END DO
    END IF
    DO is = 1,nsurf
      surf_west(is) = ssurf(is,jdum,jy,jz)
    END DO
    IF (SolidBuryX(nx) < 0.0) THEN
      CALL bdsurf(ncomp,nsurf,nsurf_sec,nbnd,surf_east)
    ELSE
      DO is = 1,nsurf
        surf_east(is) = 0.0
      END DO
    END IF
  END IF    !  Block used only if ierode = 1 (erosion or burial)

ELSE

  jdum = jx+1
  DO i = 1,ncomp
    sce(i) = s(i,jdum,jy,jz)*xgram(jdum,jy,jz)
  END DO
  IF (isaturate == 1) THEN
    DO i = 1,ncomp
      sge(i) = sgas(i,jdum,jy,jz)
    END DO
  END IF

  IF (ierode == 1) THEN
    DO i = 1,ncomp
      sex_east(i) = sch(i,jdum,jy,jz)
    END DO
    DO is = 1,nsurf
      surf_east(is) = ssurf(is,jdum,jy,jz)
    END DO
  END IF    !  Block only used if ierode = 1 (erosion or burial)

  jdum = jx-1
  DO i = 1,ncomp
    scw(i) = s(i,jdum,jy,jz)*xgram(jdum,jy,jz)
  END DO
  IF (isaturate == 1) THEN
    DO i = 1,ncomp
      sgw(i) = sgas(i,jdum,jy,jz)
    END DO
  END IF

  IF (ierode == 1) THEN
    DO i = 1,ncomp
      sex_west(i) = sch(i,jdum,jy,jz)
    END DO
    DO is = 1,nsurf
      surf_west(is) = ssurf(is,jdum,jy,jz)
    END DO
  END IF    !  Block only used if ierode = 1 (erosion or burial)

END IF

100 CONTINUE
IF (ny == 1) GO TO 200

IF (jy == 1) THEN

  jdum = jy+1
  DO i = 1,ncomp
    scn(i) = s(i,jx,jdum,jz)*xgram(jx,jdum,jz)
  END DO

 IF (nBoundaryConditionZone > 0) THEN   !! Boundary cells by grid

!!!  BC by grid cells

    jdum2 = 0
    IF (isaturate == 1) THEN
      DO i = 1,ncomp
        sgn(i) = sgas(i,jdum,jy,jz)
      END DO
      CALL bdgas_by_grid(ncomp,nspec,nrct,ngas,jx,jdum2,jz,sgs)
    END IF

!!!    IF (jc(3) == 1 .OR. qy(jx,0,jz) > 0.0) THEN
    IF (JcByGrid(jx,jy-1,jz) == 1) THEN
      CALL bdrecalc_by_grid(ncomp,nspec,jx,jdum2,jz,scs)
    END IF

  ELSE    !!  Conventional face treatment of BC

    nbnd = 3
    IF (isaturate == 1) THEN
      DO i = 1,ncomp
        sgn(i) = sgas(i,jx,jdum,jz)
      END DO
      CALL bdgas(ncomp,nspec,nrct,ngas,nbnd,sgs)
    END IF
    IF (jc(3) == 1 .OR. qy(jx,0,jz) > 0.0) THEN
      CALL bdrecalc(ncomp,nspec,nbnd,scs)
    END IF

  END IF

  IF (ierode == 1) THEN
    DO i = 1,ncomp
      sex_north(i) = sch(i,jx,jdum,jz)
    END DO
    IF (SolidBuryY(0) > 0.0) THEN
      CALL bdexchange(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,nbnd,sex_south)
    ELSE
      DO i = 1,ncomp
        sex_south(i) = 0.0
      END DO
    END IF
    DO is = 1,nsurf
      surf_north(is) = ssurf(is,jx,jdum,jz)
    END DO
    IF (SolidBuryY(0) > 0.0) THEN
      CALL bdsurf(ncomp,nsurf,nsurf_sec,nbnd,surf_south)
    ELSE
      DO is = 1,nsurf
        surf_south(is) = 0.0
      END DO
    END IF
  END IF    !   Block only used if ierode = 1 (erosion or burial)

ELSE IF (jy == ny) THEN

  jdum = jy-1
  DO i = 1,ncomp
    scs(i) = s(i,jx,jdum,jz)*xgram(jx,jdum,jz)
  END DO

IF (nBoundaryConditionZone > 0) THEN   !! Boundary cells by grid

!!!  BC by grid cells

    jdum2 = ny+1
    IF (isaturate == 1) THEN
      DO i = 1,ncomp
        sgs(i) = sgas(i,jx,jdum2,jz)
      END DO
      CALL bdgas_by_grid(ncomp,nspec,nrct,ngas,jx,jdum2,jz,sgn)
    END IF

!!!    IF (jc(4) == 1 .OR. qy(jx,ny,jz) < 0.0) THEN
    IF (JcByGrid(jx,jy+1,jz) == 1) THEN
      CALL bdrecalc_by_grid(ncomp,nspec,jx,jdum2,jz,scn)
    END IF

  ELSE    !!  Conventional face treatment of BC

    nbnd = 4
    IF (isaturate == 1) THEN
      DO i = 1,ncomp
        sgs(i) = sgas(i,jx,jdum,jz)
      END DO
      CALL bdgas(ncomp,nspec,nrct,ngas,nbnd,sgn)
    END IF
    IF (jc(4) == 1 .OR. qy(jx,ny,jz) < 0.0) THEN
      CALL bdrecalc(ncomp,nspec,nbnd,scn)
    END IF

  END IF

  IF (ierode == 1) THEN
    DO i = 1,ncomp
      sex_south(i) = sch(i,jx,jdum,jz)
    END DO
    IF (SolidBuryY(ny) < 0.0) THEN
      CALL bdexchange(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,nbnd,sex_north)
    ELSE
      DO i = 1,ncomp
        sex_north(i) = 0.0
      END DO
    END IF
    DO is = 1,nsurf
      surf_south(is) = ssurf(is,jx,jdum,jz)
    END DO
    IF (SolidBuryY(ny) < 0.0) THEN
      CALL bdsurf(ncomp,nsurf,nsurf_sec,nbnd,surf_north)
    ELSE
      DO is = 1,nsurf
        surf_north(is) = 0.0
      END DO
    END IF
  END IF     !  Block only used if ierode = 1 (erosion or burial)

ELSE
  jdum = jy+1
  DO i = 1,ncomp
    scn(i) = s(i,jx,jdum,jz)*xgram(jx,jdum,jz)
  END DO
  IF (isaturate == 1) THEN
    DO i = 1,ncomp
      sgn(i) = sgas(i,jx,jdum,jz)
    END DO
  END IF
  IF (ierode == 1) THEN
    DO i = 1,ncomp
      sex_north(i) = sch(i,jx,jdum,jz)
    END DO
    DO is = 1,nsurf
      surf_north(is) = ssurf(is,jx,jdum,jz)
    END DO
  END IF       !  Block only used if ierode = 1 (erosion or burial)

  jdum = jy-1
  DO i = 1,ncomp
    scs(i) = s(i,jx,jdum,jz)*xgram(jx,jdum,jz)
  END DO
  IF (isaturate == 1) THEN
    DO i = 1,ncomp
      sgs(i) = sgas(i,jx,jdum,jz)
    END DO
  END IF
  IF (ierode == 1) THEN
    DO i = 1,ncomp
      sex_south(i) = sch(i,jx,jdum,jz)
    END DO
    DO is = 1,nsurf
      surf_south(is) = ssurf(is,jx,jdum,jz)
    END DO
  END IF        !  Block only used if ierode = 1 (erosion or burial)

END IF

200 IF (species_diffusion) THEN
  dgradw = 0.0
  dgrade = 0.0
  dgrads = 0.0
  dgradn = 0.0
  DO i = 1,ncomp
    dgradw = dgradw + chg(i)*a_d(jx,jy,jz)*(s_dsp(i,jx-1,jy,jz) - s_dsp(i,jx,jy,jz))
    dgrade = dgrade + chg(i)*c_d(jx,jy,jz)*(s_dsp(i,jx+1,jy,jz) - s_dsp(i,jx,jy,jz))
    dgrads = dgrads + chg(i)*f_d(jx,jy,jz)*(s_dsp(i,jx,jy-1,jz) - s_dsp(i,jx,jy,jz))
    dgradn = dgradn + chg(i)*d_d(jx,jy,jz)*(s_dsp(i,jx,jy+1,jz) - s_dsp(i,jx,jy,jz))
    sigma_w(i) = 0.5*(dstar(jx,jy,jz)*s_chg(i,jx,jy,jz) + dstar(jx-1,jy,jz)*s_chg(i,jx-1,jy,jz))
    sigma_e(i) = 0.5*(dstar(jx,jy,jz)*s_chg(i,jx,jy,jz) + dstar(jx+1,jy,jz)*s_chg(i,jx+1,jy,jz))
    sigma_s(i) = 0.5*(dstar(jx,jy,jz)*s_chg(i,jx,jy,jz) + dstar(jx,jy-1,jz)*s_chg(i,jx,jy-1,jz))
    sigma_n(i) = 0.5*(dstar(jx,jy,jz)*s_chg(i,jx,jy,jz) + dstar(jx,jy+1,jz)*s_chg(i,jx,jy+1,jz))
  END DO
  sumsigma_w = 0.5*(dstar(jx,jy,jz)*sumwtchg(jx,jy,jz) + dstar(jx-1,jy,jz)*sumwtchg(jx-1,jy,jz))
  sumsigma_e = 0.5*(dstar(jx,jy,jz)*sumwtchg(jx,jy,jz) + dstar(jx+1,jy,jz)*sumwtchg(jx+1,jy,jz))
  sumsigma_s = 0.5*(dstar(jx,jy,jz)*sumwtchg(jx,jy,jz) + dstar(jx,jy-1,jz)*sumwtchg(jx,jy-1,jz))
  sumsigma_n = 0.5*(dstar(jx,jy,jz)*sumwtchg(jx,jy,jz) + dstar(jx,jy+1,jz)*sumwtchg(jx,jy+1,jz))
END IF

!!!IF (isaturate == 1) THEN
!!!  IF (gaspump(1,jx,jy,jz) /= 0.0) THEN
!!!    call GasInjection(ncomp,nspec,ngas,sgaspump,jx,jy,jz)
!!!  END IF
!!!END IF

DO i = 1,ncomp
  
  ind = (j-1)*(neqn) + i
  
  IF (nx == 1) GO TO 300
  
  IF (jx == 1) THEN
    
    IF (jc(1) == 1 .or. JcByGrid(jx-1,jy,jz) == 1 ) THEN    ! Dirichlet bdy
      xvectors = a(jx,jy,jz)*scw(i) + c(jx,jy,jz)*sce(i)
      IF (isaturate == 1) THEN
        xvecgas = ag(jx,jy,jz)*sgw(i) + cg(jx,jy,jz)*sge(i)
      ELSE
        xvecgas = 0.0
      END IF
    ELSE
      xvectors = c(jx,jy,jz)*sce(i)
      IF (isaturate == 1) THEN
        xvecgas = ag(jx,jy,jz)*sgw(i) + cg(jx,jy,jz)*sge(i)
      ELSE
        xvecgas = 0.0
      END IF
    END IF
    
  ELSE IF (jx == nx) THEN
    
    IF (jc(2) == 1 .or. JcByGrid(jx+1,jy,jz) == 1) THEN    ! Dirichlet bdy
      xvectors = a(jx,jy,jz)*scw(i) + c(jx,jy,jz)*sce(i)
      IF (isaturate == 1) THEN
        xvecgas = ag(jx,jy,jz)*sgw(i) + cg(jx,jy,jz)*sge(i)
      ELSE
        xvecgas = 0.0
      END IF
    ELSE      
      xvectors = a(jx,jy,jz)*scw(i)
      IF (isaturate == 1) THEN
        xvecgas = ag(jx,jy,jz)*sgw(i) + cg(jx,jy,jz)*sge(i)
      ELSE
        xvecgas = 0.0
      END IF
    END IF
    
  ELSE
    
    xvectors = a(jx,jy,jz)*scw(i) + c(jx,jy,jz)*sce(i)
    IF (isaturate == 1) THEN
      xvecgas = ag(jx,jy,jz)*sgw(i) + cg(jx,jy,jz)*sge(i)
    ELSE
      xvecgas = 0.0
    END IF
    
  END IF
  
  IF (ierode == 1) THEN
    xvec_ex = cbu(jx,jy,jz)*sex_east(i) + abu(jx,jy,jz)*sex_west(i)
  ELSE
    xvec_ex = 0.0
  END IF
  
  IF (species_diffusion) THEN
     
    xspecdiffw = a_d(jx,jy,jz)*(s_dsp(i,jx-1,jy,jz) - s_dsp(i,jx,jy,jz))  &
        - (sigma_w(i)/sumsigma_w) * dgradw
    xspecdiffe = c_d(jx,jy,jz)*(s_dsp(i,jx+1,jy,jz) - s_dsp(i,jx,jy,jz))  &
        - (sigma_e(i)/sumsigma_e) * dgrade

    xspecdiffs = f_d(jx,jy,jz)*(s_dsp(i,jx,jy-1,jz) - s_dsp(i,jx,jy,jz))  &
        - (sigma_s(i)/sumsigma_s) * dgrads
    xspecdiffn = d_d(jx,jy,jz)*(s_dsp(i,jx,jy+1,jz) - s_dsp(i,jx,jy,jz))  &
        - (sigma_n(i)/sumsigma_n) * dgradn


!!    IF (Migration) THEN
!!      xspecdiffn = d_d(jx,jy,jz)*(s_dsp(i,jx,jy+1,jz) - s_dsp(i,jx,jy,jz)) & 
!!        - (sigma_n(i)/sumsigma_n) * dgradn
!!    ELSE
!!      xspecdiffn = d_d(jx,jy,jz)*(s_dsp(i,jx,jy+1,jz) - s_dsp(i,jx,jy,jz))  
!!        - (sigma_n(i)/sumsigma_n) * dgradn
!!    END IF

!    xspecdiffw = a_d(jx,jy,jz)*(s_dsp(i,jx-1,jy,jz) - s_dsp(i,jx,jy,jz))  
!    xspecdiffe = c_d(jx,jy,jz)*(s_dsp(i,jx+1,jy,jz) - s_dsp(i,jx,jy,jz))  

  ELSE
    
    xspecdiffw = 0.0
    xspecdiffe = 0.0
    xspecdiffn = 0.0
    xspecdiffs = 0.0
    
  END IF
  

! ************************
  
  300   CONTINUE
  
  IF (ny == 1) GO TO 400
  
  IF (jy == 1) THEN
!!!    IF (jc(3) == 1) THEN    ! Dirichlet bdy
    IF (jc(3) == 1 .or. JcByGrid(jx,jy-1,jz) == 1) THEN    ! Dirichlet bdy
      yvectors = d(jx,jy,jz)*scn(i) + f(jx,jy,jz)*scs(i)
      IF (isaturate == 1) THEN
        yvecgas =  dg(jx,jy,jz)*sgn(i) + fg(jx,jy,jz)*sgs(i)
      ELSE
        yvecgas = 0.0
      END IF
    ELSE
      yvectors = d(jx,jy,jz)*scn(i)
      IF (isaturate == 1) THEN
        yvecgas =  dg(jx,jy,jz)*sgn(i) + fg(jx,jy,jz)*sgs(i)
 !!       yvecgas =  dg(jx,jy,jz)*sgn(i)
      ELSE
        yvecgas = 0.0
      END IF
    END IF
  ELSE IF (jy == ny) THEN
!!!    IF (jc(4) == 1) THEN
    IF (jc(4) == 1 .or. JcByGrid(jx,jy+1,jz) == 1) THEN    ! Dirichlet bdy
      yvectors = d(jx,jy,jz)*scn(i) + f(jx,jy,jz)*scs(i)
      IF (isaturate == 1) THEN
        yvecgas =  dg(jx,jy,jz)*sgn(i) + fg(jx,jy,jz)*sgs(i)
      ELSE
        yvecgas = 0.0
      END IF
    ELSE
      yvectors = f(jx,jy,jz)*scs(i)
      IF (isaturate == 1) THEN
        yvecgas =  dg(jx,jy,jz)*sgn(i) + fg(jx,jy,jz)*sgs(i)
!!        yvecgas =  fg(jx,jy,jz)*sgs(i)
      ELSE
        yvecgas = 0.0
      END IF
    END IF
    
  ELSE
    
    yvectors = d(jx,jy,jz)*scn(i) + f(jx,jy,jz)*scs(i)   
    
    IF (isaturate == 1) THEN
      yvecgas =  dg(jx,jy,jz)*sgn(i) + fg(jx,jy,jz)*sgs(i)
    ELSE
      yvecgas = 0.0
    END IF
    
  END IF
  
  IF (ierode == 1) THEN
    yvec_ex = dbu(jx,jy,jz)*sex_north(i) + fbu(jx,jy,jz)*sex_south(i)
  ELSE
    yvec_ex = 0.0
  END IF
  
  400   CONTINUE
  
  IF (nx == 1) GO TO 500
  
  IF (jx == 1 .AND. netflowx(0,jy,jz) > 0.0) THEN
    IF (jc(1) == 2 .or. JcByGrid(jx-1,jy,jz) /= 1) THEN
      xbdflux = a(jx,jy,jz)*scw(i)
    END IF
  ELSE IF (jx == nx .AND. netflowx(nx,jy,jz) < 0.0) THEN
    IF (jc(2) == 2 .or. JcByGrid(jx+1,jy,jz) /= 1) THEN
      xbdflux = c(jx,jy,jz)*sce(i)
    END IF
  ELSE
    xbdflux = 0.0
  END IF
  
  500   CONTINUE
  IF (ny == 1) GO TO 600
  
  IF (jy == 1 .AND. qy(jx,0,jz) > 0.0) THEN
    IF (jc(3) == 2 .or. JcByGrid(jx,jy-1,jz) /= 1) THEN
      ybdflux = f(jx,jy,jz)*scs(i)
    END IF
  ELSE IF (jy == ny .AND. qy(jx,ny,jz) < 0.0) THEN
    IF (jc(4) == 2 .or. JcByGrid(jx,jy+1,jz) /= 1) THEN
      ybdflux = d(jx,jy,jz)*scn(i)
    END IF
  ELSE
    ybdflux = 0.0
  END IF
  
600 CONTINUE
  
!!NOTE:  "GIMRT" source term in m**3/year
  source = 0.0d0
  IF (wells) THEN
   
    DO npz = 1,npump(jx,jy,jz)
    
      IF (qg(npz,jx,jy,jz) > 0.0) THEN    ! Injection well

        source = source + xgram(jx,jy,jz)*qg(npz,jx,jy,jz)*rotemp*scond(i,intbnd(npz,jx,jy,jz))/CellVolume

      ELSE IF (qg(npz,jx,jy,jz) < 0.0) THEN    ! Pumping well

        source = source + xgram(jx,jy,jz)*qg(npz,jx,jy,jz)*rotemp*s(i,jx,jy,jz)/CellVolume
        
      ELSE
        CONTINUE
      END IF
    
    END DO
    
  END IF

! ************************************
! Edit by Lucien Stolze, June 2023
! Extract solutes via transpiration
  IF ((transpifix .OR. transpitimeseries) .AND. Richards) THEN
        if (ny == 1 .AND. nz == 1) THEN
        A_transpi = dyy(jy) * dzz(jx,jy,jz)
        source = source - xgram(jx,jy,jz)*transpirate_cell(jx)*A_transpi*rotemp*s(i,jx,jy,jz)/CellVolume
  ENDIF
  ENDIF
! ************************************
! end of edit by Lucien Stolze, June 2023

  GasSource = 0.0

!!!  IF (isaturate == 1) THEN
!!!    IF (gaspump(1,jx,jy,jz) /= 0.0) THEN
!!!      IF (cylindrical .OR. spherical) THEN
!!        GasSource = gaspump(jx,jy,jz)*sgaspump(i)
!!!        GasSource = gaspump(1,jx,jy,jz)*sgaspump(i)/CellVolume
!!!      ELSE
!!!        GasSource = gaspump(1,jx,jy,jz)*sgaspump(i)/CellVolume
!!!      END IF
!!!    ELSE
!!!      GasSource = 0.0
!!!    END IF
!!!  END IF

  IF (jy == 1) THEN
    IF (ReadNuft .AND. infiltration /= 0) THEN
      recharge = xgram(jx,jy,jz)*qrecharge(jx,jy)*scond(i,infiltration)/dyy(jy)
    ELSE
      recharge = 0.0
    END IF
  ELSE
    recharge = 0.0
  END IF
  
  IF (i /= ikh2o) THEN
    aq_accum = xgram(jx,jy,jz)*r*portemp*                            &
      (H2Oreacted(jx,jy,jz)*rotemp*satl*s(i,jx,jy,jz) -            &
       rotempOld*satlOld*sn(i,jx,jy,jz))*(1.0 + Retardation*distrib(i) )


  ELSE
    aq_accum = xgram(jx,jy,jz)*r*portemp*                            &
      (rotemp*satl*s(i,jx,jy,jz) -            &
       rotempOld*satlOld*sn(i,jx,jy,jz))*(1.0 + Retardation*distrib(i) )
  END IF
  
  IF (isaturate == 1) THEN
!!!    gas_accum = portemp*r*(satgas*sgas(i,jx,jy,jz) - satgasOld*sgasn(i,jx,jy,jz))
        gas_accum = portemp*r*satgas*(sgas(i,jx,jy,jz) - sgasn(i,jx,jy,jz))
  ELSE
    gas_accum = 0.0
    gas_transport = 0.0
  END IF
  
  IF (ierode == 1) THEN
    ex_accum = r*(sch(i,jx,jy,jz)-sexOld(i,jx,jy,jz)- ssurfOld(i,jx,jy,jz))
  ELSE
    ex_accum = r*(sNCexch_local(i)+sNCsurf_local(i) -sexOld(i,jx,jy,jz) - ssurfOld(i,jx,jy,jz)) 
    ex_transport = 0.0
  END IF
  
  IF (nx == 1 .AND. ny == 1) THEN

    fxx(ind) = MultiplyCell*(aq_accum + gas_accum + ex_accum - source - GasSource)


  ELSE IF (nx == 1) THEN
    
    IF (ierode == 1) THEN
      ex_transport =  df*ebu(jx,jy,jz)*sch(i,jx,jy,jz)
    END IF
    IF (isaturate == 1) THEN
      gas_transport = df*eg(jx,jy,jz)*sgas(i,jx,jy,jz)
    END IF
    
    fxx(ind) = MultiplyCell*(aq_accum + gas_accum + ex_accum - recharge - source - GasSource )&
        + xgram(jx,jy,jz)*df*e(jx,jy,jz)*s(i,jx,jy,jz) + yvectors*df  &
        + df*ybdflux + yvec_ex*df  &
        + yvecgas*df + ex_transport  &
        + gas_transport 
    
  ELSE IF (ny == 1) THEN  
    
    IF (ierode == 1) THEN
      ex_transport =  df*bbu(jx,jy,jz)*sch(i,jx,jy,jz)
    END IF
    IF (isaturate == 1) THEN
      gas_transport = df*bg(jx,jy,jz)*sgas(i,jx,jy,jz)
    END IF
    
    fxx(ind) = MultiplyCell*(aq_accum + gas_accum + ex_accum - recharge - source - GasSource)  &
              + xgram(jx,jy,jz)*df*b(jx,jy,jz)*s(i,jx,jy,jz)   &  !! Diagonal aqueous transport
              + xvectors*df     &   !! Off-diagonal aqueous transport
              + df*xbdflux      &   !! Advective flux through boundary
              + xvec_ex*df      &   !! Erosion flux of exchangers
              + ex_transport    &   !! exchanger burial or transport
              + xvecgas*df      &   !! off-diagonal gas transport
              + gas_transport   &   !! diagonal gas transport
              + xspecdiffw*df   &   !! Species-dependent diffusion
              + xspecdiffe*df       !! Species-dependent diffusion
    
      CONTINUE
      
      if (jx == 1) then
        continue
      end if
      
  ELSE
      
    IF (ierode == 1) THEN
      ex_transport =  df*sch(i,jx,jy,jz)* ( bbu(jx,jy,jz) + ebu(jx,jy,jz) )
    END IF
    IF (isaturate == 1) THEN
      gas_transport = df*sgas(i,jx,jy,jz)* ( bg(jx,jy,jz)+ eg(jx,jy,jz) )
    END IF
    fxx(ind) = MultiplyCell*(aq_accum + gas_accum + ex_accum - recharge - source - GasSource ) &
        + xgram(jx,jy,jz)*df*b(jx,jy,jz)*s(i,jx,jy,jz) + xgram(jx,jy,jz)*df*e(jx,jy,jz)*s(i,jx,jy,jz)  &
        + xvectors*df + yvectors*df + df*xbdflux + df*ybdflux  &
        + xvec_ex*df + yvec_ex*df + xvecgas*df + yvecgas*df  &
        + ex_transport + gas_transport  &
        + xspecdiffw*df + xspecdiffe*df + xspecdiffs*df + xspecdiffn*df ! Species-dependent diffusion
        
  END IF
  
  IF (DABS(fxx(ind)) > fxmax(i)) THEN
    fxmax(i) = DABS(fxx(ind))
  END IF
  
END DO

DO is = 1,nsurf
  ind = (j-1)*(neqn) + is+ncomp+nexchange
  
  IF (ierode == 1) THEN  
    surf_accum = r*(ssurf(is,jx,jy,jz)-ssurfn(is,jx,jy,jz))
    xvec_surf = cbu(jx,jy,jz)*surf_east(is) + abu(jx,jy,jz)*surf_west(is)
    yvec_surf = dbu(jx,jy,jz)*surf_north(is) + fbu(jx,jy,jz)*surf_south(is)
  ELSE
    check = ssurf_local(is) - ssurfn(is,jx,jy,jz)
    surf_accum = r*check
!!!    surf_accum = r*(ssurf_local(is)-ssurfn(is,jx,jy,jz))
!!!    surf_accum = r*( ssurf_local(is) - c_surf(is,jinit(jx,jy,jz)) )
    surf_accum = 0.0d0
    xvec_surf = 0.0
    yvec_surf = 0.0
  END IF
  
  IF (nx == 1 .AND. ny == 1) THEN
    
    fxx(ind) = MultiplyCell*surf_accum
    
  ELSE IF (nx == 1) THEN
    
    IF (ierode == 1) THEN
      surf_transport = df*ebu(jx,jy,jz)*ssurf(is,jx,jy,jz)
    ELSE
      surf_transport = 0.0
    END IF
    fxx(ind) = MultiplyCell*(surf_accum) + surf_transport + yvec_surf*df
    
  ELSE IF (ny == 1) THEN
    
    IF (ierode == 1) THEN
      surf_transport = df*bbu(jx,jy,jz)*ssurf(is,jx,jy,jz)
    ELSE
      surf_transport = 0.0
    END IF
    fxx(ind) = MultiplyCell*(surf_accum) + surf_transport + xvec_surf*df
    
  ELSE
    
    IF (ierode == 1) THEN
      surf_transport = df*ssurf(is,jx,jy,jz)*( bbu(jx,jy,jz) + ebu(jx,jy,jz) )
    ELSE
      surf_transport = 0.0
    END IF
    fxx(ind) = MultiplyCell*(surf_accum) + surf_transport + xvec_surf*df + yvec_surf*df

  END IF
  
END DO

800 IF (activecell(jx,jy,jz) == 0) THEN

  DO i = 1,ncomp
    ind = (j-1)*(neqn) + i
    fxx(ind) = sp10(i,jx,jy,jz) - spcond10(i,jinit(jx,jy,jz))
  END DO
  DO ix = 1,nexchange
    ind = (j-1)*(neqn) + ix+ncomp
    fxx(ind) = spex10(ix,jx,jy,jz) - spcondex10(ix,jinit(jx,jy,jz))      &
          *xgram(jx,jy,jz)*ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
  END DO
  DO is = 1,nsurf
    ind = (j-1)*(neqn) + is+ncomp+nexchange
    fxx(ind) = spsurf10(is,jx,jy,jz) - spcondsurf10(is,jinit(jx,jy,jz))   &
          *xgram(jx,jy,jz)*ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
  END DO
END IF

RETURN
END SUBROUTINE FxTopeGlobal
!***************************************************************************
