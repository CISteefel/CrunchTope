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
    
SUBROUTINE AllocateGIMRT(ncomp,nspec,ngas,nrct,nexchange,nsurf,nsurf_sec,npot,neqn,nx,ny,nz)
USE crunchtype
USE params
USE runtime
USE crunch_interface
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE flow
USE temperature
USE io
USE ReadFlow
USE modflowModule

IMPLICIT NONE

!!  External arrays and variables

INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nspec
INTEGER(I4B), INTENT(IN)                           :: ngas
INTEGER(I4B), INTENT(IN)                           :: nrct
INTEGER(I4B), INTENT(IN)                           :: nexchange
INTEGER(I4B), INTENT(IN)                           :: nsurf
INTEGER(I4B), INTENT(IN)                           :: nsurf_sec
INTEGER(I4B), INTENT(IN)                           :: npot
INTEGER(I4B), INTENT(IN)                           :: neqn
INTEGER(I4B), INTENT(IN)                           :: nx
INTEGER(I4B), INTENT(IN)                           :: ny
INTEGER(I4B), INTENT(IN)                           :: nz

!!  Internal arrays and variables

INTEGER(I4B)                                       :: nxyz

nxyz = nx*ny*nz

IF (ALLOCATED(sce)) THEN
  DEALLOCATE(sce)
END IF
ALLOCATE(sce(ncomp))
IF (ALLOCATED(scw)) THEN
  DEALLOCATE(scw)
END IF
ALLOCATE(scw(ncomp))
IF (ALLOCATED(scn)) THEN
  DEALLOCATE(scn)
END IF
ALLOCATE(scn(ncomp))
IF (ALLOCATED(scs)) THEN
  DEALLOCATE(scs)
END IF
ALLOCATE(scs(ncomp))


sce = 0.0d0
scw = 0.0d0
scn = 0.0d0
scs = 0.0d0

IF (ALLOCATED(gam)) THEN
  DEALLOCATE(gam)
END IF
ALLOCATE(gam(ncomp+nspec,nx,ny,nz))

IF (ALLOCATED(keqaq)) THEN
  DEALLOCATE(keqaq)
END IF
ALLOCATE(keqaq(nspec,nx,ny,nz))

IF (ALLOCATED(keqgas)) THEN
  DEALLOCATE(keqgas)
END IF
ALLOCATE(keqgas(ngas,nx,ny,nz))

IF (ALLOCATED(keqsurf)) THEN
  DEALLOCATE(keqsurf)
END IF
ALLOCATE(keqsurf(nsurf_sec,nx,ny,nz))

IF (ALLOCATED(keqmin)) THEN
  DEALLOCATE(keqmin)
END IF
ALLOCATE(keqmin(nreactmax,nrct,nx,ny,nz))
  IF (ALLOCATED(silogGlobal)) THEN
    DEALLOCATE(silogGlobal)
  END IF
ALLOCATE(silogGlobal(nreactmax,nrct,nx,ny,nz))

gam = 0.0
keqaq = 500.00
keqgas = 500.00
keqsurf = 500.00
keqmin = 500.00

IF (ALLOCATED(a)) THEN
  DEALLOCATE(a)
END IF
ALLOCATE(a(nx,ny,nz))
IF (ALLOCATED(b)) THEN
  DEALLOCATE(b)
END IF
ALLOCATE(b(nx,ny,nz))
IF (ALLOCATED(c)) THEN
  DEALLOCATE(c)
END IF
ALLOCATE(c(nx,ny,nz))
IF (ALLOCATED(d)) THEN
  DEALLOCATE(d)
END IF
ALLOCATE(d(nx,ny,nz))
IF (ALLOCATED(e)) THEN
  DEALLOCATE(e)
END IF
ALLOCATE(e(nx,ny,nz))
IF (ALLOCATED(f)) THEN
  DEALLOCATE(f)
END IF
ALLOCATE(f(nx,ny,nz))
a = 0.0
b = 0.0
c = 0.0
d = 0.0
e = 0.0
f = 0.0
IF (ALLOCATED(ssurfn)) THEN
  DEALLOCATE(ssurfn)
END IF
ALLOCATE(ssurfn(nsurf,nx,ny,nz))
IF (ALLOCATED(LogTotalSurface)) THEN
  DEALLOCATE(LogTotalSurface)
END IF
ALLOCATE(LogTotalSurface(nsurf,nx,ny,nz))
IF (ALLOCATED(sexold)) THEN
  DEALLOCATE(sexold)
END IF
ALLOCATE(sexold(ncomp,nx,ny,nz))
IF (ALLOCATED(ssurfold)) THEN
  DEALLOCATE(ssurfold)
END IF
ALLOCATE(ssurfold(ncomp,nx,ny,nz))
ssurfn = 0.0
sexold = 0.0
ssurfold = 0.0
IF (ALLOCATED(fjac)) THEN
  DEALLOCATE(fjac)
END IF
ALLOCATE(fjac(ncomp,ncomp,nx,ny,nz))
fjac = 0.0
IF (ALLOCATED(fxx)) THEN
  DEALLOCATE(fxx)
END IF
ALLOCATE(fxx(neqn*nxyz))
fxx = 0.0
IF (nx > 1 .AND. ny == 1) THEN
  IF (ALLOCATED(alf)) THEN
    DEALLOCATE(alf)
  END IF
  ALLOCATE(alf(neqn,neqn,3))
ELSE 
  IF (ALLOCATED(alf)) THEN
    DEALLOCATE(alf)
  END IF
  ALLOCATE(alf(neqn,neqn,5))
END IF
alf = 0.0
IF (ALLOCATED(blockfl)) THEN
  DEALLOCATE(blockfl)
END IF
ALLOCATE(blockfl(neqn,neqn))
IF (ALLOCATED(blockl)) THEN
  DEALLOCATE(blockl)
END IF
ALLOCATE(blockl(neqn,neqn))
IF (ALLOCATED(blockm)) THEN
  DEALLOCATE(blockm)
END IF
ALLOCATE(blockm(neqn,neqn))
IF (ALLOCATED(blockr)) THEN
  DEALLOCATE(blockr)
END IF
ALLOCATE(blockr(neqn,neqn))
IF (ALLOCATED(blockfr)) THEN
  DEALLOCATE(blockfr)
END IF
ALLOCATE(blockfr(neqn,neqn))

IF (giambalvo) THEN
  IF (ALLOCATED(dflux_x)) THEN
    DEALLOCATE(dflux_x)
  END IF
  ALLOCATE(dflux_x(ncomp,ny,2))
  IF (ALLOCATED(advflux_x)) THEN
    DEALLOCATE(advflux_x)
  END IF
  ALLOCATE(advflux_x(ncomp,ny,2))
END IF

IF (ierode == 1) THEN
  IF (ALLOCATED(aabur)) THEN
    DEALLOCATE(aabur)
  END IF
  ALLOCATE(aabur(nx))
  IF (ALLOCATED(bbbur)) THEN
    DEALLOCATE(bbbur)
  END IF
  ALLOCATE(bbbur(nx))
  IF (ALLOCATED(ccbur)) THEN
    DEALLOCATE(ccbur)
  END IF
  ALLOCATE(ccbur(nx))
  IF (ALLOCATED(rrbur)) THEN
    DEALLOCATE(rrbur)
  END IF
  ALLOCATE(rrbur(nx))
  IF (ALLOCATED(uubur)) THEN
    DEALLOCATE(uubur)
  END IF
  ALLOCATE(uubur(nx))
  aabur = 0.0
  bbbur = 0.0
  ccbur = 0.0
  rrbur = 0.0
  uubur = 0.0
  IF (ALLOCATED(abu)) THEN
    DEALLOCATE(abu)
  END IF
  ALLOCATE(abu(nx,ny,nz))
  IF (ALLOCATED(bbu)) THEN
    DEALLOCATE(bbu)
  END IF
  ALLOCATE(bbu(nx,ny,nz))
  IF (ALLOCATED(cbu)) THEN
    DEALLOCATE(cbu)
  END IF
  ALLOCATE(cbu(nx,ny,nz))
  IF (ALLOCATED(dbu)) THEN
    DEALLOCATE(dbu)
  END IF
  ALLOCATE(dbu(nx,ny,nz))
  IF (ALLOCATED(ebu)) THEN
    DEALLOCATE(ebu)
  END IF
  ALLOCATE(ebu(nx,ny,nz))
  IF (ALLOCATED(fbu)) THEN
    DEALLOCATE(fbu)
  END IF
  ALLOCATE(fbu(nx,ny,nz))
  abu = 0.0
  bbu = 0.0
  cbu = 0.0
  dbu = 0.0
  ebu = 0.0
  fbu = 0.0
  IF (ALLOCATED(ctvd)) THEN
    DEALLOCATE(ctvd)
  END IF
  ALLOCATE(ctvd(-1:nx+2,-1:ny+2,1))
  ctvd = 0.0
  IF (ALLOCATED(sch)) THEN
    DEALLOCATE(sch)
  END IF
  ALLOCATE(sch(ncomp,nx,ny,nz))
  
  IF (ALLOCATED(ssurf)) THEN
    DEALLOCATE(ssurf)
  END IF
  ALLOCATE(ssurf(nsurf,nx,ny,nz))
  
  IF (ALLOCATED(fch)) THEN
    DEALLOCATE(fch)
  END IF
  ALLOCATE(fch(ncomp,neqn,nx,ny,nz))
  IF (ALLOCATED(fsurf)) THEN
    DEALLOCATE(fsurf)
  END IF
  ALLOCATE(fsurf(nsurf,nsurf+ncomp,nx,ny,nz))
  sch = 0.0
  ssurf = 0.0
  fch = 0.0
  fsurf = 0.0
  IF (ALLOCATED(fjpotncomp)) THEN
    DEALLOCATE(fjpotncomp)
  END IF
  ALLOCATE(fjpotncomp(npot,ncomp,nx,ny,nz))
  IF (ALLOCATED(fjpotnsurf)) THEN
    DEALLOCATE(fjpotnsurf)
  END IF
  ALLOCATE(fjpotnsurf(npot,nsurf,nx,ny,nz))



  IF (ALLOCATED(surfcharge)) THEN
    DEALLOCATE(surfcharge)
  END IF
  ALLOCATE(surfcharge(nrct))
  fjpotncomp = 0.0
  fjpotnsurf = 0.0
  LogPotential = 0.0
  surfcharge = 0.0
ELSE IF (ierode /= 1) THEN
  IF (ALLOCATED(sNCexch_local)) THEN
    DEALLOCATE(sNCexch_local)
  END IF
  ALLOCATE(sNCexch_local(ncomp))
  IF (ALLOCATED(sNCsurf_local)) THEN
    DEALLOCATE(sNCsurf_local)
  END IF
  ALLOCATE(sNCsurf_local(ncomp))
  IF (ALLOCATED(ssurf_local)) THEN
    DEALLOCATE(ssurf_local)
  END IF
  ALLOCATE(ssurf_local(nsurf))
  IF (ALLOCATED(fch_local)) THEN
    DEALLOCATE(fch_local)
  END IF
  ALLOCATE(fch_local(ncomp,neqn))
  IF (ALLOCATED(fsurf_local)) THEN
    DEALLOCATE(fsurf_local)
  END IF
  ALLOCATE(fsurf_local(nsurf,nsurf+ncomp))
  ssurf_local = 0.0
  fch_local = 0.0
  fsurf_local = 0.0
  sNCexch_local = 0.0
  sNCsurf_local = 0.0
  IF (ALLOCATED(fjpotncomp)) THEN
    DEALLOCATE(fjpotncomp)
  END IF
  ALLOCATE(fjpotncomp(npot,ncomp,nx,ny,nz))
  IF (ALLOCATED(fjpotnsurf)) THEN
    DEALLOCATE(fjpotnsurf)
  END IF
  ALLOCATE(fjpotnsurf(npot,nsurf,nx,ny,nz))
  IF (ALLOCATED(LogPotential)) THEN
    DEALLOCATE(LogPotential)
  END IF
  ALLOCATE(LogPotential(nsurf,nx,ny,nz))
  IF (ALLOCATED(surfcharge)) THEN
    DEALLOCATE(surfcharge)
  END IF
  ALLOCATE(surfcharge(nrct))
  fjpotncomp = 0.0
  fjpotnsurf = 0.0
  LogPotential = 0.0
  surfcharge = 0.0
ELSE
  CONTINUE
END IF

IF (isaturate == 1) THEN
  IF (ALLOCATED(sge)) THEN
    DEALLOCATE(sge)
  END IF
  ALLOCATE(sge(ncomp))
  IF (ALLOCATED(sgw)) THEN
    DEALLOCATE(sgw)
  END IF
  ALLOCATE(sgw(ncomp))
  IF (ALLOCATED(sgn)) THEN
    DEALLOCATE(sgn)
  END IF
  ALLOCATE(sgn(ncomp))
  IF (ALLOCATED(sgs)) THEN
    DEALLOCATE(sgs)
  END IF
  ALLOCATE(sgs(ncomp))
  IF (ALLOCATED(sgaspump)) THEN
    DEALLOCATE(sgaspump)
  END IF
  ALLOCATE(sgaspump(ncomp))
  IF (ALLOCATED(ag)) THEN
    DEALLOCATE(ag)
  END IF
  ALLOCATE(ag(nx,ny,nz))
  IF (ALLOCATED(bg)) THEN
    DEALLOCATE(bg)
  END IF
  ALLOCATE(bg(nx,ny,nz))
  IF (ALLOCATED(cg)) THEN
    DEALLOCATE(cg)
  END IF
  ALLOCATE(cg(nx,ny,nz))
  IF (ALLOCATED(dg)) THEN
    DEALLOCATE(dg)
  END IF
  ALLOCATE(dg(nx,ny,nz))
  IF (ALLOCATED(eg)) THEN
    DEALLOCATE(eg)
  END IF
  ALLOCATE(eg(nx,ny,nz))
  IF (ALLOCATED(fg)) THEN
    DEALLOCATE(fg)
  END IF
ALLOCATE(fg(0:nx+1,0:ny+1,nz))

  ag = 0.0
  bg = 0.0
  cg = 0.0
  dg = 0.0
  eg = 0.0
  fg = 0.0
  sge = 0.0d0
  sgw = 0.0d0
  sgn = 0.0d0
  sgs = 0.0d0
  sgaspump = 0.0d0
  IF (ALLOCATED(sgas)) THEN
    DEALLOCATE(sgas)
  END IF
  ALLOCATE(sgas(ncomp,0:nx+1,0:ny+1,0:nz+1))
  IF (ALLOCATED(sgasn)) THEN
    DEALLOCATE(sgasn)
  END IF
  ALLOCATE(sgasn(ncomp,nx,ny,nz))
  IF (ALLOCATED(fgas)) THEN
    DEALLOCATE(fgas)
  END IF
  ALLOCATE(fgas(ncomp,ncomp,nx,ny,nz))
  sgas = 0.0
  sgasn = 0.0
  fgas = 0.0
END IF

IF (species_diffusion) THEN
  IF (ALLOCATED(schg)) THEN
    DEALLOCATE(schg)
  END IF
  ALLOCATE(schg(ncomp))
  IF (ALLOCATED(sdsp)) THEN
    DEALLOCATE(sdsp)
  END IF
  ALLOCATE(sdsp(ncomp))
  IF (ALLOCATED(fjac_d)) THEN
    DEALLOCATE(fjac_d)
  END IF
  ALLOCATE(fjac_d(ncomp,ncomp,nx,ny,nz))
  IF (ALLOCATED(fjac_chg)) THEN
    DEALLOCATE(fjac_chg)
  END IF
  ALLOCATE(fjac_chg(ncomp,ncomp,nx,ny,nz))
  IF (ALLOCATED(a_d)) THEN
    DEALLOCATE(a_d)
  END IF
  ALLOCATE(a_d(nx,ny,nz))
  IF (ALLOCATED(b_d)) THEN
    DEALLOCATE(b_d)
  END IF
  ALLOCATE(b_d(nx,ny,nz))
  IF (ALLOCATED(c_d)) THEN
    DEALLOCATE(c_d)
  END IF
  ALLOCATE(c_d(nx,ny,nz))
  IF (ALLOCATED(d_d)) THEN
    DEALLOCATE(d_d)
  END IF
  ALLOCATE(d_d(nx,ny,nz))
  IF (ALLOCATED(f_d)) THEN
    DEALLOCATE(f_d)
  END IF
  ALLOCATE(f_d(nx,ny,nz))
  IF (ALLOCATED(e_d)) THEN
    DEALLOCATE(e_d)
  END IF
  ALLOCATE(e_d(nx,ny,nz))
  IF (ALLOCATED(s_dsp)) THEN
    DEALLOCATE(s_dsp)
  END IF
  ALLOCATE(s_dsp(ncomp,0:nx+1,0:ny+1,nz))
  IF (ALLOCATED(s_chg)) THEN
    DEALLOCATE(s_chg)
  END IF
  ALLOCATE(s_chg(ncomp,0:nx+1,0:ny+1,nz))
  IF (ALLOCATED(d_correct)) THEN
    DEALLOCATE(d_correct)
  END IF
  ALLOCATE(d_correct(ncomp+nspec))
  IF (ALLOCATED(sumwtchg)) THEN
    DEALLOCATE(sumwtchg)
  END IF
  ALLOCATE(sumwtchg(0:nx+1,0:ny+1,nz))
  IF (ALLOCATED(sigma_w)) THEN
    DEALLOCATE(sigma_w)
  END IF
  ALLOCATE(sigma_w(ncomp))
  IF (ALLOCATED(sigma_e)) THEN
    DEALLOCATE(sigma_e)
  END IF
  ALLOCATE(sigma_e(ncomp))
  IF (ALLOCATED(sigma_s)) THEN
    DEALLOCATE(sigma_s)
  END IF
  ALLOCATE(sigma_s(ncomp))
  IF (ALLOCATED(sigma_n)) THEN
    DEALLOCATE(sigma_n)
  END IF
  ALLOCATE(sigma_n(ncomp))
  IF (ALLOCATED(fanalyt)) THEN
    DEALLOCATE(fanalyt)
  END IF
  ALLOCATE(fanalyt(ncomp,ncomp))
  schg = 0.0d0
  sdsp = 0.0d0
  fanalyt = 0.0d0
  fjac_d = 0.0
  fjac_chg = 0.0
  a_d = 0.0
  b_d = 0.0
  c_d = 0.0
  d_d = 0.0
  f_d = 0.0
  e_d = 0.0
  s_dsp = 0.0
  s_chg = 0.0
  d_correct = 0.0
  sumwtchg = 0.0
  sigma_w = 0.0
  sigma_e = 0.0
  sigma_s = 0.0
  sigma_n = 0.0
END IF

IF (nxyz == 1) THEN
  IF (ALLOCATED(aaa)) THEN
    DEALLOCATE(aaa)
  END IF
  ALLOCATE(aaa(neqn,neqn))
  IF (ALLOCATED(bb)) THEN
    DEALLOCATE(bb)
  END IF
  ALLOCATE(bb(neqn))
  IF (ALLOCATED(xn)) THEN
    DEALLOCATE(xn)
  END IF
  ALLOCATE(xn(neqn))
  IF (ALLOCATED(indd)) THEN
    DEALLOCATE(indd)
  END IF
  ALLOCATE(indd(neqn))
  aaa = 0.0
  bb = 0.0
  xn = 0.0
  indd = 0
ELSE IF (nxyz == nx .AND. ihindmarsh == 1) THEN
  IF (ALLOCATED(aah)) THEN
    DEALLOCATE(aah)
  END IF
  ALLOCATE(aah(neqn,neqn,nx))
  IF (ALLOCATED(bbh)) THEN
    DEALLOCATE(bbh)
  END IF
  ALLOCATE(bbh(neqn,neqn,nx))
  IF (ALLOCATED(cch)) THEN
    DEALLOCATE(cch)
  END IF
  ALLOCATE(cch(neqn,neqn,nx))
  IF (ALLOCATED(yh)) THEN
    DEALLOCATE(yh)
  END IF
  ALLOCATE(yh(neqn,nx))
  IF (ALLOCATED(indexx)) THEN
    DEALLOCATE(indexx)
  END IF
  ALLOCATE(indexx(neqn,nx))
  aah = 0.0
  bbh = 0.0
  cch = 0.0
  yh = 0.0
  indexx = 0
  IF (ALLOCATED(aaa)) THEN
    DEALLOCATE(aaa)
  END IF
  ALLOCATE(aaa(neqn,neqn))
  IF (ALLOCATED(bb)) THEN
    DEALLOCATE(bb)
  END IF
  ALLOCATE(bb(neqn))
  IF (ALLOCATED(indd)) THEN
    DEALLOCATE(indd)
  END IF
  ALLOCATE(indd(neqn))
  aaa = 0.0
  bb = 0.0
  indd = 0
ELSE 
  IF (ALLOCATED(xn)) THEN
    DEALLOCATE(xn)
  END IF
  ALLOCATE(xn(neqn*nx*ny*nz))
  IF (ALLOCATED(aah)) THEN
    DEALLOCATE(aah)
  END IF
  ALLOCATE(aah(neqn,neqn,nx))
  IF (ALLOCATED(bbh)) THEN
    DEALLOCATE(bbh)
  END IF
  ALLOCATE(bbh(neqn,neqn,nx))
  IF (ALLOCATED(cch)) THEN
    DEALLOCATE(cch)
  END IF
  ALLOCATE(cch(neqn,neqn,nx))
  IF (ALLOCATED(yh)) THEN
    DEALLOCATE(yh)
  END IF
  ALLOCATE(yh(neqn,nx))
  IF (ALLOCATED(indexx)) THEN
    DEALLOCATE(indexx)
  END IF
  ALLOCATE(indexx(neqn,nx))
  aah = 0.0
  bbh = 0.0
  cch = 0.0
  yh = 0.0
  indexx = 0
END IF

!!IF (nsurf > 0) THEN
  IF (ALLOCATED(SurfaceCon)) THEN
    DEALLOCATE(SurfaceCon)
  END IF
  ALLOCATE(SurfaceCon(ncomp))
  IF (ALLOCATED(surf_east)) THEN
    DEALLOCATE(surf_east)
  END IF
  ALLOCATE(surf_east(ncomp))
  IF (ALLOCATED(surf_west)) THEN
    DEALLOCATE(surf_west)
  END IF
  ALLOCATE(surf_west(ncomp))
  IF (ALLOCATED(surf_north)) THEN
    DEALLOCATE(surf_north)
  END IF
  ALLOCATE(surf_north(ncomp))
  IF (ALLOCATED(surf_south)) THEN
    DEALLOCATE(surf_south)
  END IF
  ALLOCATE(surf_south(ncomp))
  surf_east = 0.0d0
  surf_west = 0.0d0
  surf_north = 0.0d0
  surf_south = 0.0d0
!!END IF
!!IF (nexchange > 0) THEN
  IF (ALLOCATED(ExchangeCon)) THEN
    DEALLOCATE(ExchangeCon)
  END IF
  ALLOCATE(ExchangeCon(ncomp))
  IF (ALLOCATED(fweight)) THEN
    DEALLOCATE(fweight)
  END IF
  ALLOCATE(fweight(neqn,neqn))
  IF (ALLOCATED(sex_east)) THEN
    DEALLOCATE(sex_east)
  END IF
  ALLOCATE(sex_east(ncomp))
  IF (ALLOCATED(sex_west)) THEN
    DEALLOCATE(sex_west)
  END IF
  ALLOCATE(sex_west(ncomp))
  IF (ALLOCATED(sex_north)) THEN
    DEALLOCATE(sex_north)
  END IF
  ALLOCATE(sex_north(ncomp))
  IF (ALLOCATED(sex_south)) THEN
    DEALLOCATE(sex_south)
  END IF
  ALLOCATE(sex_south(ncomp))
  sex_east = 0.0d0
  sex_west = 0.0d0
  sex_north = 0.0d0
  sex_south = 0.0d0
!!END IF

IF (ALLOCATED(IntegerArray)) THEN
  DEALLOCATE(IntegerArray)
END IF
ALLOCATE(IntegerArray(nx,ny,nz))
IF (ALLOCATED(RealArray)) THEN
  DEALLOCATE(RealArray)
END IF
ALLOCATE(RealArray(nx,ny,nz))

IF (Duan .OR. Duan2006) THEN
  IF (ALLOCATED(GasPressureTotal)) THEN
    DEALLOCATE(GasPressureTotal)
  END IF
  ALLOCATE(GasPressureTotal(nx,ny,nz))
END IF


RETURN
END SUBROUTINE AllocateGIMRT