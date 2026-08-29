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
    
SUBROUTINE AllocateOS3D(ncomp,nspec,ngas,nrct,nexchange,nsurf,nsurf_sec,npot,neqn,nx,ny,nz)
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

  IF (ALLOCATED(cflowin)) THEN
    DEALLOCATE(cflowin)
  END IF
  ALLOCATE(cflowin(ncomp))
  IF (ALLOCATED(cflowout)) THEN
    DEALLOCATE(cflowout)
  END IF
  ALLOCATE(cflowout(ncomp))
  cflowin = 0.0
  cflowout = 0.0
  IF (ALLOCATED(ftvd)) THEN
    DEALLOCATE(ftvd)
  END IF
  ALLOCATE(ftvd(0:nx,0:ny,0:nz))
  IF (ALLOCATED(gtvd)) THEN
    DEALLOCATE(gtvd)
  END IF
  ALLOCATE(gtvd(0:nx,0:ny,0:nz))
  IF (ALLOCATED(htvd)) THEN
    DEALLOCATE(htvd)
  END IF
  ALLOCATE(htvd(0:nx,0:ny,0:nz))
  IF (ALLOCATED(sorp)) THEN
    DEALLOCATE(sorp)
  END IF
  ALLOCATE(sorp(nx,ny,nz))
  sorp = 1.0
  
  IF (ALLOCATED(gamma)) THEN
    DEALLOCATE(gamma)
  END IF
  ALLOCATE(gamma(ncomp+nspec,nx,ny,nz))
  
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
  
  IF (ALLOCATED(fjac_loc)) THEN
    DEALLOCATE(fjac_loc)
  END IF
  ALLOCATE(fjac_loc(ncomp,ncomp))
  
  gamma = 0.0
  keqaq = 500.00
  keqgas = 500.00
  keqsurf = 500.00
  keqmin = 500.00
  fjac_loc = 0.0
  
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
  IF (ALLOCATED(skdold)) THEN
    DEALLOCATE(skdold)
  END IF
  ALLOCATE(skdold(ncomp,nx,ny,nz))
  ssurfn = 0.0
  sexold = 0.0
  ssurfold = 0.0
  skdold = 0.0
  IF (ALLOCATED(fxx)) THEN
    DEALLOCATE(fxx)
  END IF
  ALLOCATE(fxx(neqn))
  fxx = 0.0
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
  IF (ALLOCATED(ctvd)) THEN
    DEALLOCATE(ctvd)
  END IF
  ALLOCATE(ctvd(-1:nx+2,-1:ny+2,-1:nz+2))
  ctvd = 0.0
  IF (ALLOCATED(LogPotential)) THEN
    DEALLOCATE(LogPotential)
  END IF
  ALLOCATE(LogPotential(nsurf,nx,ny,nz))
  LogPotential = 0.0
  IF (ALLOCATED(surfcharge)) THEN
    DEALLOCATE(surfcharge)
  END IF
  ALLOCATE(surfcharge(nrct))
  surfcharge = 0.0
  IF (ALLOCATED(fjpotncomp_local)) THEN
    DEALLOCATE(fjpotncomp_local)
  END IF
  ALLOCATE(fjpotncomp_local(npot,ncomp))
  IF (ALLOCATED(fjpotnsurf_local)) THEN
    DEALLOCATE(fjpotnsurf_local)
  END IF
  ALLOCATE(fjpotnsurf_local(npot,nsurf))
  fjpotncomp_local = 0.0
  fjpotnsurf_local = 0.0
  IF (spherical) THEN
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
  END IF
  
IF (Duan .OR. Duan2006) THEN
  IF (ALLOCATED(GasPressureTotal)) THEN
    DEALLOCATE(GasPressureTotal)
  END IF
  ALLOCATE(GasPressureTotal(nx,ny,nz))
END IF

!!IF (nsurf > 0) THEN
  IF (ALLOCATED(SurfaceCon)) THEN
    DEALLOCATE(SurfaceCon)
  END IF
  ALLOCATE(SurfaceCon(ncomp))
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
!!END IF

IF (ALLOCATED(IntegerArray)) THEN
  DEALLOCATE(IntegerArray)
END IF
ALLOCATE(IntegerArray(nx,ny,nz))
IF (ALLOCATED(RealArray)) THEN
  DEALLOCATE(RealArray)
END IF
ALLOCATE(RealArray(nx,ny,nz))

RETURN
END SUBROUTINE AllocateOS3D