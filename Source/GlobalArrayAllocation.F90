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

SUBROUTINE GlobalArrayAllocation(ncomp,nspec,nkin,nrct,ngas,npot,nexchange,nexch_sec,nsurf,nsurf_sec,ikin,nx,ny,nz)
USE crunchtype
USE params
USE runtime
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE flow
USE temperature
USE io
USE strings
USE ReadFlow
USE modflowModule

IMPLICIT NONE

INTEGER(I4B), INTENT(IN)   :: ncomp,nspec,nkin,nrct,ngas,npot,nexchange,nexch_sec,nsurf,nsurf_sec,ikin,nx,ny,nz

IF (ALLOCATED(qx)) THEN
  DEALLOCATE(qx)
  ALLOCATE(qx(0:nx,ny,nz))
ELSE
  ALLOCATE(qx(0:nx,ny,nz))
END IF
IF (ALLOCATED(qy)) THEN
  DEALLOCATE(qy)
  ALLOCATE(qy(nx,0:ny,nz))
ELSE
  ALLOCATE(qy(nx,0:ny,nz))
END IF
IF (ALLOCATED(qz)) THEN
  DEALLOCATE(qz)
  ALLOCATE(qz(nx,ny,0:nz))
ELSE
  ALLOCATE(qz(nx,ny,0:nz))
END IF
IF (ALLOCATED(dspx)) THEN
  DEALLOCATE(dspx)
  ALLOCATE(dspx(nx,ny,nz))
ELSE
  ALLOCATE(dspx(nx,ny,nz))
END IF
IF (ALLOCATED(dspy)) THEN
  DEALLOCATE(dspy)
  ALLOCATE(dspy(nx,ny,nz))
ELSE
  ALLOCATE(dspy(nx,ny,nz))
END IF
IF (ALLOCATED(dspz)) THEN
  DEALLOCATE(dspz)
  ALLOCATE(dspz(nx,ny,nz))
ELSE
  ALLOCATE(dspz(nx,ny,nz))
END IF
IF (ALLOCATED(qxgas)) THEN
  DEALLOCATE(qxgas)
  ALLOCATE(qxgas(0:nx,ny,nz))
ELSE
  ALLOCATE(qxgas(0:nx,ny,nz))
END IF
IF (ALLOCATED(qygas)) THEN
  DEALLOCATE(qygas)
  ALLOCATE(qygas(nx,0:ny,nz))
ELSE
  ALLOCATE(qygas(nx,0:ny,nz))
END IF
IF (ALLOCATED(qzgas)) THEN
  DEALLOCATE(qzgas)
  ALLOCATE(qzgas(nx,ny,0:nz))
ELSE
  ALLOCATE(qzgas(nx,ny,0:nz))
END IF
IF (ALLOCATED(netflowx)) THEN
  DEALLOCATE(netflowx)
  ALLOCATE(netflowx(0:nx,ny,nz))
ELSE
  ALLOCATE(netflowx(0:nx,ny,nz))
END IF
IF (ALLOCATED(netDiffuseX)) THEN
  DEALLOCATE(netDiffuseX)
  ALLOCATE(netDiffuseX(0:nx,ny,nz))
ELSE
  ALLOCATE(netDiffuseX(0:nx,ny,nz))
END IF
IF (ALLOCATED(satliq)) THEN
  DEALLOCATE(satliq)
  ALLOCATE(satliq(-1:nx+2,-1:ny+2,-1:nz+2))
ELSE
  ALLOCATE(satliq(-1:nx+2,-1:ny+2,-1:nz+2))
END IF
IF (ALLOCATED(satliqold)) THEN
  DEALLOCATE(satliqold)
  ALLOCATE(satliqold(-1:nx+2,-1:ny+2,-1:nz+2))
ELSE
  ALLOCATE(satliqold(-1:nx+2,-1:ny+2,-1:nz+2))
END IF
IF (ALLOCATED(dstar)) THEN
  DEALLOCATE(dstar)
  ALLOCATE(dstar(-1:nx+2,-1:ny+2,nz))
ELSE
  ALLOCATE(dstar(-1:nx+2,-1:ny+2,nz))
END IF
IF (ALLOCATED(por)) THEN
  DEALLOCATE(por)
  ALLOCATE(por(-1:nx+2,-1:ny+2,-1:nz+2))
ELSE
  ALLOCATE(por(-1:nx+2,-1:ny+2,-1:nz+2))
END IF
IF (ALLOCATED(porOld)) THEN
  DEALLOCATE(porOld)
  ALLOCATE(porOld(-1:nx+2,-1:ny+2,-1:nz+2))
ELSE
  ALLOCATE(porOld(-1:nx+2,-1:ny+2,-1:nz+2))
END IF
IF (ALLOCATED(porin)) THEN
  DEALLOCATE(porin)
  ALLOCATE(porin(-1:nx+2,-1:ny+2,-1:nz+2))
ELSE
  ALLOCATE(porin(-1:nx+2,-1:ny+2,-1:nz+2))
END IF

IF (ALLOCATED(pres)) THEN
  DEALLOCATE(pres)
  ALLOCATE(pres(0:nx+1,0:ny+1,0:nz+1))
ELSE
  ALLOCATE(pres(0:nx+1,0:ny+1,0:nz+1))
END IF
IF (ALLOCATED(activecell)) THEN
  DEALLOCATE(activecell)
  ALLOCATE(activecell(nx,ny,nz))
ELSE
  ALLOCATE(activecell(nx,ny,nz))
END IF
IF (ALLOCATED(dxy)) THEN
  DEALLOCATE(dxy)
  ALLOCATE(dxy(nx,ny,nz))
ELSE
  ALLOCATE(dxy(nx,ny,nz))
END IF
IF (ALLOCATED(spno2)) THEN
  DEALLOCATE(spno2)
  ALLOCATE(spno2(nx,ny,nz))
ELSE
  ALLOCATE(spno2(nx,ny,nz))
END IF
IF (ALLOCATED(spnno2)) THEN
  DEALLOCATE(spnno2)
  ALLOCATE(spnno2(nx,ny,nz))
ELSE
  ALLOCATE(spnno2(nx,ny,nz))
END IF


IF (ny == 1 .AND. nz == 1) THEN

  IF (ALLOCATED(s)) THEN
    DEALLOCATE(s)
    ALLOCATE(s(ncomp,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(s(ncomp,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(sp)) THEN
    DEALLOCATE(sp)
    ALLOCATE(sp(ncomp+nspec,0:nx+1,ny,nz))
  ELSE
    ALLOCATE(sp(ncomp+nspec,0:nx+1,ny,nz))
  END IF
  IF (ALLOCATED(sp10)) THEN
    DEALLOCATE(sp10)
    ALLOCATE(sp10(ncomp+nspec,0:nx+1,ny,nz))
  ELSE
    ALLOCATE(sp10(ncomp+nspec,0:nx+1,ny,nz))
  END IF
  IF (ALLOCATED(spex)) THEN
    DEALLOCATE(spex)
    ALLOCATE(spex(nexchange+nexch_sec,0:nx+1,ny,nz))
  ELSE
    ALLOCATE(spex(nexchange+nexch_sec,0:nx+1,ny,nz))
  END IF
  IF (ALLOCATED(spsurf)) THEN
    DEALLOCATE(spsurf)
    ALLOCATE(spsurf(nsurf+nsurf_sec,0:nx+1,ny,nz))
  ELSE
    ALLOCATE(spsurf(nsurf+nsurf_sec,0:nx+1,ny,nz))
  END IF
  IF (ALLOCATED(spgas)) THEN
    DEALLOCATE(spgas)
    ALLOCATE(spgas(ngas,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(spgas(ngas,0:nx+1,ny,nz))
  END IF
  IF (ALLOCATED(LogPotential)) THEN
    DEALLOCATE(LogPotential)
    ALLOCATE(LogPotential(nsurf,0:nx+1,ny,nz))
  ELSE
    ALLOCATE(LogPotential(nsurf,0:nx+1,ny,nz))
  END IF
  IF (ALLOCATED(spex10)) THEN
    DEALLOCATE(spex10)
    ALLOCATE(spex10(nexchange+nexch_sec,0:nx+1,ny,nz))
  ELSE
    ALLOCATE(spex10(nexchange+nexch_sec,0:nx+1,ny,nz))
  END IF
  IF (ALLOCATED(spsurf10)) THEN
    DEALLOCATE(spsurf10)
    ALLOCATE(spsurf10(nsurf+nsurf_sec,0:nx+1,ny,nz))
  ELSE
    ALLOCATE(spsurf10(nsurf+nsurf_sec,0:nx+1,ny,nz))
  END IF
  IF (ALLOCATED(spgas10)) THEN
    DEALLOCATE(spgas10)
    ALLOCATE(spgas10(ngas,0:nx+1,ny,nz))
  ELSE
    ALLOCATE(spgas10(ngas,0:nx+1,ny,nz))
  END IF

  IF (ALLOCATED(t)) THEN
    DEALLOCATE(t)
    ALLOCATE(t(0:nx+1,ny,nz))
  ELSE
    ALLOCATE(t(0:nx+1,ny,nz))
  END IF
!!!  IF (ALLOCATED(por)) THEN
!!!    DEALLOCATE(por)
!!!    ALLOCATE(por(0:nx+1,ny,nz))
!!!  ELSE
!!!    ALLOCATE(por(0:nx+1,ny,nz))
!!!  END IF

ELSE IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN

  IF (ALLOCATED(s)) THEN
    DEALLOCATE(s)
    ALLOCATE(s(ncomp,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(s(ncomp,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(sp)) THEN
    DEALLOCATE(sp)
    ALLOCATE(sp(ncomp+nspec,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(sp(ncomp+nspec,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(sp10)) THEN
    DEALLOCATE(sp10)
    ALLOCATE(sp10(ncomp+nspec,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(sp10(ncomp+nspec,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(spex)) THEN
    DEALLOCATE(spex)
    ALLOCATE(spex(nexchange+nexch_sec,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(spex(nexchange+nexch_sec,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(spsurf)) THEN
    DEALLOCATE(spsurf)
    ALLOCATE(spsurf(nsurf+nsurf_sec,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(spsurf(nsurf+nsurf_sec,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(spgas)) THEN
    DEALLOCATE(spgas)
    ALLOCATE(spgas(ngas,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(spgas(ngas,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(LogPotential)) THEN
    DEALLOCATE(LogPotential)
    ALLOCATE(LogPotential(nsurf,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(LogPotential(nsurf,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(spex10)) THEN
    DEALLOCATE(spex10)
    ALLOCATE(spex10(nexchange+nexch_sec,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(spex10(nexchange+nexch_sec,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(spsurf10)) THEN
    DEALLOCATE(spsurf10)
    ALLOCATE(spsurf10(nsurf+nsurf_sec,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(spsurf10(nsurf+nsurf_sec,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(spgas10)) THEN
    DEALLOCATE(spgas10)
    ALLOCATE(spgas10(ngas,0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(spgas10(ngas,0:nx+1,0:ny+1,nz))
  END IF
  IF (ALLOCATED(t)) THEN
    DEALLOCATE(t)
    ALLOCATE(t(0:nx+1,0:ny+1,nz))
  ELSE
    ALLOCATE(t(0:nx+1,0:ny+1,nz))
  END IF
!!!  IF (ALLOCATED(por)) THEN
!!!    DEALLOCATE(por)
!!!    ALLOCATE(por(0:nx+1,0:ny+1,nz))
!!!  ELSE
!!!    ALLOCATE(por(0:nx+1,0:ny+1,nz))
!!!  END IF

!!!ELSE IF (nx > 1 .AND. ny == 1 .AND. nz > 1) THEN

!!!  IF (ALLOCATED(s)) THEN
!!!    DEALLOCATE(s)
!!!    ALLOCATE(s(ncomp,0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(s(ncomp,0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(sp)) THEN
!!!    DEALLOCATE(sp)
!!!    ALLOCATE(sp(ncomp+nspec,0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(sp(ncomp+nspec,0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(sp10)) THEN
!!!    DEALLOCATE(sp10)
!!!    ALLOCATE(sp10(ncomp+nspec,0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(sp10(ncomp+nspec,0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(spex)) THEN
!!!    DEALLOCATE(spex)
!!!    ALLOCATE(spex(nexchange+nexch_sec,0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(spex(nexchange+nexch_sec,0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(spsurf)) THEN
!!!    DEALLOCATE(spsurf)
!!!    ALLOCATE(spsurf(nsurf+nsurf_sec,0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(spsurf(nsurf+nsurf_sec,0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(spgas)) THEN
!!!    DEALLOCATE(spgas)
!!!    ALLOCATE(spgas(ngas,0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(spgas(ngas,0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(LogPotential)) THEN
!!!    DEALLOCATE(LogPotential)
!!!    ALLOCATE(LogPotential(nsurf,0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(LogPotential(nsurf,0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(spex10)) THEN
!!!    DEALLOCATE(spex10)
!!!    ALLOCATE(spex10(nexchange+nexch_sec,0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(spex10(nexchange+nexch_sec,0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(spsurf10)) THEN
!!!    DEALLOCATE(spsurf10)
!!!    ALLOCATE(spsurf10(nsurf+nsurf_sec,0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(spsurf10(nsurf+nsurf_sec,0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(spgas10)) THEN
!!!    DEALLOCATE(spgas10)
!!!    ALLOCATE(spgas10(ngas,0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(spgas10(ngas,0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(t)) THEN
!!!    DEALLOCATE(t)
!!!    ALLOCATE(t(0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(t(0:nx+1,ny,0:nz+1))
!!!  END IF
!!!  IF (ALLOCATED(por)) THEN
!!!    DEALLOCATE(por)
!!!    ALLOCATE(por(0:nx+1,ny,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(por(0:nx+1,ny,0:nz+1))
!!!  END IF

! Zhi Li 20200708
! ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
ELSE

  IF (ALLOCATED(s)) THEN
    DEALLOCATE(s)
    ALLOCATE(s(ncomp,0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(s(ncomp,0:nx+1,0:ny+1,0:nz+1))
  END IF
  IF (ALLOCATED(sp)) THEN
    DEALLOCATE(sp)
    ALLOCATE(sp(ncomp+nspec,0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(sp(ncomp+nspec,0:nx+1,0:ny+1,0:nz+1))
  END IF
  IF (ALLOCATED(sp10)) THEN
    DEALLOCATE(sp10)
    ALLOCATE(sp10(ncomp+nspec,0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(sp10(ncomp+nspec,0:nx+1,0:ny+1,0:nz+1))
  END IF
  IF (ALLOCATED(spex)) THEN
    DEALLOCATE(spex)
    ALLOCATE(spex(nexchange+nexch_sec,0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(spex(nexchange+nexch_sec,0:nx+1,0:ny+1,0:nz+1))
  END IF
  IF (ALLOCATED(spsurf)) THEN
    DEALLOCATE(spsurf)
    ALLOCATE(spsurf(nsurf+nsurf_sec,0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(spsurf(nsurf+nsurf_sec,0:nx+1,0:ny+1,0:nz+1))
  END IF
  IF (ALLOCATED(spgas)) THEN
    DEALLOCATE(spgas)
    ALLOCATE(spgas(ngas,0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(spgas(ngas,0:nx+1,0:ny+1,0:nz+1))
  END IF
  IF (ALLOCATED(LogPotential)) THEN
    DEALLOCATE(LogPotential)
    ALLOCATE(LogPotential(nsurf,0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(LogPotential(nsurf,0:nx+1,0:ny+1,0:nz+1))
  END IF

  IF (ALLOCATED(spex10)) THEN
    DEALLOCATE(spex10)
    ALLOCATE(spex10(nexchange+nexch_sec,0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(spex10(nexchange+nexch_sec,0:nx+1,0:ny+1,0:nz+1))
  END IF
  IF (ALLOCATED(spsurf10)) THEN
    DEALLOCATE(spsurf10)
    ALLOCATE(spsurf10(nsurf+nsurf_sec,0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(spsurf10(nsurf+nsurf_sec,0:nx+1,0:ny+1,0:nz+1))
  END IF
  IF (ALLOCATED(spgas10)) THEN
    DEALLOCATE(spgas10)
    ALLOCATE(spgas10(ngas,0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(spgas10(ngas,0:nx+1,0:ny+1,0:nz+1))
  END IF

  IF (ALLOCATED(t)) THEN
    DEALLOCATE(t)
    ALLOCATE(t(0:nx+1,0:ny+1,0:nz+1))
  ELSE
    ALLOCATE(t(0:nx+1,0:ny+1,0:nz+1))
  END IF
!!!  IF (ALLOCATED(por)) THEN
!!!    DEALLOCATE(por)
!!!    ALLOCATE(por(0:nx+1,0:ny+1,0:nz+1))
!!!  ELSE
!!!    ALLOCATE(por(0:nx+1,0:ny+1,0:nz+1))
!!!  END IF

END IF

IF (ALLOCATED(sn)) THEN
  DEALLOCATE(sn)
  ALLOCATE(sn(ncomp,nx,ny,nz))
ELSE
  ALLOCATE(sn(ncomp,nx,ny,nz))
END IF
IF (ALLOCATED(spold)) THEN
  DEALLOCATE(spold)
  ALLOCATE(spold(ncomp+nspec,nx,ny,nz))
ELSE
  ALLOCATE(spold(ncomp+nspec,nx,ny,nz))
END IF
IF (ALLOCATED(spexold)) THEN
  DEALLOCATE(spexold)
  ALLOCATE(spexold(nexchange+nexch_sec,nx,ny,nz))
ELSE
  ALLOCATE(spexold(nexchange+nexch_sec,nx,ny,nz))
END IF


IF (ALLOCATED(spgasold)) THEN
  DEALLOCATE(spgasold)
  ALLOCATE(spgasold(ngas,nx,ny,nz))
ELSE
  ALLOCATE(spgasold(ngas,nx,ny,nz))
END IF

IF (ALLOCATED(spsurfold)) THEN
  DEALLOCATE(spsurfold)
  ALLOCATE(spsurfold(nsurf+nsurf_sec,nx,ny,nz))
ELSE
  ALLOCATE(spsurfold(nsurf+nsurf_sec,nx,ny,nz))
END IF


IF (ALLOCATED(volfx)) THEN
  DEALLOCATE(volfx)
  ALLOCATE(volfx(nrct,0:nx,ny,nz))
ELSE
  ALLOCATE(volfx(nrct,0:nx,ny,nz))
END IF
IF (ALLOCATED(volSave)) THEN
  DEALLOCATE(volSave)
  ALLOCATE(volSave(nrct,0:nx,ny,nz))
ELSE
  ALLOCATE(volSave(nrct,0:nx,ny,nz))
END IF
volSave = 0.0d0
IF (ALLOCATED(volSaveByTimeStep)) THEN
  DEALLOCATE(volSaveByTimeStep)
  ALLOCATE(volSaveByTimeStep(101,nrct,0:nx,ny,nz))
ELSE
  ALLOCATE(volSaveByTimeStep(101,nrct,0:nx,ny,nz))
END IF
volSaveByTimeStep = 0.0d0
IF (ALLOCATED(VolumeLastTimeStep)) THEN
  DEALLOCATE(VolumeLastTimeStep)
  ALLOCATE(VolumeLastTimeStep(nrct,nx,ny,nz))
ELSE
  ALLOCATE(VolumeLastTimeStep(nrct,nx,ny,nz))
END IF
IF (ALLOCATED(rminSaveForDePaolo)) THEN
  DEALLOCATE(rminSaveForDePaolo)
  ALLOCATE(rminSaveForDePaolo(2,nx,ny,nz))
ELSE
  ALLOCATE(rminSaveForDePaolo(2,nx,ny,nz))
END IF
IF (ALLOCATED(dppt)) THEN
  DEALLOCATE(dppt)
  ALLOCATE(dppt(nrct,nx,ny,nz))
ELSE
  ALLOCATE(dppt(nrct,nx,ny,nz))
END IF
IF (ALLOCATED(area)) THEN
  DEALLOCATE(area)
  ALLOCATE(area(nrct,0:nx+1,ny,nz))
ELSE
  ALLOCATE(area(nrct,0:nx+1,ny,nz))
END IF
IF (ALLOCATED(jinit)) THEN
  DEALLOCATE(jinit)
  ALLOCATE(jinit(0:nx+1,0:ny+1,0:nz+1))
ELSE
  ALLOCATE(jinit(0:nx+1,0:ny+1,0:nz+1))
END IF
IF (ALLOCATED(xgram)) THEN
  DEALLOCATE(xgram)
  ALLOCATE(xgram(-1:nx+1,0:ny+1,0:nz+1))
ELSE
  ALLOCATE(xgram(-1:nx+1,0:ny+1,0:nz+1))
END IF
IF (ALLOCATED(xgramOld)) THEN
  DEALLOCATE(xgramOld)
  ALLOCATE(xgramOld(-1:nx+1,0:ny+1,0:nz+1))
ELSE
  ALLOCATE(xgramOld(-1:nx+1,0:ny+1,0:nz+1))
END IF

IF (ALLOCATED(H2Oreacted)) THEN
  DEALLOCATE(H2Oreacted)
  ALLOCATE(H2Oreacted(1:nx,1:ny,1:nz))
ELSE
  ALLOCATE(H2Oreacted(1:nx,1:ny,1:nz))
END IF

IF (Duan .OR. Duan2006) THEN
  IF (ALLOCATED(vrSave)) THEN
    DEALLOCATE(vrSave)
    ALLOCATE(vrSave(nx,ny,nz))
  ELSE
    ALLOCATE(vrSave(nx,ny,nz))
  END IF
END IF

H2Oreacted = 1.0d0

IF (ALLOCATED(told)) THEN
  DEALLOCATE(told)
  ALLOCATE(told(nx,ny,nz))
ELSE
  ALLOCATE(told(nx,ny,nz))
END IF


IF (ALLOCATED(raq_tot)) THEN
  DEALLOCATE(raq_tot)
  ALLOCATE(raq_tot(ikin,nx,ny,nz))
ELSE
  ALLOCATE(raq_tot(ikin,nx,ny,nz))
END IF
IF (ALLOCATED(ro)) THEN
  DEALLOCATE(ro)
  ALLOCATE(ro(-1:nx+2,-1:ny+2,-1:nz+2))
ELSE
  ALLOCATE(ro(-1:nx+2,-1:ny+2,-1:nz+2))
END IF
IF (ALLOCATED(roOld)) THEN
  DEALLOCATE(roOld)
  ALLOCATE(roOld(-1:nx+2,-1:ny+2,-1:nz+2))
ELSE
  ALLOCATE(roOld(-1:nx+2,-1:ny+2,-1:nz+2))
END IF
IF (ALLOCATED(rogas)) THEN
  DEALLOCATE(rogas)
  ALLOCATE(rogas(-1:nx+2,-1:ny+2,-1:nz+2))
ELSE
  ALLOCATE(rogas(-1:nx+2,-1:ny+2,-1:nz+2))
END IF
IF (ALLOCATED(exchangesites)) THEN
  DEALLOCATE(exchangesites)
  ALLOCATE(exchangesites(nexchange,nx,ny,nz))
ELSE
  ALLOCATE(exchangesites(nexchange,nx,ny,nz))
END IF

IF (ALLOCATED(MetabolicLagMineral)) THEN
  DEALLOCATE(MetabolicLagMineral)
  ALLOCATE(MetabolicLagMineral(1,nMonodBiomassMineral,nx,ny,nz))
ELSE
  ALLOCATE(MetabolicLagMineral(1,nMonodBiomassMineral,nx,ny,nz))
END IF
IF (ALLOCATED(SatLog)) THEN
  DEALLOCATE(SatLog)
  ALLOCATE(SatLog(ikin,nx,ny,nz))
ELSE
  ALLOCATE(SatLog(ikin,nx,ny,nz))
END IF
IF (ALLOCATED(tauZeroMineral)) THEN
  DEALLOCATE(tauZeroMineral)
  ALLOCATE(tauZeroMineral(1,nMonodBiomassMineral,nx,ny,nz))
ELSE
  ALLOCATE(tauZeroMineral(1,nMonodBiomassMineral,nx,ny,nz))
END IF
IF (ALLOCATED(MetabolicLagAqueous)) THEN
  DEALLOCATE(MetabolicLagAqueous)
  ALLOCATE(MetabolicLagAqueous(nMonodBiomassAqueous,nx,ny,nz))
ELSE
  ALLOCATE(MetabolicLagAqueous(nMonodBiomassAqueous,nx,ny,nz))
END IF
IF (ALLOCATED(tauZeroAqueous)) THEN
  DEALLOCATE(tauZeroAqueous)
  ALLOCATE(tauZeroAqueous(nMonodBiomassAqueous,nx,ny,nz))
ELSE
  ALLOCATE(tauZeroAqueous(nMonodBiomassAqueous,nx,ny,nz))
END IF
MetabolicLagMineral = 1.0d0
tauZeroMineral = 0.0d0
MetabolicLagAqueous = 1.0d0
tauZeroAqueous = 0.0d0

IF (KateMaher) THEN
  IF (ALLOCATED(muCalcium)) THEN
    DEALLOCATE(muCalcium)
    ALLOCATE(muCalcium(nx,ny,nz))   !!
  ELSE
    ALLOCATE(muCalcium(nx,ny,nz))   !!
  END IF
  IF (ALLOCATED(muCalciumBoundary)) THEN
    DEALLOCATE(muCalciumBoundary)
    ALLOCATE(muCalciumBoundary(6))   !!
  ELSE
    ALLOCATE(muCalciumBoundary(6))   !!
  END IF
  IF (ALLOCATED(muUranium234)) THEN
    DEALLOCATE(muUranium234)
    ALLOCATE(muUranium234(nx,ny,nz))   !!
  ELSE
    ALLOCATE(muUranium234(nx,ny,nz))   !!
  END IF
  IF (ALLOCATED(muUranium234Boundary)) THEN
    DEALLOCATE(muUranium234Boundary)
    ALLOCATE(muUranium234Boundary(6))   !!
  ELSE
    ALLOCATE(muUranium234Boundary(6))   !!
  END IF
  IF (ALLOCATED(muUranium238)) THEN
    DEALLOCATE(muUranium238)
    ALLOCATE(muUranium238(nx,ny,nz))   !!
  ELSE
    ALLOCATE(muUranium238(nx,ny,nz))   !!
  END IF
  IF (ALLOCATED(muUranium238Boundary)) THEN
    DEALLOCATE(muUranium238Boundary)
    ALLOCATE(muUranium238Boundary(6))   !!
  ELSE
    ALLOCATE(muUranium238Boundary(6))   !!
  END IF

  IF (ALLOCATED(muCalciumBulk)) THEN
    DEALLOCATE(muCalciumBulk)
    ALLOCATE(muCalciumBulk(nx,ny,nz))   !!
  ELSE
    ALLOCATE(muCalciumBulk(nx,ny,nz))   !!
  END IF
  IF (ALLOCATED(muUranium234Bulk)) THEN
    DEALLOCATE(muUranium234Bulk)
    ALLOCATE(muUranium234Bulk(nx,ny,nz))   !!
  ELSE
    ALLOCATE(muUranium234Bulk(nx,ny,nz))   !!
  END IF
  IF (ALLOCATED(muUranium238Bulk)) THEN
    DEALLOCATE(muUranium238Bulk)
    ALLOCATE(muUranium238Bulk(nx,ny,nz))   !!
  ELSE
    ALLOCATE(muUranium238Bulk(nx,ny,nz))   !!
  END IF

END IF

IF (ALLOCATED(decay_correct)) THEN
  DEALLOCATE(decay_correct)
  ALLOCATE(decay_correct(ncomp,nrct))
ELSE
  ALLOCATE(decay_correct(ncomp,nrct))
END IF
decay_correct = 1.0
IF (ALLOCATED(si)) THEN
  DEALLOCATE(si)
  ALLOCATE(si(nreactmax,nrct))
ELSE
  ALLOCATE(si(nreactmax,nrct))
END IF
IF (ALLOCATED(silog)) THEN
  DEALLOCATE(silog)
  ALLOCATE(silog(nreactmax,nrct))
ELSE
  ALLOCATE(silog(nreactmax,nrct))
END IF
IF (ALLOCATED(siln)) THEN
  DEALLOCATE(siln)
  ALLOCATE(siln(nreactmax,nrct))
ELSE
  ALLOCATE(siln(nreactmax,nrct))
END IF
IF (ALLOCATED(snorm)) THEN
  DEALLOCATE(snorm)
  ALLOCATE(snorm(nreactmax,nrct))
ELSE
  ALLOCATE(snorm(nreactmax,nrct))
END IF
IF (ALLOCATED(actenergy)) THEN
  DEALLOCATE(actenergy)
  ALLOCATE(actenergy(nreactmax,nrct))
ELSE
  ALLOCATE(actenergy(nreactmax,nrct))
END IF
IF (ALLOCATED(pre_rmin)) THEN
  DEALLOCATE(pre_rmin)
  ALLOCATE(pre_rmin(nreactmax,nrct))
ELSE
  ALLOCATE(pre_rmin(nreactmax,nrct))
END IF
IF (ALLOCATED(rmin)) THEN
  DEALLOCATE(rmin)
  ALLOCATE(rmin(nreactmax,nrct))
ELSE
  ALLOCATE(rmin(nreactmax,nrct))
END IF




IF (ALLOCATED(jac_rmin)) THEN
  DEALLOCATE(jac_rmin)
  ALLOCATE(jac_rmin(ncomp+nexchange+nsurf,nreactmax,nrct))
ELSE
  ALLOCATE(jac_rmin(ncomp+nexchange+nsurf,nreactmax,nrct))
END IF
IF (ALLOCATED(pre_raq)) THEN
  DEALLOCATE(pre_raq)
  ALLOCATE(pre_raq(nreactkinmax,ikin))
ELSE
  ALLOCATE(pre_raq(nreactkinmax,ikin))
END IF
IF (ALLOCATED(rdkin)) THEN
  DEALLOCATE(rdkin)
  ALLOCATE(rdkin(ikin,ncomp))
ELSE
  ALLOCATE(rdkin(ikin,ncomp))
END IF
IF (ALLOCATED(raq)) THEN
  DEALLOCATE(raq)
  ALLOCATE(raq(nreactkinmax,ikin))
ELSE
  ALLOCATE(raq(nreactkinmax,ikin))
END IF
IF (ALLOCATED(d_sp)) THEN
  DEALLOCATE(d_sp)
  ALLOCATE(d_sp(ncomp+nspec))
ELSE
  ALLOCATE(d_sp(ncomp+nspec))
END IF


END SUBROUTINE GlobalArrayAllocation
