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
    
SUBROUTINE AllocateALL(ncomp,nspec,ngas,nrct,nexchange,nsurf,nsurf_sec,npot,neqn,nx,ny,nz)
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
USE NanoCrystal
USE isotope, ONLY: IsotopeMineralRare,IsotopeMineralCommon

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

IF (ALLOCATED(sion)) THEN
  DEALLOCATE(sion)
END IF
ALLOCATE(sion(nx,ny,nz))
sion = 0.0
IF (ALLOCATED(fxmax)) THEN
  DEALLOCATE(fxmax)
END IF
ALLOCATE(fxmax(neqn))
fxmax = 0.0
IF (ALLOCATED(sumrd)) THEN
  DEALLOCATE(sumrd)
END IF
ALLOCATE(sumrd(ncomp+nexchange+nsurf))
sumrd = 0.0
IF (ALLOCATED(sumjackin)) THEN
  DEALLOCATE(sumjackin)
END IF
ALLOCATE(sumjackin(ncomp))
sumjackin = 0.0
IF (ALLOCATED(surf)) THEN
  DEALLOCATE(surf)
END IF
ALLOCATE(surf(nreactmax,nrct))
surf = 0.0
IF (ALLOCATED(jac_sat)) THEN
  DEALLOCATE(jac_sat)
END IF
ALLOCATE(jac_sat(ncomp+nexchange+nsurf))
jac_sat = 0.0
IF (ALLOCATED(jac_rateFactor)) THEN
  DEALLOCATE(jac_rateFactor)
END IF
ALLOCATE(jac_rateFactor(ncomp+nexchange+nsurf))
jac_rateFactor = 0.0
IF (ALLOCATED(jac_pre)) THEN
  DEALLOCATE(jac_pre)
END IF
ALLOCATE(jac_pre(ncomp+nexchange+nsurf,nreactmax))
jac_pre = 0.0
IF (ALLOCATED(jac_preKin)) THEN
  DEALLOCATE(jac_preKin)
END IF
ALLOCATE(jac_preKin(ncomp,nreactkinmax))
jac_preKin = 0.0

IF (ALLOCATED(sppTMP)) THEN
  DEALLOCATE(sppTMP)
END IF
ALLOCATE(sppTMP(ncomp+nspec))
sppTMP = 0.0
IF (ALLOCATED(sppTMP10)) THEN
  DEALLOCATE(sppTMP10)
END IF
ALLOCATE(sppTMP10(ncomp+nspec))
sppTMP10 = 0.0

IF (ALLOCATED(sppTMPperturb)) THEN
  DEALLOCATE(sppTMPperturb)
END IF
ALLOCATE(sppTMPperturb(ncomp+nspec))
sppTMPperturb = 0.0
IF (ALLOCATED(sppTMP10perturb)) THEN
  DEALLOCATE(sppTMP10perturb)
END IF
ALLOCATE(sppTMP10perturb(ncomp+nspec))
sppTMP10perturb = 0.0


!!  Allocation of Crystal Size Distribution stuff

!!CALL AllocateCSD(ncomp,nrct,nx,ny,nz)

!!IF (JennyRifle) THEN
!!
!!  IF (ALLOCATED(tauZero)) THEN
!!    DEALLOCATE(tauZero)
!!  END IF
!!  ALLOCATE(tauZero(nx,ny,nz))
!!  tauZero = 0.0d0
!!
!!  IF (ALLOCATED(MetabolicLag)) THEN
!!    DEALLOCATE(MetabolicLag)
!!  END IF
!!  ALLOCATE(MetabolicLag(nx,ny,nz))
!!  MetabolicLag = 0.0d0
!!
!!END IF

RETURN
END SUBROUTINE AllocateALL