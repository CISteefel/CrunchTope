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
subroutine CrunchPETScInitializePressure(nx,ny,nz,userP,ierr,xvecP,bvecP,amatP)
USE crunchtype
USE flow, ONLY:  XvecCrunchP, BvecCrunchP

#include "petsc/finclude/petscmat.h"
#include "petsc/finclude/petscksp.h"
USE petscmat
USE petscksp
 
IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                              :: nx
INTEGER(I4B), INTENT(IN)                                              :: ny
INTEGER(I4B), INTENT(IN)                                              :: nz
INTEGER(I4B), INTENT(IN OUT)                                          :: ierr

!  Internal variables and arrays

INTEGER(I4B)                                                          :: nxyz
INTEGER(I4B)                                                          :: linefil

! ******************** PETSC declarations ********************************
PetscFortranAddr     userP(*)
Mat                  amatP
Vec                  bvecP,xvecP
!!SLES                 sles
PC                   pc
KSP                  ksp
!!Scalar               zeroPetsc
! ************************end PETSc declarations of PETSc variables ******

nxyz = nx*ny*nz

IF (ALLOCATED(XvecCrunchP)) THEN
  DEALLOCATE(XvecCrunchP)
END IF
ALLOCATE(XvecCrunchP(0:nxyz-1))
IF (ALLOCATED(BvecCrunchP)) THEN
  DEALLOCATE(BvecCrunchP)
END IF
ALLOCATE(BvecCrunchP(0:nxyz-1))

IF (ny == 1 .AND. nz == 1) THEN       !  1D problem (assumes an X coordinate direction
  linefil = 3
ELSE
  IF (ny > 1 .AND. nz > 1) THEN      !  3D problem
    linefil = 7
  ELSE
    linefil = 5
  END IF
END IF

call MatCreateSeqAIJ(PETSC_COMM_SELF,nxyz,nxyz,linefil,PETSC_NULL_INTEGER,amatP,ierr)
!!call MatSetOption(amatP,MAT_COLUMN_ORIENTED,ierr)
call MatSetOption(amatP,MAT_ROW_ORIENTED,PETSC_FALSE,ierr)
call VecCreateSeqWithArray(PETSC_COMM_SELF,1,nxyz,BvecCrunchP,bvecP,ierr)
call VecCreateSeqWithArray(PETSC_COMM_SELF,1,nxyz,XvecCrunchP,xvecP,ierr)
call KSPCreate(PETSC_COMM_SELF,ksp,ierr)
call KSPSetFromOptions(ksp,ierr)
call KSPGetPC(ksp,pc,ierr)
!!call SLESGetKSP(sles,ksp,ierr)

userP(1) = amatP%v
userP(2) = bvecP%v
userP(3) = xvecP%v
userP(5) = pc%v
userP(6) = ksp%v

RETURN
END subroutine CrunchPETScInitializePressure
