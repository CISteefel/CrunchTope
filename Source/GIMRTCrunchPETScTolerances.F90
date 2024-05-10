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
    
subroutine GIMRTCrunchPETScTolerances(user,rtolksp,atolksp,dtolksp,maxitsksp,ierr)
USE crunchtype
USE solver
USE mpi
 
#include "petsc/finclude/petsc.h"
USE petscmat
USE petscksp

IMPLICIT NONE

!  External variables and arrays
REAL(DP), INTENT(IN)                                                  :: rtolksp
REAL(DP), INTENT(IN)                                                  :: atolksp
REAL(DP), INTENT(IN)                                                  :: dtolksp

INTEGER(I4B), INTENT(IN)                                              :: maxitsksp
INTEGER(I4B), INTENT(IN OUT)                                          :: ierr

!  Internal variables and arrays

! ******************** PETSC declarations ********************************
PetscFortranAddr     user(*)
!!SLES                 sles
PC                   pc
KSP                  ksp
! ************************end PETSc declarations of PETSc variables ******

pc%v = user(5)
ksp%v = user(6)
 
! Tolerances for linear solver set here

!!!call KSPSetTolerances(ksp,rtolksp,atolksp,dtolksp,maxitsksp,ierr)
call KSPSetTolerances(ksp,rtolksp,PETSC_DEFAULT_REAL,PETSC_DEFAULT_REAL,maxitsksp,ierr)

! Choose linear solver

IF (GIMRT_SolverMethod == 'bcgs') THEN
  CALL KSPSetType(ksp,KSPBCGS,ierr)
ELSE 
  CALL KSPSetType(ksp,KSPLGMRES,ierr)
!!!  CALL KSPGMRESSetRestart(ksp,30)
!!!  CALL KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
!!!  CALL KSPGMRESSetHapTol(ksp,1.0D-12)
END IF

! Choose preconditioning method

IF (GIMRT_PCMethod == 'ilu') THEN
  call KSPGetPC(ksp,pc,ierr)
  CALL PCSetType(pc,PCILU,ierr)
  CALL PCFactorSetLevels(pc,GIMRTlevel,ierr)
ELSE IF (GIMRT_PCMethod == 'lu') THEN
  call KSPGetPC(ksp,pc,ierr)
  CALL PCSetType(pc,PCLU,ierr)
  CALL PCFactorSetLevels(pc,GIMRTlevel,ierr)
ELSE
  call KSPGetPC(ksp,pc,ierr)
  CALL PCSetType(pc,PCBJACOBI,ierr)
!!!  CALL PCBJacobiSetLocalBlocks(pc,40,ierr)
  CALL PCFactorSetLevels(pc,GIMRTlevel,ierr)
END IF

RETURN
END subroutine GIMRTCrunchPETScTolerances

