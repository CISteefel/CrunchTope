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

    
SUBROUTINE ReadFlowField(nx,ny,nz,                                          &
  userC,userD,userP,user,amatpetsc,amatD,amatP,bvec,xvec,bvecD,xvecD,bvecP,xvecP)
USE crunchtype
USE params
USE runtime
USE crunch_interface
USE transport
USE medium
USE flow
USE temperature
USE ReadFlow

#include "petsc/finclude/petscmat.h"
USE petscmat

IMPLICIT NONE

!!  External variables and arrays

INTEGER(I4B), INTENT(IN)                               :: nx
INTEGER(I4B), INTENT(IN)                               :: ny
INTEGER(I4B), INTENT(IN)                               :: nz

!!  Internal variables and arrays

LOGICAL(LGT)                                               :: ext

INTEGER(I4B)                                               :: ls
INTEGER(I4B)                                               :: ierr
INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz

REAL(DP), PARAMETER                                        :: eps=1.D-12

REAL(DP)                                                   :: minsat

! ******************** PETSC declarations ********************************
PetscFortranAddr     userC(6),userD(6),userP(6),user(6)
Mat                  amatpetsc,amatD,amatP
Vec                  bvec,xvec,bvecD,xvecD,bvecP,xvecP
! ************************end PETSc declarations of PETSc variables ******

IF (ALLOCATED(satliqNuft)) THEN
  DEALLOCATE(satliqNuft)
END IF
ALLOCATE(satliqNuft(-1:nx+2,-1:ny+2,-1:nz+2))

IF (ALLOCATED(satliqNuftold)) THEN
  DEALLOCATE(satliqNuftold)
END IF
ALLOCATE(satliqNuftold(-1:nx+2,-1:ny+2,-1:nz+2))

IF (ALLOCATED(qxNuft)) THEN
  DEALLOCATE(qxNuft)
END IF
ALLOCATE(qxNuft(0:nx,ny,nz))

IF (ALLOCATED(qyNuft)) THEN
  DEALLOCATE(qyNuft)
END IF
ALLOCATE(qyNuft(nx,0:ny,nz))

IF (ALLOCATED(qzNuft)) THEN
  DEALLOCATE(qzNuft)
END IF
ALLOCATE(qzNuft(nx,ny,0:nz))

IF (ALLOCATED(qxNuftOld)) THEN
  DEALLOCATE(qxNuftOld)
END IF
ALLOCATE(qxNuftOld(0:nx,ny,nz))

IF (ALLOCATED(qyNuftOld)) THEN
  DEALLOCATE(qyNuftOld)
END IF
ALLOCATE(qyNuftOld(nx,0:ny,nz))

IF (ALLOCATED(qzNuftOld)) THEN
  DEALLOCATE(qzNuftOld)
END IF
ALLOCATE(qzNuftOld(nx,ny,0:nz))

IF (ALLOCATED(qgNuft)) THEN
  DEALLOCATE(qgNuft)
END IF
ALLOCATE(qgNuft(nx,ny,nz))

IF (ALLOCATED(qgNuftOld)) THEN
  DEALLOCATE(qgNuftOld)
END IF
ALLOCATE(qgNuftOld(nx,ny,nz))

IF (ALLOCATED(qxgasNuft)) THEN
  DEALLOCATE(qxgasNuft)
END IF
ALLOCATE(qxgasNuft(0:nx,ny,nz))

IF (ALLOCATED(qygasNuft)) THEN
  DEALLOCATE(qygasNuft)
END IF
ALLOCATE(qygasNuft(nx,0:ny,nz))

IF (ALLOCATED(qzgasNuft)) THEN
  DEALLOCATE(qzgasNuft)
END IF
ALLOCATE(qzgasNuft(nx,ny,0:nz))

IF (ALLOCATED(qxgasNuftOld)) THEN
  DEALLOCATE(qxgasNuftOld)
END IF
ALLOCATE(qxgasNuftOld(0:nx,ny,nz))

IF (ALLOCATED(qygasNuftOld)) THEN
  DEALLOCATE(qygasNuftOld)
END IF
ALLOCATE(qygasNuftOld(nx,0:ny,nz))

IF (ALLOCATED(qzgasNuftOld)) THEN
  DEALLOCATE(qzgasNuftOld)
END IF
ALLOCATE(qzgasNuftOld(nx,ny,0:nz))

IF (ALLOCATED(roNuft)) THEN
  DEALLOCATE(roNuft)
END IF
ALLOCATE(roNuft(-1:nx+2,-1:ny+2,-1:nz+2))

IF (ALLOCATED(roNuftOld)) THEN
  DEALLOCATE(roNuftOld)
END IF
ALLOCATE(roNuftOld(-1:nx+2,-1:ny+2,-1:nz+2))

qxNuft = 0.0
qyNuft = 0.0
qzNuft = 0.0
qxNuftOld = 0.0
qyNuftOld = 0.0
qzNuftOld = 0.0
qgNuft = 0.0
qgNuftOld = 0.0
qxgasNuft = 0.0
qygasNuft = 0.0
qzgasNuft = 0.0
qxgasNuftOld = 0.0
qygasNuftOld = 0.0
qzgasNuftOld = 0.0
roNuft = 0.0
roNuftOld = 0.0

CALL stringlen(NuftFile,ls)

INQUIRE(FILE=NuftFile,EXIST=ext)
IF (.NOT. ext) THEN
  CALL stringlen(NuftFile,ls)
  WRITE(*,*) 
  WRITE(*,*) ' NUFT velocity file not found: ', NuftFile(1:ls)
  WRITE(*,*)
!         ***** PETSc closeout**************
      IF (petscon) then
        call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
      END IF
      IF (CalculateFlow) then
        call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
      END IF
      IF (os3dpetsc) then
        call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
      END IF
      call PetscFinalize(ierr)
!       ****** PETSc closeout finished *********
    READ(*,*)
    STOP
END IF

!  Check for size of NUFT domain to make sure it matches

iunitNuft= 151
INQUIRE(FILE=NuftFile,EXIST=ext)
IF (.NOT. ext) THEN
  CALL stringlen(NuftFile,ls)
  WRITE(*,*) 
  WRITE(*,*) ' Nuft flow file not found: ', NuftFile(1:ls)
  WRITE(*,*)
!         ***** PETSc closeout**************
      IF (petscon) then
        call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
      END IF
      IF (CalculateFlow) then
        call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
      END IF
      IF (os3dpetsc) then
        call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
      END IF
      call PetscFinalize(ierr)
!       ****** PETSc closeout finished *********
    READ(*,*)
    STOP
END IF

OPEN(UNIT=iunitNuft,FILE=NuftFile,FORM='unformatted',STATUS='old',POSITION='append')
BACKSPACE iunitNuft
READ(iunitNuft) TotNuftSteps
NumNuftSteps = 1
REWIND(iunitNuft)

READ(iunitNuft) nxNuft,nyNuft,nzNuft
IF (nxNuft /= nx .OR. nyNuft /= ny .OR. nzNuft /= nz) THEN
  WRITE(*,*)
  WRITE(*,*) ' Dimensions in NUFT output file dont match CRUNCH dimensions'
  WRITE(*,*) ' Reading: ', NuftFile(1:ls)
  WRITE(*,*)
  WRITE(*,*) ' NUFT dimensions: '
  WRITE(*,*) '     NX = ',nxNuft
  WRITE(*,*) '     NY = ',nyNuft
  WRITE(*,*) '     NZ = ',nzNuft
  WRITE(*,*)
  WRITE(*,*) ' CRUNCH dimensions: '
  WRITE(*,*) '     NX = ',nx
  WRITE(*,*) '     NY = ',ny
  WRITE(*,*) '     NZ = ',nz
  WRITE(*,*)
!         ***** PETSc closeout**************
      IF (petscon) then
        call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
      END IF
      IF (CalculateFlow) then
        call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
      END IF
      IF (os3dpetsc) then
        call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
      END IF
      call PetscFinalize(ierr)
!       ****** PETSc closeout finished *********
    READ(*,*)
    STOP
END IF

! Now check to see which field variables will be read from NUFT

READ(iunitNuft) NumNuftVariables
IF (ALLOCATED(NuftVariable)) THEN
  DEALLOCATE(NuftVariable)
END IF
ALLOCATE(NuftVariable(NumNuftVariables))
READ(iunitNuft) NuftVariable

! Search for the various input options

ReadNuftLiqFlux = .FALSE.
ReadNuftSaturation = .FALSE.
ReadNuftPorosity = .FALSE.
ReadNuftTemperature = .FALSE.
ReadNuftLiqDensity = .FALSE.
ReadNuftGasFlux = .FALSE.
ReadNuftGasDensity = .FALSE.

DO i = 1,NumNuftVariables
  IF (NuftVariable(i) == 'qint.liquid') THEN
    ReadNuftLiqFlux = .TRUE.  
  END IF     
  IF (NuftVariable(i) == 'S.liquid') THEN
    ReadNuftSaturation = .TRUE.  
  END IF   
  IF (NuftVariable(i) == 'porosity') THEN
    ReadNuftPorosity = .TRUE.  
  END IF   
  IF (NuftVariable(i) == 'temperature') THEN
    ReadNuftTemperature = .TRUE.  
  END IF   
  IF (NuftVariable(i) == 'rho.liquid') THEN
    ReadNuftLiqDensity = .TRUE.  
  END IF   
  IF (NuftVariable(i) == 'rho.gas') THEN
    ReadNuftGasDensity = .TRUE.  
  END IF  
  IF (NuftVariable(i) == 'qint.gas') THEN
    ReadNuftGasFlux = .TRUE.  
  END IF   
END DO

! Check for source terms

READ(iunitNuft) NumSourceTerms
IF (ALLOCATED(jxNuftSource)) THEN
  DEALLOCATE(jxNuftSource)
END IF
ALLOCATE(jxNuftSource(NumSourceTerms))
IF (ALLOCATED(jyNuftSource)) THEN
  DEALLOCATE(jyNuftSource)
END IF
ALLOCATE(jyNuftSource(NumSourceTerms))
IF (ALLOCATED(jzNuftSource)) THEN
  DEALLOCATE(jzNuftSource)
END IF
ALLOCATE(jzNuftSource(NumSourceTerms))
IF (NumSourceTerms > 0) THEN
  READ(iunitNuft) jxNuftSource
  READ(iunitNuft) jyNuftSource
  READ(iunitNuft) jzNuftSource
END IF

READ(iunitNuft) timeNuft   

IF (timeNuft /= 0.0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Initial time from NUFT should be 0 '
  WRITE(*,*)
!         ***** PETSc closeout**************
      IF (petscon) then
        call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
      END IF
      IF (CalculateFlow) then
        call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
      END IF
      IF (os3dpetsc) then
        call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
      END IF
      call PetscFinalize(ierr)
!       ****** PETSc closeout finished *********
    READ(*,*)
    STOP
END IF

IF (NumSourceTerms > 0) THEN
  DO i = 1,NumSourceTerms
    READ(iunitNuft) qgNuft(jxNuftSource(i),jyNuftSource(i),jzNuftSource(i))
  END DO
  qgNuft = qgNuft*secyr             !  NUFT source term in kg/s --> convert to years
END IF
qgNuftOld = qgNuft
!!!qg = qgNuft

IF (ReadNuftLiqFlux) THEN
  READ(iunitNuft) qxNuft
  READ(iunitNuft) qyNuft
  READ(iunitNuft) qzNuft
  qxNuft = qxNuft*secyr             !  NUFT source term in kg/m**2/s --> convert to years
  qyNuft = qyNuft*secyr             !  NUFT source term in kg/m**2/s --> convert to years
  qzNuft = qzNuft*secyr             !  NUFT source term in kg/m**2/s --> convert to years
  qx = qxNuft
  qy = qyNuft
  qz = qzNuft
  qxNuftOld = qxNuft
  qyNuftOld = qyNuft
  qzNuftOld = qzNuft
END IF
IF (ReadNuftSaturation) THEN
  READ(iunitNuft) satliqNuft
  satliqNuftOld = satliqNuft
  satliq = satliqNuft
  DO jz = -1,nz+2
    DO jy = -1,ny+2
      DO jx = -1,nx+2
        IF (satliq(jx,jy,jz) < 0.005) THEN
          satliq(jx,jy,jz) = 0.005
        END IF
      END DO
    END DO
  END DO
  minsat = minval(DABS(satliq(1:nx,1:ny,1:nz)))
  WRITE(*,*)
  WRITE(*,*) ' Minsat = ',minsat
  WRITE(*,*)
  IF (minsat == 0.0) THEN
!         ***** PETSc closeout**************
      IF (petscon) then
        call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
      END IF
      IF (CalculateFlow) then
        call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
      END IF
      IF (os3dpetsc) then
        call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
      END IF
      call PetscFinalize(ierr)
!       ****** PETSc closeout finished *********
      READ(*,*)
      STOP
  END IF
  IF (minsat+eps < 1.0d0) THEN
    WRITE(*,*)
    WRITE(*,*) ' NUFT read indicates this is an unsaturated problem'
    WRITE(*,*)
    isaturate = 1
  END IF
END IF
IF (ReadNuftPorosity) THEN
  READ(iunitNuft) por
END IF
IF (ReadNuftTemperature) THEN
  READ(iunitNuft) t
END IF
IF (ReadNuftLiqDensity) THEN
  READ(iunitNuft) roNuft      
  ro = roNuft
  roNuftOld = roNuft
END IF
IF (ReadNuftGasDensity) THEN
  READ(iunitNuft) rogas
END IF
IF (ReadNuftGasFlux) THEN
  READ(iunitNuft) qxgasNuft
  READ(iunitNuft) qygasNuft
  READ(iunitNuft) qzgasNuft
END IF
qxgasNuft = qxgasNuft*secyr
qygasNuft = qygasNuft*secyr
qzgasNuft = qzgasNuft*secyr
qxgas = qxgasNuft
qygas = qygasNuft
qzgas = qzgasNuft
qxgasNuftOld = qxgasNuft
qygasNuftOld = qygasNuft
qzgasNuftOld = qzgasNuft


RETURN
END SUBROUTINE ReadFlowField
