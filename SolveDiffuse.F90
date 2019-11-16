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

    
SUBROUTINE SolveDiffuse(nx,ny,nz,nn,icomp,delt,user,amatD)
USE crunchtype
USE params
USE concentration
USE solver
USE medium
USE transport
USE temperature

#include "petsc/finclude/petscmat.h"
USE petscmat
 
IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                              :: nx
INTEGER(I4B), INTENT(IN)                                              :: ny
INTEGER(I4B), INTENT(IN)                                              :: nz
INTEGER(I4B), INTENT(IN)                                              :: nn
INTEGER(I4B), INTENT(IN)                                              :: icomp

REAL(DP), INTENT(IN)                                                  :: delt

!  Internal variables and arrays

INTEGER(I4B)                                                          :: jx
INTEGER(I4B)                                                          :: jy
INTEGER(I4B)                                                          :: jz
INTEGER(I4B)                                                          :: j
INTEGER(I4B)                                                          :: i
INTEGER(I4B)                                                          :: ierr
INTEGER(I4B)                                                          :: itsiterate
INTEGER(I4B)                                                          :: nxyz
     
REAL(DP)                                                              :: AccumulationTerm
REAL(DP)                                                              :: RightHandSide
REAL(DP)                                                              :: DiagonalTerm

! *******************begin PETSc declarations of f90 variables***********
INTEGER(I4B)             ::numprocs
INTEGER(I4B)             ::irank
INTEGER(I4B)             ::linefil
INTEGER(I4B), PARAMETER  ::maxitsksp=100
REAL(DP), PARAMETER      ::zero=0.0d0

!*********************end PETSc declarations ******************************

! ******************** PETSC declarations ********************************
PetscFortranAddr     user(6)
Mat                  amatD
! ************************end PETSc declarations of PETSc variables ******

IF (nn == 0) THEN
  RETURN
END IF

IF (icomp > 1) GOTO 500               !  No need to calculate matrix, since matrix has been calculated for component 1

call MatZeroEntries(amatD,ierr)

IF (nx > 1 .AND. ny ==1 .AND. nz == 1) THEN           ! 1D problem assuming jx is coordinate

  jy = 1
  jz = 1
  DO jx = 2,nx-1
    j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1      
    IF (activecell(jx,jy,jz) == 0) THEN
      DiagonalTerm = 1.0d0
      CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
      CALL MatSetValues(amatD,1,j,1,j-1,0.0d0,INSERT_VALUES,ierr)
      CALL MatSetValues(amatD,1,j,1,j+1,0.0d0,INSERT_VALUES,ierr)
    ELSE
      AccumulationTerm = dxy(jx,jy,jz)*ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)/delt
      DiagonalTerm = bDD(jx,jy,jz) + AccumulationTerm
      CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
      CALL MatSetValues(amatD,1,j,1,j-1,aDD(jx,jy,jz),INSERT_VALUES,ierr)
      CALL MatSetValues(amatD,1,j,1,j+1,cDD(jx,jy,jz),INSERT_VALUES,ierr)
    END IF
  END DO

  jx = 1
  j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
  IF (activecell(jx,jy,jz) == 0) THEN
    DiagonalTerm = 1.0d0
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j+1,0.0d0,INSERT_VALUES,ierr)
  ELSE
    AccumulationTerm = dxy(jx,jy,jz)*ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)/delt
    DiagonalTerm = bDD(jx,jy,jz) + AccumulationTerm
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j+1,cDD(jx,jy,jz),INSERT_VALUES,ierr)
  END IF

  jx = nx
  j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
  IF (activecell(jx,jy,jz) == 0) THEN
    DiagonalTerm = 1.0d0
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j-1,0.0d0,INSERT_VALUES,ierr)
  ELSE
    AccumulationTerm = dxy(jx,jy,jz)*ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)/delt
    DiagonalTerm = bDD(jx,jy,jz) + AccumulationTerm
    CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)  
    CALL MatSetValues(amatD,1,j,1,j-1,aDD(jx,jy,jz),INSERT_VALUES,ierr)
  END IF

ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN       ! 3D problem

  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 2,nx-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecell(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatD,1,j,1,j-1,0.0d0,INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+1,0.0d0,INSERT_VALUES,ierr)
        ELSE
          CALL MatSetValues(amatD,1,j,1,j-1,aDD(jx,jy,jz),INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+1,cDD(jx,jy,jz),INSERT_VALUES,ierr)
        END IF
      END DO
      jx = 1
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j+1,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j+1,cDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
      jx = nx
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j-1,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j-1,aDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF

    END DO
  END DO

  DO jz = 1,nz
    DO jx = 1,nx
      DO jy = 2,ny-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecell(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatD,1,j,1,j-nx,0.0d0,INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+nx,0.0d0,INSERT_VALUES,ierr)
        ELSE
          CALL MatSetValues(amatD,1,j,1,j-nx,fDD(jx,jy,jz),INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+nx,dDD(jx,jy,jz),INSERT_VALUES,ierr)
        END IF
      END DO
      jy = 1
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j+nx,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j+nx,dDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
      jy = ny
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j-nx,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j-nx,fDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
    END DO
  END DO

  DO jy = 1,ny
    DO jx = 1,nx
      DO jz = 2,nz-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecell(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatD,1,j,1,j-nx*ny,0.0d0,INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+nx*ny,0.0d0,INSERT_VALUES,ierr)
        ELSE
          CALL MatSetValues(amatD,1,j,1,j-nx*ny,iDD(jx,jy,jz),INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+nx*ny,gDD(jx,jy,jz),INSERT_VALUES,ierr)
        END IF
      END DO
      jz = 1
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j+nx*ny,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j+nx*ny,gDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
      jz = nz
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j-nx*ny,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j-nx*ny,iDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
    END DO
  END DO

  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecell(jx,jy,jz) == 0) THEN
          DiagonalTerm = 1.0d0
          CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)
        ELSE
          AccumulationTerm = dxy(jx,jy,jz)*ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)/delt
          DiagonalTerm = bDD(jx,jy,jz) + eDD(jx,jy,jz) + hDD(jx,jy,jz) + AccumulationTerm
          CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)
        END IF

      END DO
    END DO
  END DO

ELSE                                                !  2D problem

  IF (nz > 1 .AND. ny == 1) THEN

    jy = 1
    DO jz = 1,nz
      DO jx = 2,nx-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecell(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatD,1,j,1,j-1,0.0d0,INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+1,0.0d0,INSERT_VALUES,ierr)
        ELSE
          CALL MatSetValues(amatD,1,j,1,j-1,aDD(jx,jy,jz),INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+1,cDD(jx,jy,jz),INSERT_VALUES,ierr)
        END IF
      END DO
      jx = 1
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j+1,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j+1,cDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
      jx = nx
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j-1,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j-1,aDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
    END DO

    jy = 1
    DO jx = 1,nx
      DO jz = 2,nz-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecell(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatD,1,j,1,j-nx*ny,0.0d0,INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+nx*ny,0.0d0,INSERT_VALUES,ierr)
        ELSE
          CALL MatSetValues(amatD,1,j,1,j-nx*ny,iDD(jx,jy,jz),INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+nx*ny,gDD(jx,jy,jz),INSERT_VALUES,ierr)
        END IF
      END DO
      jz = 1
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j+nx*ny,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j+nx*ny,gDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
      jz = nz
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j-nx*ny,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j-nx*ny,iDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
    END DO

    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
          IF (activecell(jx,jy,jz) == 0) THEN
            DiagonalTerm = 1.0d0
            CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)
          ELSE
            AccumulationTerm = dxy(jx,jy,jz)*ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)/delt
            DiagonalTerm = bDD(jx,jy,jz) + hDD(jx,jy,jz) + AccumulationTerm
            CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)
          END IF
        END DO
      END DO
    END DO
 
  ELSE

    jz = 1
    DO jy = 1,ny
      DO jx = 2,nx-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecell(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatD,1,j,1,j-1,0.0d0,INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+1,0.0d0,INSERT_VALUES,ierr)
        ELSE
          CALL MatSetValues(amatD,1,j,1,j-1,aDD(jx,jy,jz),INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+1,cDD(jx,jy,jz),INSERT_VALUES,ierr)
        END IF
      END DO
      jx = 1
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j+1,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j+1,cDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
      jx = nx
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j-1,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j-1,aDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
    END DO

    jz = 1
    DO jx = 1,nx
      DO jy = 2,ny-1
        j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
        IF (activecell(jx,jy,jz) == 0) THEN
          CALL MatSetValues(amatD,1,j,1,j-nx,0.0d0,INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+nx,0.0d0,INSERT_VALUES,ierr)
        ELSE
          CALL MatSetValues(amatD,1,j,1,j-nx,fDD(jx,jy,jz),INSERT_VALUES,ierr)
          CALL MatSetValues(amatD,1,j,1,j+nx,dDD(jx,jy,jz),INSERT_VALUES,ierr)
        END IF
      END DO
      jy = 1
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j+nx,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j+nx,dDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
      jy = ny
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
      IF (activecell(jx,jy,jz) == 0) THEN
        CALL MatSetValues(amatD,1,j,1,j-nx,0.0d0,INSERT_VALUES,ierr)
      ELSE
        CALL MatSetValues(amatD,1,j,1,j-nx,fDD(jx,jy,jz),INSERT_VALUES,ierr)
      END IF
    END DO

    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1 
          IF (activecell(jx,jy,jz) == 0) THEN
            DiagonalTerm = 1.0d0
            CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)
          ELSE
            AccumulationTerm = dxy(jx,jy,jz)*ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)/delt
            DiagonalTerm = bDD(jx,jy,jz) + eDD(jx,jy,jz) + AccumulationTerm
            CALL MatSetValues(amatD,1,j,1,j,DiagonalTerm,INSERT_VALUES,ierr)
          END IF
        END DO
      END DO
    END DO
  
  END IF

END IF

500 CONTINUE

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
      IF (activecell(jx,jy,jz) == 0) THEN
        RightHandSide = sn(icomp,jx,jy,jz)
      ELSE
        RightHandSide = dxy(jx,jy,jz)*ro(jx,jy,jz)*por(jx,jy,jz)*satliq(jx,jy,jz)*sn(icomp,jx,jy,jz)/delt
      END IF
      BvecCrunchD(j) = RightHandSide
    END DO
  END DO
END DO

CALL MatAssemblyBegin(amatD,MAT_FINAL_ASSEMBLY,ierr)
CALL MatAssemblyEnd(amatD,MAT_FINAL_ASSEMBLY,ierr)

RETURN
END SUBROUTINE SolveDiffuse
