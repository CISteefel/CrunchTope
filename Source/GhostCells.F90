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
    
SUBROUTINE GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,xxx,TEXT)
USE crunchtype

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(IN)                                    :: lowX
INTEGER(I4B), INTENT(IN)                                    :: lowY
INTEGER(I4B), INTENT(IN)                                    :: lowZ
INTEGER(I4B), INTENT(IN)                                    :: highX
INTEGER(I4B), INTENT(IN)                                    :: highY
INTEGER(I4B), INTENT(IN)                                    :: highZ

REAL(DP), DIMENSION(lowX:highX,lowY:highY,lowZ:highZ), INTENT(INOUT)                   :: xxx

CHARACTER (LEN=15), INTENT(IN)                              :: text

!  Internal variables and arrays

INTEGER(I4B)                                                :: nxcheck
INTEGER(I4B)                                                :: nycheck
INTEGER(I4B)                                                :: nzcheck
INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                     :: xdummy

IF (lowX /= -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Lower bounds in X direction should be -1 for: ', text
  WRITE(*,*) 
  STOP
END IF
IF (lowY /= -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Lower bounds in Y direction should be -1 for: ', text
  WRITE(*,*) 
  STOP
END IF
IF (lowZ /= -1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Lower bounds in Z direction should be -1 for: ', text
  WRITE(*,*) 
  STOP
END IF
IF (highX /= nx+2) THEN
  WRITE(*,*)
  WRITE(*,*) ' Upper bounds in X direction should be nx+2 for: ', text
  WRITE(*,*) 
  STOP
END IF
IF (highY /= ny+2) THEN
  WRITE(*,*)
  WRITE(*,*) ' Upper bounds in Y direction should be ny+2 for: ', text
  WRITE(*,*) 
  STOP
END IF
IF (highZ /= nz+2) THEN
  WRITE(*,*)
  WRITE(*,*) ' Upper bounds in X direction should be nz+2 for: ', text
  WRITE(*,*) 
  STOP
END IF

xxx(-1,:,:) = xxx(1,:,:)
xxx(0,:,:) = xxx(1,:,:)
xxx(nx+2,:,:) = xxx(nx,:,:)
xxx(nx+1,:,:) = xxx(nx,:,:)

xxx(:,-1,:) = xxx(:,1,:)
xxx(:,0,:) = xxx(:,1,:)
xxx(:,ny+2,:) = xxx(:,ny,:)
xxx(:,ny+1,:) = xxx(:,ny,:)

xxx(:,:,-1) = xxx(:,:,1)
xxx(:,:,0) = xxx(:,:,1)
xxx(:,:,nz+2) = xxx(:,:,nz)
xxx(:,:,nz+1) = xxx(:,:,nz)

RETURN

  END SUBROUTINE GhostCells
  
  
  
  
SUBROUTINE GhostCells_Richards(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,xxx, text)
USE crunchtype

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(IN)                                    :: lowX
INTEGER(I4B), INTENT(IN)                                    :: lowY
INTEGER(I4B), INTENT(IN)                                    :: lowZ
INTEGER(I4B), INTENT(IN)                                    :: highX
INTEGER(I4B), INTENT(IN)                                    :: highY
INTEGER(I4B), INTENT(IN)                                    :: highZ

REAL(DP), DIMENSION(lowX:highX,lowY:highY,lowZ:highZ), INTENT(INOUT)                   :: xxx

CHARACTER (LEN=15), INTENT(IN)                              :: text

!  Internal variables and arrays

INTEGER(I4B)                                                :: nxcheck
INTEGER(I4B)                                                :: nycheck
INTEGER(I4B)                                                :: nzcheck
INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                     :: xdummy


IF (lowX /= 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Lower bounds in X direction should be 0 for: ', text
    WRITE(*,*) 
    STOP
END IF

IF (highX /= nx+1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Upper bounds in X direction should be nx+1 for: ', text
  WRITE(*,*) 
  STOP
END IF

xxx(0,:,:) = xxx(1,:,:)
xxx(nx+1,:,:) = xxx(nx,:,:)
  
IF (nx > 1 .AND. ny > 1) THEN ! more than two-dimensional problem
  
  IF (lowY /= 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Lower bounds in Y direction should be 0 for: ', text
    WRITE(*,*) 
    STOP
  END IF
  
  IF (highY /= ny+1) THEN
    WRITE(*,*)
    WRITE(*,*) ' Upper bounds in Y direction should be ny+1 for: ', text
    WRITE(*,*) 
    STOP
  END IF
  
  xxx(:,0,:) = xxx(:,1,:)
  xxx(:,ny+1,:) = xxx(:,ny,:)

END IF

IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
  
  IF (lowZ /= 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Lower bounds in Z direction should be 0 for: ', text
    WRITE(*,*) 
    STOP
  END IF
  
  
  IF (highZ /= nz+1) THEN
    WRITE(*,*)
    WRITE(*,*) ' Upper bounds in X direction should be nz+1 for: ', text
    WRITE(*,*) 
    STOP
  END IF
  
  xxx(:,:,0) = xxx(:,:,1)
  xxx(:,:,nz+1) = xxx(:,:,nz)
END IF

RETURN

END SUBROUTINE GhostCells_Richards
