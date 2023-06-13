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
    
SUBROUTINE AllocateGasesGimrt(nx,ny,nz,ncomp)

USE crunchtype
USE params
USE runtime
USE concentration
USE transport
USE solver
USE ReadFlow

IMPLICIT NONE

!!  External variables and arrays

INTEGER(I4B), INTENT(IN)                               :: nx
INTEGER(I4B), INTENT(IN)                               :: ny
INTEGER(I4B), INTENT(IN)                               :: nz
INTEGER(I4B), INTENT(IN)                               :: ncomp

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

RETURN
END SUBROUTINE AllocateGasesGimrt