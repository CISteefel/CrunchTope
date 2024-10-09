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

    
SUBROUTINE restart(time,nn,nint,nexchange,nsurf,nrct,nx,ny,nz,nstop,nstopsave, &
     delt,dtold,tstep,deltmin,dtmaxcour,dtmax,userC,userD,userP,user,     &
     amatpetsc,amatD,amatP,bvec,xvec,bvecD,xvecD,bvecP,xvecP)
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
USE modflowModule

#include "petsc/finclude/petscmat.h"
USE petscmat

IMPLICIT NONE

REAL(DP), INTENT(OUT)                         :: time
INTEGER(I4B), INTENT(OUT)                     :: nn
INTEGER(I4B), INTENT(OUT)                     :: nint
INTEGER(I4B), INTENT(IN)                      :: nexchange
INTEGER(I4B), INTENT(IN)                      :: nsurf
INTEGER(I4B), INTENT(IN)                      :: nrct
INTEGER(I4B), INTENT(IN)                      :: nx
INTEGER(I4B), INTENT(IN)                      :: ny
INTEGER(I4B), INTENT(IN)                      :: nz
INTEGER(I4B), INTENT(INOUT)                   :: nstop
INTEGER(I4B), INTENT(OUT)                     :: nstopsave
REAL(DP), INTENT(INOUT)                       :: delt
REAL(DP), INTENT(INOUT)                       :: dtold
REAL(DP), INTENT(INOUT)                       :: tstep
REAL(DP), INTENT(INOUT)                       :: deltmin
REAL(DP), INTENT(INOUT)                       :: dtmaxcour
REAL(DP), INTENT(INOUT)                       :: dtmax

INTEGER(I4B)                                :: ierr
INTEGER(I4B)                                :: iures
INTEGER(I4B)                                :: ls
INTEGER(I4B)                                :: ncount

LOGICAL(LGT)                                :: ext
LOGICAL(LGT)                                :: TrueFalse

REAL(DP)                                    :: DummyReal
REAL(DP), DIMENSION(:), ALLOCATABLE         :: tempreal
REAL(DP), DIMENSION(:), ALLOCATABLE         :: RealDummyArray
INTEGER(I4B), DIMENSION(:), ALLOCATABLE         :: IntegerDummyArray
INTEGER(I4B)                                :: nxyz


! ******************** PETSC declarations ********************************
PetscFortranAddr     userC(6),userD(6),userP(6),user(6)
Mat                  amatpetsc,amatD,amatP
Vec                  bvec,xvec,bvecD,xvecD,bvecP,xvecP
! ************************end PETSc declarations of PETSc variables ******

INQUIRE(FILE=restartfile,EXIST=ext)
IF (.NOT. ext) THEN
  CALL stringlen(restartfile,ls)
  WRITE(*,*) 
  WRITE(*,*) ' Restart file not found: ', restartfile(1:ls)
  WRITE(*,*)
!         ***** PETSc Closeout*****************
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
!    ****** PETSc closeout finished *********
      READ(*,*)
      STOP
    END IF
!!   *******  PETSc Closeout   *****************

iures = 131
INQUIRE(FILE=restartfile,OPENED=truefalse)
IF (truefalse) then
  CLOSE(UNIT=iures)
  OPEN(UNIT=iures,FILE=restartfile,FORM='unformatted',STATUS='old',ERR=6001)
ELSE
  OPEN(UNIT=iures,FILE=restartfile,FORM='unformatted',STATUS='old',ERR=6001)
END IF

    READ(iures) time
    READ(iures) nn
    READ(iures) nint
    READ(iures) DummyReal,dtold,DummyReal,DummyReal,DummyReal,dtmax
    READ(iures) keqaq
    READ(iures) keqgas
    READ(iures) keqsurf
    READ(iures) xgram
    READ(iures) spnO2
    READ(iures) spnnO2
    READ(iures) sp
    READ(iures) s
    READ(iures) sn
    READ(iures) sp10
    READ(iures) spold
    READ(iures) spex
    READ(iures) spex10
    READ(iures) lngamma
    READ(iures) exchangesites
    READ(iures) spexold
    READ(iures) spgas
    READ(iures) spgasold
    READ(iures) spgas10
    IF (isaturate==1) then
      READ(iures) sgas
      READ(iures) sgasn
    ENDIF
    if (ierode==1) then
      READ(iures) ssurf
      READ(iures) ssurfn
    endif
  
    
    READ(iures) sexold
    READ(iures) ssurfold
    READ(iures) spsurf
    READ(iures) spsurf10
    READ(iures) spsurfold 
    READ(iures) raq_tot
    READ(iures) sion
    IF (ALLOCATED(IntegerDummyArray)) THEN
      DEALLOCATE(IntegerDummyArray)
    END IF
    ALLOCATE(IntegerDummyArray(nxyz))
    READ(iures) IntegerDummyArray
!!!    READ(iures) jinit
    READ(iures) keqmin
    READ(iures) volfx
    READ(iures) dppt
    READ(iures) area
    IF (ALLOCATED(RealDummyArray)) THEN
      DEALLOCATE(RealDummyArray)
    END IF
    ALLOCATE(RealDummyArray(nrct*nxyz))
    READ(iures) RealDummyArray
    READ(iures) RealDummyArray
    READ(iures) RealDummyArray
!!!    READ(iures) areainByGrid
!!!    READ(iures) volinByGrid
!!!    READ(iures) specificByGrid
    READ(iures) LogPotential
    READ(iures) t
    READ(iures) told
    READ(iures) ro   
    READ(iures) por 
    READ(iures) satliq
    READ(iures) qxgas
    READ(iures) qygas
    READ(iures) qzgas
    READ(iures) pres
    READ(iures) ActiveCell
    READ(iures) VolSaveByTimeStep
    READ(iures) Volsave
    READ(iures) ncounter
    
  
  porin = por

CLOSE(UNIT=iures,STATUS='keep')

IF (ALLOCATED(RealDummyArray)) THEN
  DEALLOCATE(RealDummyArray)
END IF
IF (ALLOCATED(IntegerDummyArray)) THEN
  DEALLOCATE(IntegerDummyArray)
END IF

! Update file counter 

IF (AppendRestart) THEN
  IF (ALLOCATED(tempreal)) THEN
    DEALLOCATE(tempreal)
  END IF
  ALLOCATE(tempreal(nstop))
  tempreal = prtint
  DEALLOCATE(prtint)
  
  nstopsave = nstop
  nstop = nstop + nint-1
  ALLOCATE(prtint(nstop))

  DO ncount = 1,nstopsave
    prtint(nint-1+ncount) = tempreal(ncount)
  END DO   

  DEALLOCATE(tempreal) 

  IF (time > prtint(nint)) THEN
    WRITE(*,*) 
    WRITE(*,*) ' You have specified an output time < the restart time'
    WRITE(*,*) ' Set parameter "spatial_profile_at_time" greater than restart time'
    WRITE(*,*)
    WRITE(*,*) ' Restart time     : ', time
    WRITE(*,*) ' First output time: ', prtint(nint)
    WRITE(*,*)
!       ***** PETSc closeout**************
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
!     ****** PETSc closeout finished ********* 
    READ(*,*)
    STOP
  END IF
ELSE
  time = 0.0
  nn = 0
  nint = 1
END IF

IF (time > OutputTime(OutputTimeCounter) .AND. OutputTime(OutputTimeCounter) /= 0.0) THEN
  WRITE(*,*) 
  WRITE(*,*) ' You have specified an output time for a time series < the restart time'
  WRITE(*,*) ' Set parameter "time_series_output" greater than restart time'
  WRITE(*,*)
  WRITE(*,*) ' Restart time     : ', time
  WRITE(*,*) ' First time series output: ', OutputTime(OutputTimeCounter)
  WRITE(*,*)
!       ***** PETSc closeout**************
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
!    ****** PETSc closeout finished ********* 
  READ(*,*)
  STOP
END IF

IF (time+delt > OutputTime(OutputTimeCounter) .AND.                                     &
       OutputTime(OutputTimeCounter) /= time .AND. NumOutputTimes > OutputTimeCounter   &  
       .AND. OutputTime(OutputTimeCounter) /= 0.0) THEN
  delt = OutputTime(OutputTimeCounter) - time
  WRITE(*,*) ' Adjusting time step to match time series output'
  WRITE(*,5085) delt*OutputTimeScale
  WRITE(*,*)
END IF

RETURN

6001 CALL stringlen(restartfile,ls)
WRITE(*,*)
WRITE(*,*) ' Error in opening restart file: ', restartfile(1:ls)
WRITE(*,*)
!         ***** PETSc closeout**************
           IF (petscon) THEN
             call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
           END IF
           IF (CalculateFlow) THEN
             call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
           END IF
           IF (os3dpetsc) THEN
             call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
           END IF
           call PetscFinalize(ierr)
!       ****** PETSc closeout finished *********
READ(*,*)
STOP

5085 FORMAT(1X,' ---> New time step = ',1PE11.4)

RETURN
END SUBROUTINE restart
