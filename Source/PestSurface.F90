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

    
SUBROUTINE PestSurface(nout,ncomp,nchem)
USE CrunchType
USE CrunchFunctions
USE params
USE RunTime, ONLY: PestSurfaceOutputFile
USE strings
USE concentration

IMPLICIT NONE

INTERFACE
  SUBROUTINE PrimarySpeciesCheck(ncomp,npar,DummyStringArray)
  USE CrunchType
  USE concentration
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                        :: ncomp
  INTEGER(I4B), INTENT(IN)                                        :: npar
  CHARACTER (LEN=mls), DIMENSION(:)                               :: DummyStringArray
  END SUBROUTINE PrimarySpeciesCheck
END INTERFACE

!! External arrays and variables

INTEGER(I4B), INTENT(IN)                                      :: nout
INTEGER(I4B), INTENT(IN)                                      :: ncomp
INTEGER(I4B), INTENT(IN)                                      :: nchem

!! Internal arrays and variables
INTEGER(I4B)                                                  :: id
INTEGER(I4B)                                                  :: iff
INTEGER(I4B)                                                  :: ids
INTEGER(I4B)                                                  :: ls
INTEGER(I4B)                                                  :: lzs
INTEGER(I4B)                                                  :: lchar
INTEGER(I4B)                                                  :: l_string
INTEGER(I4B)                                                  :: lenarray
INTEGER(I4B)                                                  :: npar
INTEGER(I4B)                                                  :: nco
INTEGER(I4B)                                                  :: nl

CHARACTER (LEN=mls)                                           :: dumstring
CHARACTER (LEN=mls)                                           :: dumstring2
CHARACTER (LEN=mls)                                           :: dumstring3
CHARACTER (LEN=mls)                                           :: OpenBracket
CHARACTER (LEN=mls)                                           :: CloseBracket
CHARACTER (LEN=mls)                                           :: CharLengthString
CHARACTER (LEN=mls)                                           :: parchar
CHARACTER (LEN=mls)                                           :: section
CHARACTER (LEN=mls)                                           :: parfind

LOGICAL(LGT)                                                  :: continuation
LOGICAL(LGT)                                                  :: ext
LOGICAL(LGT)                                                  :: FileOpen 

continuation = .FALSE.
UnitPestSurface = 52

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
lzs=ls
CALL convan(ssch,lzs,res)
IF (ssch == 'surface') THEN
  parfind = parchar
  lchar = ls
  GO TO 200
ELSE
  GO TO 100
END IF
300 RETURN

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)                   !!  Now, look for units for surface complexation
l_string = ls
IF(ls /= 0) THEN
  npar = 0
  lzs=ls  
  IF (ssch == 'mol/g' .OR. ssch == 'mole/g' .OR. ssch == 'moles/g') THEN
    PestSurfaceUnits = 'mol/g'
  ELSE IF (ssch == 'millimol/g' .OR. ssch == 'millimol/g' .OR. ssch == 'millimoles/g' .OR.  &
                  ssch == 'mmol/g' .OR. ssch == 'mmole/g' .OR. ssch == 'mmoles/g') THEN
    PestSurfaceUnits = 'mmol/g'
  ELSE IF (ssch == 'micromol/g' .OR. ssch == 'micromole/g' .OR. ssch == 'micromoles/g' .OR. &
             ssch == 'umol/g' .OR. ssch == 'umole/g' .OR. ssch == 'umoles/g') THEN
    PestSurfaceUnits = 'umol/g'
  ELSE IF (ssch == 'log' .OR. ssch == 'logs' .OR. ssch == 'logarithm' .OR. ssch == 'logarithms') THEN
    PestSurfaceUnits = 'log'
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Pest Surface units not recognized in PEST block'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  500 CONTINUE
  IF (continuation) THEN
    id = 1
    continuation = .FALSE.
  ELSE
    id = ids + ls
  END IF

  CALL sschaine(zone,id,iff,ssch,ids,ls)

  IF (ssch == '&') THEN
    READ(nout,'(a)') zone
    continuation = .TRUE.
    GOTO 500
  END IF

  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
    npar = npar + 1
    stringarray(npar) = ssch
    GO TO 500
  ELSE
    CONTINUE
  END IF

  lenarray = npar

  IF (npar == 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' No species specified after "surface" keyword in PEST block'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  NPestSurface = npar
  IF (ALLOCATED(PestSurfaceList)) THEN
    DEALLOCATE(PestSurfaceList)
    ALLOCATE(PestSurfaceList(NPestSurface))
  ELSE
    ALLOCATE(PestSurfaceList(NPestSurface))
  END IF
  PestSurfaceList(1:npar) = StringArray(1:NPestSurface)
    
  CALL PrimarySpeciesCheck(ncomp,NPestSurface,PestSurfaceList)

501 FORMAT('pif @')
!!__GNUC__ 502 FORMAT(a2,a2,a<ls>,a<lchar>,a1,a4)

  INQUIRE(FILE=PestSurfaceOutputFile,OPENED=FileOpen)
  IF (FileOpen) THEN
    CLOSE(UNIT=UnitPestSurface,STATUS='delete')
    OPEN(UNIT=UnitPestSurface,FILE=PestSurfaceOutputFile,STATUS='new')
  ELSE
    OPEN(UNIT=UnitPestSurface,FILE=PestSurfaceOutputFile,STATUS='unknown')
  END IF

ELSE
  WRITE(*,*)
  WRITE(*,*) ' No trailing string found after'
  WRITE(*,*) parchar
  WRITE(*,*) ' In section ',section
  WRITE(*,*)
  parfind = ' '
END IF

RETURN
END SUBROUTINE PestSurface
