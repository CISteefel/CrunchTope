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


SUBROUTINE read_concentration(nout,i,isolution,constraint,  &
    ncomp,nspec,nrct,ngas,speciesfound)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                     :: nout
INTEGER(I4B), INTENT(IN)                                     :: i
INTEGER(I4B), INTENT(IN)                                     :: isolution
CHARACTER (LEN=mls), DIMENSION(:,:), INTENT(IN OUT)          :: constraint
INTEGER(I4B), INTENT(IN)                                     :: ncomp
INTEGER(I4B), INTENT(IN)                                     :: nspec
INTEGER(I4B), INTENT(IN)                                     :: nrct
INTEGER(I4B), INTENT(IN)                                     :: ngas
LOGICAL(LGT), INTENT(IN OUT)                                 :: speciesfound

!  Internal variables and arrays

CHARACTER (LEN=mls)                                          :: tempstring
LOGICAL(LGT)                                                 :: constraingas
LOGICAL(LGT)                                                 :: constrainmin

INTEGER(I4B)                                                 :: id
INTEGER(I4B)                                                 :: iff
INTEGER(I4B)                                                 :: ids
INTEGER(I4B)                                                 :: ls
INTEGER(I4B)                                                 :: lzs
INTEGER(I4B)                                                 :: nlen1
INTEGER(I4B)                                                 :: k
INTEGER(I4B)                                                 :: ids_save
INTEGER(I4B)                                                 :: ll
INTEGER(I4B)                                                 :: idsave
INTEGER(I4B)                                                 :: ls_save
INTEGER(I4B)                                                 :: lcond

CHARACTER (LEN=mls)                                          :: dumstring
CHARACTER (LEN=mls)                                          :: stringspecies

speciesfound = .false.
REWIND nout
constraingas = .false.
constrainmin = .false.

dumstring = condlabel(isolution)
CALL stringlen(dumstring,lcond)

100 READ(nout,'(a)',END=300) zone
nlen1 = LEN(zone)
!      call majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
!        call convan(ssch,lzs,res)
  CALL stringtype(ssch,lzs,res)
!!!  IF (res /= 'a') THEN
!!!    WRITE(*,*)
!!!    WRITE(*,*) ' Geochemical input should start with an ASCII string'
!!!    WRITE(*,*) '   In condition: ', dumstring(1:lcond)
!!!    WRITE(*,*) '   String found', ssch(1:30)
!!!    WRITE(*,*) '       ABORTING RUN  '
!!!    WRITE(*,*)
!!!    READ(*,*)
!!!    STOP
!!!    WRITE(*,*)
!!!    READ(*,*)
!!!    STOP
!!!  END IF
END IF
IF (ssch == ulab(i)) THEN
  speciesfound = .true.     ! primary species found
ELSE
  GO TO 100
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
!        call convan(ssch,lzs,res)
  CALL stringtype(ssch,lzs,res)
  IF (res == 'n') THEN
!  Read the concentration
    ctot(i,isolution) = DNUM(ssch)
    IF (ctot(i,isolution) == 0.0) THEN
      ctot(i,isolution) = 1.e-30
    END IF
    itype(i,isolution) = 1
    
!  Check to see if total concentration should include surface species (exchange and surface complexes)
    
    idsave = id
    ls_save = ls
    ids_save = ids
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'equilibrate_surface') THEN
        equilibrate(i,isolution) = .true.
      ELSE
        id = idsave
        ls = ls_save
        ids = ids_save
      END IF
    ELSE
      id = idsave
      ls = ls_save
      ids = ids_save
    END IF
    
  ELSE
!  An ascii string--
!  First, check to see if it is "charge"
    IF (ssch == 'charge' .OR. ssch == 'Charge' .OR. ssch == 'CHARGE') THEN
      itype(i,isolution) = 2
      GO TO 200
    END IF
    
!    then, check to see if it is a mineral
    DO k = 1,nrct
      tempstring = umin(k)
      IF (ssch == tempstring) THEN
        constraint(i,isolution) = umin(k)
        constrainmin = .true.
        GO TO 200
      END IF
    END DO
!  next, check to see if it is a gas
    DO ll = 1,ngas
      tempstring = namg(ll)
      IF (ssch == tempstring) THEN
        constraint(i,isolution) = namg(ll)
        constraingas = .true.
        GO TO 200
      END IF
    END DO

    WRITE(*,*)
    WRITE(*,*) ' Mineral or gas constraint not found'
    stringspecies = ulab(i)
    CALL stringlen(stringspecies,ls)
    WRITE(*,*) '   For species:  ',stringspecies(1:ls)
    WRITE(*,*) '   In condition: ', dumstring(1:lcond)
    READ(*,*)
    STOP
  END IF
ELSE
  WRITE(*,*)
  WRITE(*,*)   ' No information given following primary species'
  stringspecies = ulab(i)
  CALL stringlen(stringspecies,ls)
  WRITE(*,*) '   For species:  ',stringspecies(1:ls)
  WRITE(*,*)   '   In condition: ', dumstring(1:lcond)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)

IF (constraingas) THEN    ! Using gas constraint
  itype(i,isolution) = 4
  IF(ls /= 0) THEN
    lzs=ls
!          call convan(ssch,lzs,res)
    CALL stringtype(ssch,lzs,res)
    IF (res == 'n') THEN
      gaspp(i,isolution) = DNUM(ssch)
      IF (gaspp(i,isolution) == 0.0) THEN
        gaspp(i,isolution) = 1.e-30
      END IF

    ELSE    ! Gas constraint, but trailing string not a number
      WRITE(*,*)
      WRITE(*,*) ' Input following a gas constraint should be the partial pressure'
      stringspecies = ulab(i)
      CALL stringlen(stringspecies,ls)
      WRITE(*,*) '   For species:  ',stringspecies(1:ls)
      WRITE(*,*) '   In condition: ', dumstring(1:lcond)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE      !  Gas constraint, but no trailing string found
    WRITE(*,*)
    WRITE(*,*) ' Gas constraint, no partial pressure given'
    stringspecies = ulab(i)
    CALL stringlen(stringspecies,ls)
    WRITE(*,*) '   For species:  ',stringspecies(1:ls)
    WRITE(*,*) '   In condition: ', dumstring(1:lcond)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
!  Now, check again to see if there is a guess following the specification
!    of the gas partial pressure
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
!          call convan(ssch,lzs,res)
    CALL stringtype(ssch,lzs,res)
    IF (res == 'n') THEN
      guess(i,isolution) = DNUM(ssch)
      RETURN
    ELSE          ! ASCII string when expecting a number--ignore it
      RETURN
    END IF
  ELSE        !  No guess provided
    RETURN
  END IF
END IF

IF (constrainmin) THEN    ! Using mineral constraint
  itype(i,isolution) = 3
  IF(ls /= 0) THEN
    lzs=ls
!          call convan(ssch,lzs,res)
    CALL stringtype(ssch,lzs,res)
    IF (res == 'n') THEN
      guess(i,isolution) = DNUM(ssch)
      RETURN
    ELSE    ! Mineral constraint, but trailing string not a number
      WRITE(*,*)
      WRITE(*,*) ' String following a mineral constraint should be a guess'
      WRITE(*,*) '   Ignoring trailing string'
      stringspecies = ulab(i)
      CALL stringlen(stringspecies,ls)
      WRITE(*,*) '   For species:  ',stringspecies(1:ls)
      WRITE(*,*) '   In condition: ', dumstring(1:lcond)
      WRITE(*,*)
      RETURN
    END IF
  ELSE      !  Mineral constraint, no guess provided
    RETURN
  END IF
END IF

! Case where a numerical value (concentration or activity) provided following
!   the primary species name--check for constraint specification

IF(ls /= 0) THEN
  lzs=ls
!        call convan(ssch,lzs,res)
  CALL stringtype(ssch,lzs,res)
  IF (res == 'a') THEN       ! ASCII string
    IF (ssch == 'total') THEN
      itype(i,isolution) = 1
    ELSE IF (ssch == 'species') THEN
      itype(i,isolution) = 8
    ELSE IF (ssch == 'activity') THEN
      itype(i,isolution) = 7
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Dont understand trailing string'
      WRITE(*,*) '   Ignoring trailing string'
      stringspecies = ulab(i)
      CALL stringlen(stringspecies,ls)
      WRITE(*,*) '   For species:  ',stringspecies(1:ls)
      WRITE(*,*) '   In condition: ', dumstring(1:lcond)
    END IF
    
!  Now, check for a guess
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
!            call convan(ssch,lzs,res)
      CALL stringtype(ssch,lzs,res)
      IF (res == 'n') THEN
        guess(i,isolution) = DNUM(ssch)
      ELSE          ! ASCII string when expecting a number--ignore it
        CONTINUE
      END IF
    ELSE        !  No guess provided
      RETURN
    END IF
  ELSE
    itype(i,isolution) = 1     ! No trailing ASCII string, so assume
!                                      default--total concentration
    IF(ls /= 0) THEN
      guess(i,isolution) = DNUM(ssch)
    ELSE        !  No guess provided
      RETURN
    END IF
  END IF
END IF

5050 FORMAT(1X,'Condition number ',i2,' in input file')

300  RETURN
END SUBROUTINE read_concentration
