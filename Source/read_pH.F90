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


SUBROUTINE read_ph(nout,ph,guessph,i,isolution,constraint,nrct,phfound)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
REAL(DP), DIMENSION(:), INTENT(IN OUT)                      :: ph
REAL(DP), DIMENSION(:), INTENT(IN OUT)                      :: guessph
INTEGER(I4B), INTENT(IN)                                    :: i
INTEGER(I4B), INTENT(IN)                                    :: isolution
CHARACTER (LEN=mls), DIMENSION(:,:), INTENT(IN OUT)         :: constraint
INTEGER(I4B), INTENT(IN)                                    :: nrct
LOGICAL(LGT), INTENT(IN OUT)                                :: phfound

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: tempstring

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: k


REWIND nout

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
!!!    WRITE(*,*) ' Geochemical input should start with a string'
!!!    WRITE(*,*) ssch
!!!    WRITE(*,5050) isolution
!!!    WRITE(*,*)
!!!    READ(*,*)
!!!    STOP
!!!  END IF
END IF
IF (ssch == 'pH' .OR. ssch == 'ph' .OR. ssch == 'PH') THEN
  phfound = .true.
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
    
    itype(i,isolution) = 7        !!! Identifies it as an "activity"
    ph(isolution) = DNUM(ssch)
    
  ELSE
    
!  An ascii string, so check for charge balance, then a mineral name
    IF (ssch == 'charge' .OR. ssch == 'Charge' .OR. ssch == 'CHARGE') THEN
      itype(i,isolution) = 2
      GO TO 200
    END IF
    DO k = 1,nrct
      tempstring = umin(k)
      IF (ssch == tempstring) THEN
        constraint(i,isolution) = umin(k)
        itype(i,isolution) = 3
        GO TO 200
      END IF
    END DO
    WRITE(*,*)
    WRITE(*,*) ' Mineral constraint for pH not found'
    WRITE(*,5050) isolution
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF
!  Now, check to see if there is an optional guess for pH
200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res == 'n') THEN
!  Read the guess for the pH value
    guessph(isolution) = DNUM(ssch)
  ELSE
!  Anything else, ignore it
    CONTINUE
  END IF
  
END IF

5050 FORMAT(1X,'Condition number ',i2,' in input file')

300  RETURN
END SUBROUTINE read_ph
