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


SUBROUTINE read_par(nout,lchar,parchar,parfind,realjunk,section)
USE crunchtype
USE CrunchFunctions
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
INTEGER(I4B), INTENT(OUT)                                       :: lchar

CHARACTER (LEN=mls), INTENT(IN)                                 :: section
CHARACTER (LEN=mls), INTENT(IN)                                 :: parchar
CHARACTER (LEN=mls), INTENT(IN OUT)                             :: parfind
REAL(DP), INTENT(IN OUT)                                        :: realjunk

!  Internal variables and arrays

INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: iff
INTEGER(I4B)                                                    :: ids
INTEGER(I4B)                                                    :: ls
INTEGER(I4B)                                                    :: lzs
INTEGER(I4B)                                                    :: lensection
INTEGER(I4B)                                                    :: lenparchar
INTEGER(I4B)                                                    :: nlen1

CHARACTER (LEN=mls)                                             :: LowerCase
CHARACTER (LEN=mls)                                             :: LowerCase_ssch


REWIND nout

CALL stringlen(section,lensection)
CALL stringlen(parchar,lenparchar)

LowerCase = parchar
call majuscules(LowerCase,lenparchar)

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
END IF

LowerCase_ssch = ssch
call majuscules(LowerCase_ssch,lzs)
IF (LowerCase_ssch == LowerCase) THEN
!!IF (ssch == parchar .OR. ssch == LowerCase) THEN
  parfind = parchar
  lchar = ls
  GO TO 200
ELSE
  GO TO 100
END IF
300 RETURN

200 id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (res /= 'n') THEN
    WRITE(*,*)
    WRITE(*,*) ' Parameter should be followed by a numeric value'
    WRITE(*,*) '   In section: ',section(1:lensection)
    WRITE(*,*) '   Following parameter: ',parchar(1:lenparchar)
    WRITE(*,*) ' ABORTING RUN'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  realjunk = DNUM(ssch)
  RETURN
ELSE
  realjunk = 0.00
END IF

RETURN
END SUBROUTINE read_par
