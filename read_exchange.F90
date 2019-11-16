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

SUBROUTINE read_exchange(nout,ncomp,nexchange,data1,nexch_sec,nkin)
USE crunchtype
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
INTEGER(I4B), INTENT(IN)                                        :: ncomp
INTEGER(I4B), INTENT(OUT)                                       :: nexchange
CHARACTER (LEN=mls), INTENT(IN)                                 :: data1
INTEGER(I4B), INTENT(OUT)                                       :: nexch_sec
INTEGER(I4B), INTENT(IN)                                        :: nkin

!  Internal variables and arrays

CHARACTER (LEN=mls)                                             :: dummy1
CHARACTER (LEN=mls)                                             :: namtmp

INTEGER(I4B)                                                    :: id
INTEGER(I4B)                                                    :: iff
INTEGER(I4B)                                                    :: ids
INTEGER(I4B)                                                    :: ls
INTEGER(I4B)                                                    :: lzs
INTEGER(I4B)                                                    :: ltemp
INTEGER(I4B)                                                    :: lmin
INTEGER(I4B)                                                    :: k
INTEGER(I4B)                                                    :: lexc
INTEGER(I4B)                                                    :: n
INTEGER(I4B)                                                    :: ick
INTEGER(I4B)                                                    :: i
INTEGER(I4B)                                                    :: ix
INTEGER(I4B)                                                    :: ii

REAL(DP)                                                        :: exchange_tmp
REAL(DP)                                                        :: bfit_tmp


CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                  :: nam
REAL(DP), DIMENSION(:), ALLOCATABLE                             :: sto

ALLOCATE(nam(50))
ALLOCATE(sto(50))

REWIND nout

nexchange = 0
nexch_sec = 0
iexc = 1

100 READ(nout,'(a)',END=111) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'exchange') THEN
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
!            call convan(ssch,lzs,res)
      CALL stringtype(ssch,lzs,res)
      ltemp = lzs
      IF (res == 'a') THEN
        namtmp = ssch
        nexchange = nexchange + 1
        namexc(nexchange) = namtmp
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Looking for ASCII string for name of exchanger'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF (ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (ssch /= 'on') THEN
          WRITE(*,*) ' Name of exchanger should be followed by "on" or nothing at all'
          WRITE(*,*) ' Exchange species: ',namexc(nexchange)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        id = ids + ls
        CALL sschaine(zone,id,iff,ssch,ids,ls)
        lzs=ls
        lmin = ls
        CALL stringtype(ssch,lzs,res)
        IF (res /= 'a') THEN
          WRITE(*,*)
          WRITE(*,*) ' Looking for a mineral name, not a number'
          WRITE(*,*) ' In ion exchange block'
          WRITE(*,*) ' Exchange species: ',namexc(nexchange)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        namtmp = ssch
        DO k = 1,nkin
          IF (namtmp == umin(k)) THEN
            kexch(nexchange) = k
            GO TO 50
          END IF
        END DO
        WRITE(*,*)
        WRITE(*,*) ' Mineral substrate listed in ion exchange block'
        WRITE(*,*) '     not found in minerals list'
        WRITE(*,*) ' Exchange species: ',namexc(nexchange)
        WRITE(*,*) ' Looking for mineral: ',namtmp(1:lmin)
        WRITE(*,*)
        READ(*,*)
        STOP
        
        50 iexchange(nexchange) = 1     ! Exchange on a specific mineral
        GO TO 100
        
      ELSE
        iexchange(nexchange) = 0     ! Exchange on bulk sediment
      END IF
      
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Exchange keyword specified, but no exchanger name given'
      WRITE(*,*) ' In ion exchange block'
      WRITE(*,*) ' Exchanger name must be specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
  END IF
END IF
GO TO 100

!  Now, look for the activity convention to be used

111 REWIND nout

200 READ(nout,'(a)',END=222) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'convention') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      lexc = lzs
      namtmp = ssch
      IF (namtmp == 'gaines-thomas') THEN
        iexc = 1
      ELSE IF (namtmp == 'vanselow') THEN
        iexc = 2
      ELSE IF (namtmp == 'gapon') THEN
        iexc = 3
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Dont understand exchanger activity convention'
        WRITE(*,*) namtmp(1:lexc)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No convention for exchanger activity given'
      WRITE(*,*) ' Assuming Gaines-Thomas convention'
      iexc = 1
      WRITE(*,*)
    END IF
  END IF
END IF
GO TO 200


222 OPEN(UNIT=18,FILE=data1,STATUS='old')

!  Find the beginning of the exchange section

300 READ(18,'(a)') dummy1
IF (dummy1 == 'Begin exchange') THEN
  GO TO 400
ELSE
  GO TO 300
END IF

400 READ(18,'(a)',END=444) dummy1
IF (dummy1 == 'End of exchange') THEN
  CLOSE(UNIT=18,STATUS='keep')
  DEALLOCATE(nam)
  DEALLOCATE(sto)
  RETURN
ELSE
  BACKSPACE 18
  READ(18,*,ERR=6003) nam(1),n,(sto(i+1),nam(i+1), i = 1, n), exchange_tmp,bfit_tmp
  
  DO ick = 1,n
    
    DO i = 1,ncomp
      IF (nam(ick+1) == ulab(i)) THEN
        GO TO 555
      END IF
    END DO
    DO ix = 1,nexchange
      IF (nam(ick+1) == namexc(ix)) THEN
        GO TO 555
      END IF
    END DO
!  Species not found in primary species or exchanger lists, skip reaction
    GO TO 400
    555     CONTINUE   ! Goto to next species in reaction
    
  END DO
  
  nexch_sec = nexch_sec + 1
  nam_exchsec(nexch_sec) = nam(1)
  DO i = 1,ncomp
    DO ick = 1,n
      IF (nam(ick+1) == ulab(i)) THEN
        muexc(nexch_sec,i) = sto(ick+1)
      END IF
    END DO
  END DO
  DO ix = 1,nexchange
    ii = ix + ncomp
    DO ick=1,n
      IF (nam(ick+1) == namexc(ix)) THEN
        ixlink(nexch_sec) = ix
        muexc(nexch_sec,ii) = sto(ick+1)
      END IF
    END DO
  END DO
  DO i = 1,ncomp
    DO ick=1,n
      IF (nam(ick+1) == ulab(i)) THEN
        nclink(nexch_sec) = i
      END IF
    END DO
  END DO
  keqexc(nexch_sec) = exchange_tmp*clg
  bfit(nexch_sec) = bfit_tmp*clg
  
END IF

GO TO 400

444 CLOSE(UNIT=18,STATUS='keep')

DEALLOCATE(nam)
DEALLOCATE(sto)

RETURN

6003 WRITE(*,*) ' Error in reading ion exchange reactions'
READ(*,*)
STOP

END SUBROUTINE read_exchange
