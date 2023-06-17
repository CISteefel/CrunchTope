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


SUBROUTINE read_tempts(nout,tslength)
USE crunchtype
USE CrunchFunctions
USE params
USE flow
USE temperature
USE strings
USE medium

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(OUT)                                   :: tslength

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
REAL(DP)                                                    :: temp ! region ID allocated to each time series
REAL(DP)                                                    :: regid
INTEGER(I4B)                                   :: lfile
CHARACTER (LEN=mls), ALLOCATABLE, DIMENSION(:)                           :: temptsfile
LOGICAL(LGT)                                               :: ext
INTEGER(I4B)                                                  :: FileNameLength
CHARACTER (LEN=mls)                                           :: FileTemp
INTEGER(I4B)                                                   :: tp
INTEGER(I4B)                                                     :: ierr
REAL(DP)                                                      :: t_dum
REAL(DP)                                                      :: temp_dum
INTEGER(I4B)                                                     :: i
CHARACTER (LEN=mls)                                           :: stringdum
REAL(DP), DIMENSION(:), ALLOCATABLE          :: check1 ! temperature of the time series
REAL(DP), DIMENSION(:,:), ALLOCATABLE          :: check2 ! time of the time series

!temptsfile = ' '

IF (ALLOCATED(reg_temp_fix)) THEN
  DEALLOCATE(reg_temp_fix)
END IF
ALLOCATE(reg_temp_fix(nbreg))

IF (ALLOCATED(temp_fix)) THEN
  DEALLOCATE(temp_fix)
END IF
ALLOCATE(temp_fix(nbreg))

IF (ALLOCATED(reg_temp_ts)) THEN
  DEALLOCATE(reg_temp_ts)
END IF
ALLOCATE(reg_temp_ts(nbreg))

IF (ALLOCATED(temptsfile)) THEN
  DEALLOCATE(temptsfile)
END IF
ALLOCATE(temptsfile(nbreg))


REWIND nout

10  READ(nout,'(a)',END=1000) zone
nlen1 = LEN(zone)
!!CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
  IF (ssch == 'tempts') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      !!  *****************************************************************     
      !!check for zone ID
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        regid = DNUM(ssch)
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' A region ID must be provided for the time series '
        WRITE(*,*)
        READ(*,*)
      STOP
      ENDIF


      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
          lzs=ls
          !!  *****************************************************************     
          !!check for numerical values -> temperature fixed
          CALL convan(ssch,lzs,res)
          IF (res == 'n') THEN
            temp = DNUM(ssch)
            nb_temp_fix = nb_temp_fix + 1
            reg_temp_fix(nb_temp_fix) = regid
            temp_fix(nb_temp_fix) = temp
            
              !STOP
              !!  *****************************************************************
              !!check for the name of timeseries file
          ELSEIF (res /= 'n') THEN
              CALL stringtype(ssch,lzs,res)
              lfile = ls
              RunTempts=.true.
              nb_temp_ts = nb_temp_ts + 1
              temptsfile(nb_temp_ts) = ssch
              reg_temp_ts(nb_temp_ts)=regid !store the region ID where the time series is applied
              
              id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
        tslength = int(DNUM(ssch))

        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Must provide a numerical value for length time series temperature'
          WRITE(*,*)
          READ(*,*)
          STOP  
        ENDIF
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Must provide a length for time series temperature '
        WRITE(*,*)
        READ(*,*)
        STOP  
        ENDIF
              
          ENDIF
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Must provide either value or time series for temperature '
        WRITE(*,*)
        READ(*,*)
        STOP
      ENDIF
      !!  *****************************************************************

    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No info provided on temperature time series '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
  ELSE
    GO TO 10
  END IF
  GO TO 10
  
ELSE         ! No string found
  GO TO 10
END IF

1000 IF (nb_temp_ts>0) THEN

  IF (ALLOCATED(t_temp_ts)) THEN
    DEALLOCATE(t_temp_ts)
  END IF
  ALLOCATE(t_temp_ts(tslength))

  IF (ALLOCATED(temp_ts)) THEN
    DEALLOCATE(temp_ts)
  END IF
  ALLOCATE(temp_ts(nb_temp_ts,tslength))


 
  DO i=1,nb_temp_ts
    stringdum=temptsfile(i)
  INQUIRE(FILE=stringdum,EXIST=ext)
        IF (.NOT. ext) THEN
      WRITE(*,*)
      WRITE(*,*) ' time series file not found: ',stringdum(1:len(stringdum))
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

    OPEN(UNIT=23,FILE=stringdum,STATUS='old',ERR=8005)
    FileTemp = stringdum
    CALL stringlen(FileTemp,FileNameLength)

    DO tp = 1,tslength
      READ(23,*,iostat=IERR) t_temp_ts(tp),temp_ts(i,tp)
    END DO

  END DO
  check1=t_temp_ts
  check2=temp_ts
 !STOP
ENDIF
RETURN

8005   WRITE(*,*) ' Error opening timeseries file', temptsfile(1:lfile)
        READ(*,*)
        STOP 

END SUBROUTINE read_tempts
