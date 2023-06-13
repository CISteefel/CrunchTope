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


SUBROUTINE read_pumptimeseries(nout,nx,ny,nz)

USE crunchtype
USE CrunchFunctions
USE medium
USE params
USE flow
USE strings
USE io

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lfile
LOGICAL(LGT)                                                :: ext
CHARACTER (LEN=mls)                                           :: FileTemp
INTEGER(I4B)                                                   :: tp
INTEGER(I4B)                                                  :: FileNameLength

CHARACTER (LEN=mls)                                         :: pumptimeseriesfile
REAL(DP)                                                    :: tslength


pumptimeseriesfile = ' '

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
  
  IF (ssch == 'read_pumptimeseries') THEN
    pumptimeseries = .true.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Looking for time series length
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        tslength = DNUM(ssch)
      ELSE                
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "read_pumptimeseries"'
        WRITE(*,*) ' Looking for length time_series' 
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Looking for file name
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL stringtype(ssch,lzs,res)
      pumptimeseriesfile = ssch
      lfile = ls
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No file name provided in "read_pumptimeseries"'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
      

!!  *****************************************************************

    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No information provided following "read_pumptimeseries" '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  IF (pumptimeseries) THEN

    WRITE(iunit2,*)
    WRITE(iunit2,*) '  Reading pumptimeseries from file: ',pumptimeseriesfile(1:lfile)
    WRITE(iunit2,*)

    IF (ALLOCATED(tpump)) THEN
      DEALLOCATE(tpump)
      ALLOCATE(tpump(1:int(tslength)))
    ELSE
      ALLOCATE(tpump(1:int(tslength)))
    END IF
  
    IF (ALLOCATED(qgt)) THEN
      DEALLOCATE(qgt)
      ALLOCATE(qgt(1:int(tslength)))
    ELSE
      ALLOCATE(qgt(1:int(tslength)))
    END IF


    INQUIRE(FILE=pumptimeseriesfile,EXIST=ext)
    IF (.NOT. ext) THEN
      WRITE(*,*)
      WRITE(*,*) ' 3D pump file not found: ',pumptimeseriesfile(1:lfile)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

    OPEN(UNIT=23,FILE=pumptimeseriesfile,STATUS='old',ERR=8005)
    FileTemp = pumptimeseriesfile
    CALL stringlen(FileTemp,FileNameLength)
          DO tp = 1,int(tslength)
            READ(23,*,END=1020) tpump(tp),qgt(tp)
            qgt(tp)=((qgt(tp))/1000.0d0)*dxx(nx)*dzz(1,1,1) !! Converting from mm/year to m3/year
          END DO
            !qg(1,nx,ny,nz)=qgt(1)

    CLOSE(UNIT=23,STATUS='keep')
  
  ENDIF



  ELSE
    GO TO 10
  END IF
  GO TO 10
  
ELSE         ! No string found
  GO TO 10
END IF

GO TO 10
!!!!!!!!!!!!!!!!!!!!

RETURN

1020  WRITE(*,*) ' End of file during read'
WRITE(*,*) ' Trying to read the file: ', FileTemp(1:FileNameLength)
READ(*,*)
STOP

8005   WRITE(*,*) ' Error opening pumptimeseries file'
        READ(*,*)
STOP

1000 RETURN

END SUBROUTINE read_pumptimeseries
