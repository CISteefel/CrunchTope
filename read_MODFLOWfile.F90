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


SUBROUTINE read_MODFLOWfile(nout,lfile,mxwell,mxrivr,mxdrn)
USE crunchtype
USE params
USE strings
USE runtime

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout

INTEGER(I4B), INTENT(OUT)                                   :: lfile
INTEGER(I4B), INTENT(OUT)                                   :: mxwell
INTEGER(I4B), INTENT(OUT)                                   :: mxrivr
INTEGER(I4B), INTENT(OUT)                                   :: mxdrn

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs

CHARACTER (LEN=mls)                                         :: file

LOGICAL(LGT)                                                :: ext

MODFLOWfile = ' '
modflow = .FALSE.

REWIND nout

10  READ(nout,'(a)',END=1000) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res) 
  IF (ssch == 'read_modflow') THEN

!  Switch on MODFLOW option only if file found

    modflow = .TRUE.
    
!  Looking for MODFLOW file name (expecting prefix only)
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      MODFLOWfile = ssch
      lfile = ls
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No file name following "read_modflow" in MODFLOW section '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
  ELSE
    GO TO 10
  END IF
  
ELSE         ! No string found
  GO TO 10
END IF

IF (MODFLOWfile == ' ') THEN
  RETURN
ELSE
  file(1:lfile)=MODFLOWfile(1:lfile)
  file(lfile+1:lfile+1+3)=".hff"
  INQUIRE(FILE=file,EXIST=ext)
  IF (.NOT. ext) THEN
    WRITE(*,*) 
    WRITE(*,*) ' MODFLOW "*.hff" file not found: ', file(1:lfile+4)
    WRITE(*,*) 
    READ(*,*)
    STOP
  ELSE
    OPEN(1,file=file,form='unformatted',status='old')
  END IF

  file = ' '
  file(1:lfile)=MODFLOWfile(1:lfile)
  file(lfile+1:lfile+1+3)=".bas"
  INQUIRE(FILE=file,EXIST=ext)
  IF (.NOT. ext) THEN
    WRITE(*,*) 
    WRITE(*,*) ' MODFLOW "*.bas" file not found: ', file(1:lfile+4)
    WRITE(*,*)
    READ(*,*) 
    STOP
  ELSE
    OPEN(52,file=file,status='old')
  END IF

  file = ' '
  file(1:lfile)=MODFLOWfile(1:lfile)
  file(lfile+1:lfile+1+3)=".wel"
  INQUIRE(FILE=file,EXIST=ext)
  IF (.NOT. ext) THEN
    mxwell = 0
  ELSE
    OPEN(53,file=file,status='old')
    READ(53,*) mxwell
    CLOSE(53) 
  END IF

  file = ' '
  file(1:lfile)=MODFLOWfile(1:lfile)
  file(lfile+1:lfile+1+3)=".drn"
  INQUIRE(FILE=file,EXIST=ext)
  IF (.NOT. ext) THEN
    mxdrn = 0
  ELSE
    OPEN(53,file=file,status='old')
    READ(53,*) mxdrn
    CLOSE(53) 
  END IF

  file = ' '
  file(1:lfile)=MODFLOWfile(1:lfile)
  file(lfile+1:lfile+1+3)=".riv"
  INQUIRE(FILE=file,EXIST=ext)
  IF (.NOT. ext) THEN
    mxrivr = 0
  ELSE
    OPEN(53,file=file,status='old')
    READ(53,*) mxrivr
    CLOSE(53) 
  END IF

!!  file = ' '
!!  file(1:lfile)=MODFLOWfile(1:lfile)
!!  file(lfile+1:lfile+1+6)=".hffout"
!!  INQUIRE(FILE=file,EXIST=ext)
!!  IF (.NOT. ext) THEN
!!    WRITE(*,*) 
!!    WRITE(*,*) ' MODFLOW "*.hffout" file not found: ', file(1:lfile+7)
!!    WRITE(*,*) 
!!    STOP
!!  ELSE
!!    OPEN(3,file=file)
!!  END IF
END IF

1000 RETURN
END SUBROUTINE read_MODFLOWfile

