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


SUBROUTINE read_pumplocations(nout,nx,ny,nz,nchem2)

USE crunchtype
USE CrunchFunctions
USE medium
USE params
USE flow
USE strings
USE concentration
USE io

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(IN)                                   :: nchem2

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lenformat
INTEGER(I4B)                                                :: lfile
CHARACTER (LEN=mls)                                         :: FileFormatType
CHARACTER (LEN=mls)                                         :: pumplocfile
LOGICAL(LGT)                                                :: ext
REAL(DP)                                                    :: nbpump
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: jx
INTEGER(I4B)                                                :: jy
INTEGER(I4B)                                                :: jz
real(DP),dimension(:), ALLOCATABLE     :: xdum1
real(DP),dimension(:), ALLOCATABLE     :: ydum1
real(DP),dimension(:), ALLOCATABLE     :: zdum1
INTEGER(I4B)                                                  :: FileNameLength
CHARACTER (LEN=mls)                                           :: FileTemp
INTEGER(I4B)                                                :: nco
INTEGER(I4B)                                                :: intbnd_tmp

pumplocfile = ' '

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
  
  IF (ssch == 'read_pumplocations') THEN
    
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!Looking for nb of pumps
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
          nbpump = DNUM(ssch)
        ELSE                !  An ascii string--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' Cant interpret string following read_pumplocations'
          WRITE(*,*) ' Looking for numerical value'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Look for geochemical conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        id = ids + ls
        CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF(ls /= 0) THEN
          lzs=ls
          CALL convan(ssch,lzs,res)
  
          DO nco = 1,nchem2
            IF (ssch == condlabel(nco)) THEN
              GO TO 50
            END IF
          END DO
          WRITE(*,*)
          WRITE(*,*) ' Geochemical condition in read_pumplocations not found'
          WRITE(*,*) ' Label = ',ssch
          WRITE(*,*)
          READ(*,*)
          STOP
          50         continue
          intbnd_tmp = nco
        ELSE         !  Blank string
          WRITE(*,*)
          WRITE(*,*) ' No geochemical condition in read_pumplocations provided'
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
      pumplocfile = ssch
      lfile = ls
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No file name provided in "read_pumplocations"'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Check for file format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      lenformat = ls
      CALL majuscules(ssch,lenformat)
      IF (ls /= 0) THEN
        IF (ssch == 'singlefile3d') THEN
          FileFormatType = 'SingleFile3D'
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' File format following "read_pumplocations" not recognized '
          WRITE(*,*) ' Only SingleFile3D accepted '
          WRITE(*,*)
          READ(*,*)
          STOP
        ENDIf
      ELSE    !! No file format provided, so assume default
        !!FileFormatType = 'SingleColumn'
        WRITE(*,*)
        WRITE(*,*) ' No file format following "read_pumplocations" '
        WRITE(*,*)
        READ(*,*)
      STOP
      ENDIF

    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No information provided following "read_pumplocations"  '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

!!  *****************************************************************    

!! READ PUMP LOCATIONS:

    IF (ALLOCATED(intbnd)) THEN
      DEALLOCATE(intbnd)
      ALLOCATE(intbnd(1,1:nx,1:ny,1:nz))
    ELSE
      ALLOCATE(intbnd(1,1:nx,1:ny,1:nz))
    END IF
  
    IF (ALLOCATED(npump)) THEN
      DEALLOCATE(npump)
      ALLOCATE(npump(1:nx,1:ny,1:nz))
    ELSE
      ALLOCATE(npump(1:nx,1:ny,1:nz))
    END IF
  
    IF (ALLOCATED(qg)) THEN
      DEALLOCATE(qg)
      ALLOCATE(qg(1,nx,ny,nz))
    ELSE
      ALLOCATE(qg(1,nx,ny,nz))
    END IF
  
    IF (ALLOCATED(xdum1)) THEN
      DEALLOCATE(xdum1)
      ALLOCATE(xdum1(1:int(nbpump)))
    ELSE
      ALLOCATE(xdum1(1:int(nbpump)))
    END IF
  
    IF (ALLOCATED(ydum1)) THEN
      DEALLOCATE(ydum1)
      ALLOCATE(ydum1(1:int(nbpump)))
    ELSE
      ALLOCATE(ydum1(1:int(nbpump)))
    END IF
  
    IF (ALLOCATED(zdum1)) THEN
      DEALLOCATE(zdum1)
      ALLOCATE(zdum1(1:int(nbpump)))
    ELSE
      ALLOCATE(zdum1(1:int(nbpump)))
    END IF

    IF (FileFormatType == 'SingleFile3D') THEN

      INQUIRE(FILE=pumplocfile,EXIST=ext)
      IF (.NOT. ext) THEN
        WRITE(*,*)
        WRITE(*,*) ' 3D pump file not found: ',pumplocfile(1:lfile)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    
      OPEN(UNIT=23,FILE=pumplocfile,STATUS='old',ERR=8005)
      FileTemp = pumplocfile
      CALL stringlen(FileTemp,FileNameLength)
      DO i = 1,int(nbpump)
      !!READ(23,*,END=1020) xdum,ydum,zdum
      READ(23,*,END=1020) xdum1(i),ydum1(i),zdum1(i)
    END DO
    
    jz=1
    
    DO jy = 1,ny
      DO jx = 1,nx
        npump(jx,jy,jz)=0
      END DO 
    END DO
    
        DO jy = 1,ny
          DO jx = 1,nx
            DO i = 1,int(nbpump)
            if (jx==int(xdum1(i)) .and. jy==int(ydum1(i))) then
              npump(jx,jy,jz)=1
              intbnd(1,jx,jy,jz)=intbnd_tmp
            end if
              
          END DO
          END DO
        END DO
    
      CLOSE(UNIT=23,STATUS='keep')
    END IF
!!  *****************************************************************    

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
END SUBROUTINE read_pumplocations
