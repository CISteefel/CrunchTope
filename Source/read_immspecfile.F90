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


SUBROUTINE read_immspecfile(nout,nx,ny,nz,readimmspec,immspec_index,immspec_id,immspec_name,immspec_name_length,concfile,lfile_conc,FileFormatType_conc)
USE crunchtype
USE params
USE flow
USE strings
USE concentration
USE transport

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                                               :: nout
INTEGER(I4B), INTENT(IN)                                                               :: nx
INTEGER(I4B), INTENT(IN)                                                               :: ny
INTEGER(I4B), INTENT(IN)                                                               :: nz
LOGICAL(LGT), INTENT(IN OUT)                                                           :: readimmspec
INTEGER(I4B), INTENT(IN OUT)                                                           :: immspec_index
INTEGER(I4B),  INTENT(IN OUT)                                   :: immspec_id(50)
INTEGER(I4B),  INTENT(IN OUT)                                   :: lfile_conc(50)
CHARACTER (LEN=mls),  INTENT(IN OUT)                            :: immspec_name(50)
INTEGER(I4B),  INTENT(IN OUT)                             :: immspec_name_length(50)
CHARACTER (LEN=mls),  INTENT(IN OUT)                            :: concfile(50)
CHARACTER (LEN=mls),  INTENT(IN OUT)                           :: FileFormatType_conc(50)

!  Internal variables and arrays

INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: immspecname_lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lenformat
CHARACTER (LEN=mls)                                         :: immspecname
INTEGER(I4B)                                                :: dummyboolean
INTEGER(I4B)                                                :: tracer
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                     :: immspec_id_dum
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                     :: lfile_conc_dum
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                     :: immspecname_dum
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                     :: immspecname_lzs_dum
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                     :: concfile_dum
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                     :: FileFormatType_conc_dum
INTEGER(I4B)                                                :: dummy

nxyz = nx*ny*nz
readimmspec = .false.
!volfracfile = ' '
!bsafile = ' '
!mineral_index = 0
tracer = 0
dummy = size(ulab)

IF (ALLOCATED(immspec_id_dum)) THEN
  DEALLOCATE(immspec_id_dum)
ENDIF
ALLOCATE(immspec_id_dum(50))

IF (ALLOCATED(immspecname_dum)) THEN
  DEALLOCATE(immspecname_dum)
ENDIF
ALLOCATE(immspecname_dum(50))

IF (ALLOCATED(immspecname_lzs_dum)) THEN
  DEALLOCATE(immspecname_lzs_dum)
ENDIF
ALLOCATE(immspecname_lzs_dum(50))

IF (ALLOCATED(concfile_dum)) THEN
  DEALLOCATE(concfile_dum)
ENDIF
ALLOCATE(concfile_dum(50))

IF (ALLOCATED(lfile_conc_dum)) THEN
  DEALLOCATE(lfile_conc_dum)
ENDIF
ALLOCATE(lfile_conc_dum(50))

IF (ALLOCATED(FileFormatType_conc_dum)) THEN
  DEALLOCATE(FileFormatType_conc_dum)
ENDIF
ALLOCATE(FileFormatType_conc_dum(50))


REWIND nout

10  READ(nout,'(a)',END=1000) zone
tracer = tracer + 1
nlen1 = LEN(zone)
!!CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)

IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
  IF (ssch == 'read_immobile_species_file' .OR. ssch =='read_immobile_species_dataset') THEN
  
  readimmspec = .true.

!  Looking for a immspec name
    
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
    immspecname = ssch
    immspecname_lzs = lzs
    
    i = 0
    DO WHILE (i /= size(ulab))
    i = i+1
    IF (immspecname == ulab(i)) THEN
    !verify if primary species is part of immobile species
    IF (immobile_species(i) == 1) THEN
    CONTINUE
    ELSE
    WRITE(*,*)
    WRITE(*,*) ' Species ', immspecname(1: immspecname_lzs), ' provided in read_immobile_species_file is not part of the immobile species.'
    WRITE(*,*)
    STOP
    ENDIF
    dummyboolean = 1
    exit
    ENDIF
    END DO
    IF (dummyboolean == 1) THEN
    immspec_index = immspec_index + 1
    immspec_id_dum(immspec_index) = i
    immspecname_dum(immspec_index) = immspecname
    immspecname_lzs_dum(immspec_index) = lzs
    dummyboolean = 0
    WRITE(*,*)
    WRITE(*,*) ' Immobile species distribution of ', immspecname(1: immspecname_lzs), ' initialized from file.'
    WRITE(*,*)
    ELSE
    WRITE(*,*)
    WRITE(*,*) ' Could not find immobile species ', immspecname(1: immspecname_lzs), ' specified in read_immobile_species_file...'
    WRITE(*,*)
    STOP
    ENDIF
  
  !  Looking for immobile species file name
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL stringtype(ssch,lzs,res)
      concfile_dum(immspec_index) = ssch
      lfile_conc_dum(immspec_index) = ls
  
!!  *****************************************************************
!!    Now, check for a file format for immobile species
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      lenformat = ls
      CALL majuscules(ssch,lenformat)
      IF (ls /= 0) THEN
        IF (ssch == 'singlecolumn') THEN
          FileFormatType_conc_dum(immspec_index) = 'SingleColumn'
        ELSE IF (ssch == 'continuousread') THEN
          FileFormatType_conc_dum(immspec_index) = 'ContinuousRead'
        ELSE IF (ssch == 'unformatted') THEN
          FileFormatType_conc_dum(immspec_index) = 'Unformatted' 
        ELSE IF (ssch == 'distanceplusvariable' .OR. ssch == 'fullform' .OR. ssch == 'full') THEN
          FileFormatType_conc_dum(immspec_index) = 'FullForm'    
        ELSE IF (ssch == 'singlefile3d') THEN
          FileFormatType_conc_dum(immspec_index) = 'SingleFile3D'
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' File format not recognized: ', ssch(1:lenformat)
          WRITE(*,*)
          READ(*,*)
          STOP
        ENDIf
      ELSE    !! No file format provided, so assume default
        FileFormatType_conc_dum(immspec_index) = 'SingleColumn'
      ENDIF

  immspec_id = immspec_id_dum
  concfile = concfile_dum
  lfile_conc = lfile_conc_dum
  FileFormatType_conc = FileFormatType_conc_dum
  immspec_name = immspecname_dum
  immspec_name_length = immspecname_lzs_dum





    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No file name for conc for ', immspecname(1: immspecname_lzs),' in "read_immobile_species_file" in initial_condition section '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
  ELSE
  WRITE(*,*)
  WRITE(*,*) ' No information provided in "read_immobile_species_file" in initial_condition section '
  WRITE(*,*)
  READ(*,*)
  ENDIF  

  ELSE
    GO TO 10
  END IF
  GO TO 10

  
  
ELSE         ! No string found
  ! WRITE(*,*)
  ! WRITE(*,*) ' Information on VG alpha must be provided. '
  ! WRITE(*,*)
  ! READ(*,*)
  ! STOP
  GO TO 10
END IF



1000 RETURN
END SUBROUTINE read_immspecfile
