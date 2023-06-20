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


SUBROUTINE read_het(nout,nchem,nhet,nx,ny,nz)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nchem
INTEGER(I4B), INTENT(OUT)                                   :: nhet
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: nco
INTEGER(I4B)                                                :: ls_a
INTEGER(I4B)                                                :: ls_b

nxyz = nx*ny*nz

ReadInitialConditions = .FALSE.

REWIND nout

nhet = 0

10  READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  
  lzs=ls
  CALL convan(ssch,lzs,res)
    
  IF (ssch == 'readinitialconditions' .or. ssch == 'readinitialcondition' .or. ssch=='ReadInitialConditions') THEN         
     !! Read initial condition distribution from file
    
    ReadInitialConditions = .TRUE.
    
!!  Looking for initial conditions file name 
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL stringtype(ssch,lzs,res)
      InitialConditionsFile = ssch
!!!      lfile = ls
    
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No file name following "ReadInitialConditions" in INITIAL_CONDITIONS section '
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    RETURN
  END IF

    
!     Check to see that heterogeneity label matches one of the labels
!       for geochemical conditions (condlabel)
  
  IF (ReadInitialConditions) THEN
      GO TO 50
  END IF
  
  DO nco = 1,nchem
    IF (ssch == condlabel(nco)) THEN
      GO TO 50
    END IF

  END DO
  IF (ssch(1:ls) == 'read_mineral_file') THEN
  GO TO 10
  ELSE
  WRITE(*,*)
  WRITE(*,*) ' Label for heterogeneity not found in list of condition labels'
  WRITE(*,*) ' Label = ',ssch(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
  ENDIF
  
50 nhet = nhet + 1
  
  ndist(nhet) = nco
  IF (nxyz == 1) THEN      ! If NXYZ = 1
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (ssch == 'fix') THEN
        jjfix(nhet) = 1
      ELSE
        jjfix(nhet) = 0
      END IF
    END IF
    jxxlo(nhet) = 1
    jxxhi(nhet) = 1
    jyylo(nhet) = 1
    jyyhi(nhet) = 1
    jzzlo(nhet) = 1
    jzzhi(nhet) = 1
    RETURN
  END IF
ELSE
  GO TO 10
END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF(ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jxxlo(nhet) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow heterogeneity label'
    WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
    WRITE(*,*) ' Heterogeneity number ',nhet
    WRITE(*,*) ' Dont know what to do with this string'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jxxhi(nhet) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow heterogeneity label'
      WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
      WRITE(*,*) ' Heterogeneity number ',nhet
      WRITE(*,*) ' Dont know what to do with this string'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    jxxhi(nhet) = jxxlo(nhet)   !  Assume jxxhi=jxxlo if no info provided
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' A grid location should follow heterogeneity label'
  WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
  WRITE(*,*) ' Heterogeneity number ',nhet
  WRITE(*,*) ' Zero length string following label'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

!!!IF (ny == 1 .AND. nz == 1) THEN
!!!  id = ids + ls
!!!  CALL sschaine(zone,id,iff,ssch,ids,ls)
!!!  IF (ls /= 0) THEN
!!!    lzs=ls
!!!    CALL convan(ssch,lzs,res)
!!!    IF (ssch == 'fix') THEN
!!!      jjfix(nhet) = 1
!!!    ELSE
!!!      jjfix(nhet) = 0
!!!    END IF
!!!  END IF
!!!  jyylo(nhet) = 1
!!!  jyyhi(nhet) = 1
!!!  jzzlo(nhet) = 1
!!!  jzzhi(nhet) = 1
!!!  GO TO 10
!!!END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF(ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jyylo(nhet) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow heterogeneity label'
    WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
    WRITE(*,*) ' Heterogeneity number ',nhet
    WRITE(*,*) ' Dont know what to do with this string'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jyyhi(nhet) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow heterogeneity label'
      WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
      WRITE(*,*) ' Heterogeneity number ',nhet
      WRITE(*,*) ' Dont know what to do with this string'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    jyyhi(nhet) = jyylo(nhet)   !  Assume jyyhi=jyylo if no info provided
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' A grid location should follow heterogeneity label'
  WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
  WRITE(*,*) ' Heterogeneity number ',nhet
  WRITE(*,*) ' Zero length string following label'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

!!!IF (nz == 1) THEN
!!!  id = ids + ls
!!!  CALL sschaine(zone,id,iff,ssch,ids,ls)
!!!  IF (ls /= 0) THEN
!!!    lzs=ls
!!!    CALL convan(ssch,lzs,res)
!!!    IF (ssch == 'fix') THEN
!!!      jjfix(nhet) = 1
!!!    ELSE
!!!      jjfix(nhet) = 0
!!!    END IF
!!!  END IF
!!!  jzzlo(nhet) = 1
!!!  jzzhi(nhet) = 1
!!!  GO TO 10
!!!END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF(ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jzzlo(nhet) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow heterogeneity label'
    WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
    WRITE(*,*) ' Heterogeneity number ',nhet
    WRITE(*,*) ' Dont know what to do with this string'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jzzhi(nhet) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow heterogeneity label'
      WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
      WRITE(*,*) ' Heterogeneity number ',nhet
      WRITE(*,*) ' Dont know what to do with this string'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    jzzhi(nhet) = jzzlo(nhet)   !  Assume jzzhi=jzzlo if no info provided
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' A grid location should follow heterogeneity label'
  WRITE(*,*) ' For heterogeneity labelled ',condlabel(nco)
  WRITE(*,*) ' Heterogeneity number ',nhet
  WRITE(*,*) ' Zero length string following label'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF (ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'fix') THEN
    jjfix(nhet) = 1
  ELSE
    jjfix(nhet) = 0
  END IF
END IF


GO TO 10

500  RETURN
END SUBROUTINE read_het
