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
    
SUBROUTINE read_surface(nout,ncomp,nkin,nsurf)
USE crunchtype
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: nkin
INTEGER(I4B), INTENT(OUT)                                   :: nsurf

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: namtmp
CHARACTER (LEN=mls)                                         :: namsurfsave

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: lsurf
INTEGER(I4B)                                                :: lmin
INTEGER(I4B)                                                :: k

REWIND nout

nsurf = 0

100 READ(nout,'(a)',END=111) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  IF (ssch(1:1) == '>') THEN
    lsurf = ls
    nsurf = nsurf + 1
    IF (nsurf > msurf) THEN
      WRITE(*,*)
      WRITE(*,*) ' Number of surface sites dimensioned too small'
      WRITE(*,*) ' Msurf = ',msurf
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    namsurf(nsurf) = ssch
  ELSE
    GO TO 100
  END IF
  namsurfsave = namsurf(nsurf)

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF (ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (ssch /= 'on') THEN
      WRITE(*,*) ' Surface complex should be followed by "on" '
      WRITE(*,*) ' Surface complex: ',namsurfsave(1:lsurf)
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
      WRITE(*,*) ' In surface complexation block'
      WRITE(*,*) ' Surface complex: ',namsurfsave(1:lsurf)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    namtmp = ssch
    DO k = 1,nkin
      IF (namtmp == umin(k)) THEN
        ksurf(nsurf) = k
        GO TO 200
      END IF
    END DO
    WRITE(*,*)
    WRITE(*,*) ' Mineral substrate listed in surface complexation block'
    WRITE(*,*) '     not found in minerals list'
    WRITE(*,*) ' Surface complex: ',namsurfsave(1:lsurf)
    WRITE(*,*) ' Looking for mineral: ',namtmp(1:lmin)
    WRITE(*,*)
    READ(*,*)
    STOP
    
    200    CONTINUE
    
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Mineral host for surface hydroxyl site must be specified'
    WRITE(*,*) ' In surface complexation block'
    WRITE(*,*) ' Surface complex: ',namsurfsave(1:lsurf)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF

  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF (ls /= 0) THEN
    lzs=ls
    CALL convan(ssch,lzs,res)
    IF (ssch == '-no_edl' .OR. ssch == '--no_edl') THEN
      iedl(nsurf) = 1
    END IF
  ELSE
    GO TO 100
  END IF
  
END IF

GO TO 100

111  CONTINUE

RETURN

END SUBROUTINE read_surface
