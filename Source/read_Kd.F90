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

    
SUBROUTINE read_kd(nout,ncomp,ncnt)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings
 
IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(OUT)                                   :: ncnt

!  Internal variables and arrays

LOGICAL(LGT)                                                :: speciesnotfound
CHARACTER (LEN=mls), DIMENSION(ncomp)                       :: nam_retard
CHARACTER (LEN=mls)                                         :: namtemp
CHARACTER (LEN=mls)                                          :: dumstring
REAL(DP), DIMENSION(ncomp)                                  :: tempkd

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lspecies
INTEGER(I4B)                                                :: kk
INTEGER(I4B)                                                :: lss

speciesnotfound = .true.
REWIND nout

ncnt = 0
100 READ(nout,'(a)',END=300) zone
dumstring = zone(1:13)
CALL majuscules(dumstring,13)
IF (dumstring == 'solid_density') GOTO 100
nlen1 = LEN(zone)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  ncnt = ncnt + 1
  IF (ncnt > ncomp) THEN
    WRITE(*,*)
    WRITE(*,*) ' Too many species for retardation specified'
    WRITE(*,*) ' Number of retarded species should not exceed # of components'
    WRITE(*,*) ' Number of components = ',ncomp
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  lzs=ls
  nam_retard(ncnt) = ssch
  namtemp = nam_retard(ncnt)
  CALL stringlen(namtemp,lspecies)
  
!  Now read the Kd (L/kg)
  
  id = ids + ls
  CALL sschaine(zone,id,iff,ssch,ids,ls)
  IF(ls /= 0) THEN
    lzs=ls
    CALL stringtype(ssch,lzs,res)
    IF (res == 'n') THEN
      tempkd(ncnt) = DNUM(ssch)
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Looking for a number for the Kd in RETARDATION section'
      WRITE(*,*) ' Following primary species ',namtemp(1:lspecies)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' No Kd following primary species name in section RETARDATION'
    WRITE(*,*) ' Primary species = ',namtemp(1:lspecies)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  
END IF

GO TO 100  ! Keep reading from file until the end of list is found

!  Go here when the list of retarded species is complete--check against primary species list

300 DO kk = 1,ncnt
  namtemp = nam_retard(kk)
  CALL stringlen(namtemp,lss)
  speciesnotfound = .true.
  DO i = 1,ncomp
    IF (ulab(i) == nam_retard(kk)) THEN
      distrib(i) = tempkd(kk)
      speciesnotfound = .false.
    END IF
  END DO
  IF (speciesnotfound) THEN
    WRITE(*,*)
    WRITE(*,*) ' ERROR in RETARDATION BLOCK'
    WRITE(*,*) ' Species associated with distribution coefficient not found in'
    WRITE(*,*) '   primary species list'
    WRITE(*,*) ' Looking for species: ',namtemp(1:lss)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END DO


END SUBROUTINE read_kd
