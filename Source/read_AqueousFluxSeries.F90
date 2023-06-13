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

    
SUBROUTINE read_AqueousFluxSeries(nout,ncomp,nx,ny,nz)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings
USE runtime

IMPLICIT NONE

INTERFACE
  SUBROUTINE read_multstring(nout,lchar,parchar,parfind, stringarray,lenarray,section)
  USE crunchtype
  USE params
  USE strings
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                                    :: nout
  INTEGER(I4B), INTENT(OUT)                                   :: lchar
  CHARACTER (LEN=mls), INTENT(IN)                             :: parchar
  CHARACTER (LEN=mls), INTENT(IN OUT)                         :: parfind
  CHARACTER (LEN=mls), DIMENSION(:), INTENT(IN OUT)           :: stringarray
  INTEGER(I4B), INTENT(OUT)                                   :: lenarray
  CHARACTER (LEN=mls), INTENT(IN)                             :: section
  END SUBROUTINE read_multstring
END INTERFACE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
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
INTEGER(I4B)                                                :: i
INTEGER(I4B)                                                :: ik
INTEGER(I4B)                                                :: lchar
INTEGER(I4B)                                                :: lenarray
INTEGER(I4B)                                                :: ll
INTEGER(I4B)                                                :: ls_a
INTEGER(I4B)                                                :: ls_b
INTEGER(I4B)                                                :: nlength

CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: workchar1

!!INTEGER(I4B)                                               :: lchar
CHARACTER (LEN=mls)                                        :: parchar
CHARACTER (LEN=mls)                                        :: parfind
!!CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE             :: stringarray
!!INTEGER(I4B)                                               :: lenarray
CHARACTER (LEN=mls)                                        :: section
CHARACTER (LEN=mls)                                        :: dumstring

nAqueousFluxSeriesFile = 0

IF (ALLOCATED(AqueousFluxSeriesFile)) THEN
  DEALLOCATE(AqueousFluxSeriesFile)
END IF
ALLOCATE(AqueousFluxSeriesFile(1))

IF (.NOT. ALLOCATED(stringarray)) THEN
  ALLOCATE(stringarray(100))
END IF

IF (ALLOCATED(jxAqueousFluxSeries_lo)) THEN
  DEALLOCATE(jxAqueousFluxSeries_lo)
END IF
ALLOCATE(jxAqueousFluxSeries_lo(100))
IF (ALLOCATED(jxAqueousFluxSeries_hi)) THEN
  DEALLOCATE(jxAqueousFluxSeries_hi)
END IF
ALLOCATE(jxAqueousFluxSeries_hi(100))
IF (ALLOCATED(jyAqueousFluxSeries_lo)) THEN
  DEALLOCATE(jyAqueousFluxSeries_lo)
END IF
ALLOCATE(jyAqueousFluxSeries_lo(100))
IF (ALLOCATED(jyAqueousFluxSeries_hi)) THEN
  DEALLOCATE(jyAqueousFluxSeries_hi)
END IF
ALLOCATE(jyAqueousFluxSeries_hi(100))
IF (ALLOCATED(jzAqueousFluxSeries_lo)) THEN
  DEALLOCATE(jzAqueousFluxSeries_lo)
END IF
ALLOCATE(jzAqueousFluxSeries_lo(100))
IF (ALLOCATED(jzAqueousFluxSeries_hi)) THEN
  DEALLOCATE(jzAqueousFluxSeries_hi)
END IF
ALLOCATE(jzAqueousFluxSeries_hi(100))


nxyz = nx*ny*nz

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
  
!  Check to see if initial substring is "WriteAqueousFlux" (could be more than 1)
  
  IF (ssch == 'writeaqueousflux') THEN
    nAqueousFluxSeriesFile = nAqueousFluxSeriesFile + 1
  ELSE
    GO TO 10
  END IF
ELSE         ! No string found
  GO TO 10
END IF

id = ids + ls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  i = size(AqueousFluxSeriesFile,1)
  ALLOCATE(workchar1(i))
  workchar1 = AqueousFluxSeriesFile
  DEALLOCATE(AqueousFluxSeriesFile)
  ALLOCATE(AqueousFluxSeriesFile(nAqueousFluxSeriesFile))
  IF(nAqueousFluxSeriesFile /= 0) AqueousFluxSeriesFile(1:nAqueousFluxSeriesFile-1) = workchar1(1:nAqueousFluxSeriesFile-1)
  DEALLOCATE(workchar1)
  AqueousFluxSeriesFile(nAqueousFluxSeriesFile) = ssch
ELSE
  IF (nAqueousFluxSeriesFile == 1) THEN
    AqueousFluxSeriesFile(1) = 'AqueousFlux.out'
  ELSE
    WRITE(*,*) 
    WRITE(*,*) ' File name for AqueousFluxSeriesFile required when more than one is used'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

IF (nxyz == 1) THEN
  jxAqueousFluxSeries_lo(nAqueousFluxSeriesFile) = 1
  jyAqueousFluxSeries_lo(nAqueousFluxSeriesFile) = 1
  jzAqueousFluxSeries_lo(nAqueousFluxSeriesFile) = 1
  jxAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = 1
  jyAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = 1
  jzAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF (ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jxAqueousFluxSeries_lo(nAqueousFluxSeriesFile) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow specification of AqueousFluxSeriesFile'
    WRITE(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jxAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow specification of AqueousFluxSeriesFile'
      WRITE(*,*)
      STOP
    END IF
  ELSE
    jxAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = jxAqueousFluxSeries_lo(nAqueousFluxSeriesFile)   !  Assume jxAqueousFluxSeries_hi=jxAqueousFluxSeries_lo
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' No X or Y grid location given for nAqueousFluxSeriesFile'
  WRITE(*,*) ' AqueousFluxSeries zone ',nAqueousFluxSeriesFile
  WRITE(*,*)
  STOP
END IF

!!!  Got to here

IF (ny == 1) THEN
  jyAqueousFluxSeries_lo(nAqueousFluxSeriesFile) = 1
  jzAqueousFluxSeries_lo(nAqueousFluxSeriesFile) = 1
  jyAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = 1
  jzAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF (ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jyAqueousFluxSeries_lo(nAqueousFluxSeriesFile) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow specification of AqueousFluxSeriesFile'
    WRITE(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jyAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow specification of AqueousFluxSeriesFile'
      WRITE(*,*)
      STOP
    END IF
  ELSE
    jyAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = jyAqueousFluxSeries_lo(nAqueousFluxSeriesFile)   !  Assume jxAqueousFluxSeries_hi=jxAqueousFluxSeries_lo
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' No X or Y grid location given for nAqueousFluxSeriesFile'
  WRITE(*,*) ' AqueousFluxSeries zone ',nAqueousFluxSeriesFile
  WRITE(*,*)
  STOP
END IF

IF (nz == 1) THEN
  jzAqueousFluxSeries_lo(nAqueousFluxSeriesFile) = 1
  jzAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = 1
  GO TO 500
END IF

id = ids + ls
CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
IF (ls /= 0) THEN
  lzs=ls_a
  CALL convan(ssch_a,lzs,res)
  IF (res == 'n') THEN
    jzAqueousFluxSeries_lo(nAqueousFluxSeriesFile) = JNUM(ssch_a)
  ELSE                !  An ascii string--so bag it.
    WRITE(*,*)
    WRITE(*,*) ' A grid location should follow specification of AqueousFluxSeriesFile'
    WRITE(*,*)
    STOP
  END IF
  IF (ls_b /= 0) THEN
    lzs=ls_b
    CALL convan(ssch_b,lzs,res)
    IF (res == 'n') THEN
      jzAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = JNUM(ssch_b)
    ELSE                !  An ascii string--so bag it.
      WRITE(*,*)
      WRITE(*,*) ' A grid location should follow specification of AqueousFluxSeriesFile'
      WRITE(*,*)
      STOP
    END IF
  ELSE
    jzAqueousFluxSeries_hi(nAqueousFluxSeriesFile) = jzAqueousFluxSeries_lo(nAqueousFluxSeriesFile)   !  Assume jxAqueousFluxSeries_hi=jxAqueousFluxSeries_lo
  END IF
ELSE                  ! Zero length trailing string
  WRITE(*,*)
  WRITE(*,*) ' No X or Y grid location given for nAqueousFluxSeriesFile'
  WRITE(*,*) ' AqueousFluxSeries zone ',nAqueousFluxSeriesFile
  WRITE(*,*)
  STOP
END IF

!! Go to next line to see if there is another AqueousFluxSeriesFile file
500 GO TO 10

1000  CONTINUE

REWIND(NOUT)

!! Now read the species names to be written to the AqueousFluxSeries file

parchar = 'FluxSpecies'
parfind = ' '
lenarray = 0
nplotAqueousFlux = 0
section = 'OUTPUT'
  
CALL read_multstring(nout,lchar,parchar,parfind,stringarray,lenarray,section)
IF (parfind == 'FluxSpecies' .OR. parfind == 'fluxspecies') THEN   

!!  Check to see that strings match species names
      
  nplotAqueousFlux = lenarray

  ALLOCATE(iplotAqueousflux(nplotAqueousFlux))

  DO ll = 1,nplotAqueousFlux
    iplotAqueousFlux(ll) = 0
    IF (stringarray(ll) == 'pH' .OR. stringarray(ll) == 'ph') THEN
      WRITE(*,*) '       ERROR '
      WRITE(*,*) ' Fluxes should be written in terms of total concentrations (not pH)'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    DO ik = 1,ncomp
      IF (ulab(ik) == stringarray(ll)) THEN
        iplotAqueousFlux(ll) = ik
      END IF
    END DO
    IF (iplotAqueousFlux(ll) == 0) THEN
      dumstring = stringarray(ll)
      CALL stringlen(dumstring,nlength)
      WRITE(*,*)
      WRITE(*,*) ' Aqueous flux species not found in list: ',dumstring(1:nlength)
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  END DO

ELSE

  nplotAqueousFlux = 0

END IF
  
IF (ALLOCATED(AqueousFluxSeriesSpecies)) THEN
  DEALLOCATE(AqueousFluxSeriesSpecies)
END IF
ALLOCATE(AqueousFluxSeriesSpecies(nplotAqueousFlux))
AqueousFluxSeriesSpecies = stringarray(1:nplotAqueousFlux)

!!DEALLOCATE(stringarray)

RETURN
END SUBROUTINE read_AqueousFluxSeries
