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

    
SUBROUTINE read_toperatio(nout,ncomp,ndecay,idecayall,idpoint,kd,nchem,topefound)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE mineral
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: ncomp
INTEGER(I4B), INTENT(IN)                                    :: ndecay
INTEGER(I4B), INTENT(IN)                                    :: idecayall
INTEGER(I4B), INTENT(IN)                                    :: idpoint
INTEGER(I4B), INTENT(IN OUT)                                :: kd
INTEGER(I4B), INTENT(IN)                                    :: nchem
LOGICAL(LGT), INTENT(IN OUT)                                :: topefound

!  Internal variables and arrays

CHARACTER (LEN=mls)                                         :: dumstring
CHARACTER (LEN=mls)                                         :: toperatio

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: ipoint
INTEGER(I4B)                                                :: nisot
INTEGER(I4B)                                                :: lcond
INTEGER(I4B)                                                :: ltope
INTEGER(I4B)                                                :: isotope

REAL(DP)                                                    :: sum
REAL(DP)                                                    :: eps


eps = 1.e-10

ipoint = idecay(idpoint)
nisot = nisotope(idpoint)

dumstring = condlabel(nchem)
CALL stringlen(dumstring,lcond)
toperatio = ulab(ipoint)
CALL stringlen(toperatio,ltope)

allocate(namtope(nisot))
allocate(temp_ratio(nisot))

REWIND nout

100 READ(nout,'(a)',END=300) zone
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'isotope_ratio') THEN    !  "isotope_ratio" label found
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF (ls /= 0) THEN
      lzs = ls
      IF (ssch == toperatio) THEN   !  Decay element found
        id = ids + ls
        CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF(ls /= 0) THEN
          lzs=ls
          IF (idecayall == 1) THEN
            IF (ssch == 'all' .OR. ssch == 'All' .OR. ssch == 'ALL') THEN
              OPEN(UNIT=nscratch,STATUS='scratch')
              REWIND nscratch
              WRITE(nscratch,*) zone(ids+ls+1:mls)
              REWIND nscratch
              READ(nscratch,*,ERR=500) (namtope(isotope),temp_ratio(isotope),  &
                  isotope=1,nisot)
              CLOSE(UNIT=nscratch)
              
!  Check first to see that ratios add up to 1
              
              sum = 0.0
              DO isotope = 1,nisot
                sum = sum + temp_ratio(isotope)
              END DO
              IF (sum < (1.0-eps) .OR. sum > (1.0+eps)) THEN
                WRITE(*,*)
                WRITE(*,*) ' Isotope ratios dont sum to 1'
                WRITE(*,*) ' Isotope labels dont match for: ',toperatio(1:ltope)
                WRITE(*,*) ' In condition: ',dumstring(1:lcond)
                WRITE(*,*)
                STOP
              END IF
              
!  Check to see that names of isotopes match what was given in the "DECAY" block
!  Assumes (for now) that the order is the same, otherwise...
              
              DO isotope=1,nisot
                IF (namtope(isotope) /= decay_label(isotope,idpoint)) THEN
                  WRITE(*,*)
                  WRITE(*,*) ' Isotope labels dont match for: ',toperatio(1:ltope)
                  WRITE(*,*) ' In condition: ',dumstring(1:lcond)
                  WRITE(*,*)
                  STOP
                END IF
              END DO
              topefound = .true.
              DO kd = 1,nmindecay(idpoint)
                DO isotope=1,nisot
                  ratio_isotope_init(isotope,kd,idpoint,nchem) = temp_ratio(isotope)
                END DO
              END DO
              
            ELSE     !  Label is not "all"
              GO TO 100
            END IF
            
          ELSE    !  "idecayall" equals 0, so don't search for "all" string
            
!  Check that we have the right mineral name
            
            IF (ssch == umin(kdecay(kd,idpoint))) THEN
              
              OPEN(UNIT=nscratch,STATUS='scratch')
              REWIND nscratch
              WRITE(nscratch,*) zone(ids+ls+1:mls)
              REWIND nscratch
              READ(nscratch,*,ERR=500) (namtope(isotope),temp_ratio(isotope),  &
                  isotope=1,nisot)
              CLOSE(UNIT=nscratch)
              
              sum = 0.0
              DO isotope = 1,nisot
                sum = sum + temp_ratio(isotope)
              END DO
              IF (sum < (1.0-eps) .OR. sum > (1.0+eps)) THEN
                WRITE(*,*)
                WRITE(*,*) ' Isotope ratios dont sum to 1'
                WRITE(*,*) ' Isotope labels dont match for: ',toperatio(1:ltope)
                WRITE(*,*) ' In condition: ',dumstring(1:lcond)
                WRITE(*,*)
                STOP
              END IF
              
!  Check to see that names of isotopes match what was given in the "DECAY" block
!  Assumes (for now) that the order is the same, otherwise...
              
              DO isotope=1,nisot
                IF (namtope(isotope) /= decay_label(isotope,idpoint)) THEN
                  WRITE(*,*)
                  WRITE(*,*) ' Isotope labels dont match for: ',toperatio(1:ltope)
                  WRITE(*,*) ' In condition: ',dumstring(1:lcond)
                  WRITE(*,*)
                  STOP
                END IF
              END DO
              topefound = .true.
              DO isotope=1,nisot
                ratio_isotope_init(isotope,kd,idpoint,nchem) = temp_ratio(isotope)
              END DO
            ELSE
              GO TO 100
            END IF
            
          END IF
          
          deallocate(namtope)
          deallocate(temp_ratio)
          RETURN    !  Information found, so exit subroutine
          
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Cannot find isotopic ratio for: ',toperatio(1:ltope)
          WRITE(*,*) ' In condition: ',dumstring(1:lcond)
          WRITE(*,*)
          STOP
        END IF
      END IF
      GO TO 100     !  Continue looking
      
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Label "isotope_ratio" should be followed by radioactive element name'
      WRITE(*,*) ' In condition: ',dumstring(1:lcond)
      WRITE(*,*)
      STOP
    END IF
    
    
  ELSE
    GO TO 100
  END IF
ELSE
  GO TO 100
END IF


300 deallocate(namtope)
deallocate(temp_ratio)
RETURN

500 WRITE(*,*)
WRITE(*,*) ' Error reading isotopic ratio for: ',toperatio(1:ltope)
WRITE(*,*) ' In condition: ',dumstring(1:lcond)
WRITE(*,*)
STOP

END SUBROUTINE read_toperatio
