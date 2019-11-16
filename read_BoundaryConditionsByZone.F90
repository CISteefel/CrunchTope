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

    
SUBROUTINE read_BoundaryConditionByZone(nout,nx,ny,nz,nBoundaryConditionZone)
USE crunchtype
USE CrunchFunctions
USE runtime
USE params
USE transport
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(OUT)                                   :: nBoundaryConditionZone

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: ls_a
INTEGER(I4B)                                                :: ls_b
INTEGER(I4B)                                                :: l
INTEGER(I4B)                                                :: ncond


CHARACTER (LEN=mls)                                         :: GeochemicalConditionName

INTEGER(I4B), PARAMETER                                     :: mBoundaryConditionZone=500

IF (ALLOCATED(BoundaryConditionName)) THEN
  DEALLOCATE(BoundaryConditionName)
END IF
ALLOCATE( BoundaryConditionName(0:mBoundaryConditionZone) )

IF (ALLOCATED(jxxBC_lo)) THEN
  DEALLOCATE(jxxBC_lo)
END IF
ALLOCATE(jxxBC_lo(mBoundaryConditionZone))
IF (ALLOCATED(jyyBC_lo)) THEN
  DEALLOCATE(jyyBC_lo)
END IF
ALLOCATE(jyyBC_lo(mBoundaryConditionZone))
IF (ALLOCATED(jzzBC_lo)) THEN
  DEALLOCATE(jzzBC_lo)
END IF
ALLOCATE(jzzBC_lo(mBoundaryConditionZone))
IF (ALLOCATED(jxxBC_hi)) THEN
  DEALLOCATE(jxxBC_hi)
END IF
ALLOCATE(jxxBC_hi(mBoundaryConditionZone))
IF (ALLOCATED(jyyBC_hi)) THEN
  DEALLOCATE(jyyBC_hi)
END IF
ALLOCATE(jyyBC_hi(mBoundaryConditionZone))
IF (ALLOCATED(jzzBC_hi)) THEN
  DEALLOCATE(jzzBC_hi)
END IF
ALLOCATE(jzzBC_hi(mBoundaryConditionZone))

nxyz = nx*ny*nz
BoundaryConditionName = ' '

REWIND nout

nBoundaryConditionZone = 0

DO nCond = 1,mBoundaryConditionZone

  READ(nout,'(a)',END=500) zone

  nlen1 = LEN(zone)
  CALL majuscules(zone,nlen1)
  id = 1
  iff = mls
  CALL sschaine(zone,id,iff,ssch,ids,ls)

  lzs=ls
  CALL convan(ssch,lzs,res)

  IF (ssch == 'boundarycondition' .OR. ssch == 'BoundaryCondition' .OR. ssch == 'BOUNDARYCONDITION') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)

      GeochemicalConditionName = ssch
      
!     Now look for ASCII string indicating location of permeability
      
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)

      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'a') THEN
          IF (ssch == 'zone') THEN
            
!  "Zone" specified, so look for locations
            
            nBoundaryConditionZone = nBoundaryConditionZone + 1

            IF (nBoundaryConditionZone > mBoundaryConditionZone) THEN

              WRITE(*,*)
              WRITE(*,*)  ' Number of BoundaryCondition zones dimensioned too small'
              WRITE(*,*)  ' Number of BoundaryCondition zones = ', nBoundaryConditionZone
              WRITE(*,*)  ' Dimension of BoundaryCondition zones = ', mBoundaryConditionZone
              WRITE(*,*)  ' Contact C.I. Steefel at Berkeley Lab: "CISteefel@lbl.gov"'
              WRITE(*,*)
              STOP

            END IF
            
            BoundaryConditionName(nBoundaryConditionZone) = GeochemicalConditionName
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jxxBC_lo(nBoundaryConditionZone) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' A grid location should follow zone specification'
                WRITE(*,*) ' Dont know what to do with this string'
                WRITE(*,*)
                STOP
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jxxBC_hi(nBoundaryConditionZone) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "boundarycondition"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jxxBC_hi(nBoundaryConditionZone) = jxxBC_lo(nBoundaryConditionZone) 
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No X or Y grid location given for BoundaryCondition'
              WRITE(*,*) ' BoundaryCondition zone ',nBoundaryConditionZone
              WRITE(*,*)
              STOP
            END IF
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jyyBC_lo(nBoundaryConditionZone) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' No Y location for BoundaryCondition '
                WRITE(*,*)
                STOP
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jyyBC_hi(nBoundaryConditionZone) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "BoundaryCondition"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jyyBC_hi(nBoundaryConditionZone) = jyyBC_lo(nBoundaryConditionZone)  
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No Y location given for BoundaryCondition zone'
              WRITE(*,*) ' BoundaryCondition zone number ',nBoundaryConditionZone
              WRITE(*,*)
              STOP
            END IF    
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jzzBC_lo(nBoundaryConditionZone) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' No Z location for BoundaryCondition '
                WRITE(*,*)
                STOP
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jzzBC_hi(nBoundaryConditionZone) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "BoundaryCondition"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jzzBC_hi(nBoundaryConditionZone) = jzzBC_lo(nBoundaryConditionZone)   !  Assume jxxpermx_hi=jxxpermx_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No Z location given for BoundaryCondition zone'
              WRITE(*,*) ' BoundaryCondition zone number ',nBoundaryConditionZone
              WRITE(*,*)
              STOP
            END IF    

          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Dont understand string following BoundaryCondition specification'
            WRITE(*,*)   ssch(1:ls)
            WRITE(*,*)
            STOP
          END IF
          
        ELSE                !  A number--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' Cant interpret string following BoundaryCondition value'
          WRITE(*,*) ' Looking for an ASCII string'
          WRITE(*,*)
          STOP
        END IF

      ELSE   ! If no BoundaryCondition input, then assume X_begin and X-end format

        BoundaryConditionName(0) = 'none'

      END IF

    END IF

  END IF
  

END DO

500 DO l = 1,nBoundaryConditionZone
  IF (jxxBC_hi(l) > nx+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JX > NX+1'
    WRITE(*,*)
    STOP
  END IF
  IF (jyyBC_hi(l) > ny+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JY > NY'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzBC_hi(l) > nz+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JZ > NZ'
    WRITE(*,*)
    STOP
  END IF
  IF (jxxBC_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JX < 0'
    WRITE(*,*)
    STOP
  END IF
  IF (jyyBC_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JY < 1'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzBC_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a BoundaryCondition at JZ < 1'
    STOP
  END IF
END DO

RETURN
END SUBROUTINE read_BoundaryConditionByZone
