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

    
SUBROUTINE read_permx(nout,nx,ny,nz,npermx)
USE crunchtype
USE CrunchFunctions
USE params
USE flow
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(OUT)                                   :: npermx

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

REAL(DP)                                                    :: permx_tmp

nxyz = nx*ny*nz

permzonex(0) = 0.0
REWIND nout

npermx = 0
10 READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'permeability_x') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        permx_tmp = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "permeability"'
        WRITE(*,*) ' Looking for numerical value'
        WRITE(*,*)
        STOP
      END IF
      
! Now look for ASCII string indicating location of permeability
      
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'a') THEN
          IF (ssch == 'default' .OR. ssch == 'all') THEN
            permzonex(0) = permx_tmp
          ELSE IF (ssch == 'zone') THEN
            
!  "Zone" specified, so look for locations
            
            npermx = npermx + 1
            IF (npermx > mperm) THEN
              WRITE(*,*)
              WRITE(*,*)  ' Number of permeability zones dimensioned too small'
              WRITE(*,*)  ' Number of permeability zones = ',npermx
              WRITE(*,*)  ' Dimension of permeability zones = ',mperm
              WRITE(*,*)  ' Contact the code developer at CISteefel@lbl.gov '
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
            
            permzonex(npermx) = permx_tmp
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jxxpermx_lo(npermx) = JNUM(ssch_a)
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
                  jxxpermx_hi(npermx) = JNUM(ssch_b)

                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "permeability"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jxxpermx_hi(npermx) = jxxpermx_lo(npermx)   !  Assume jxxpermx_hi=jxxpermx_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No X or Y grid location given for permeability'
              WRITE(*,*) ' Permeability zone ',npermx
              WRITE(*,*)
              STOP
            END IF
            
            WRITE(*,*)
            WRITE(*,*) ' X permeability zone number ',npermx
            WRITE(*,*) ' Jxxpermx_lo = ', jxxpermx_lo(npermx)
            WRITE(*,*) ' Jxxpermx_hi = ',jxxpermx_hi(npermx)
            WRITE(*,*)
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jyypermx_lo(npermx) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' No Y location for permeability '
                WRITE(*,*)
                STOP
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jyypermx_hi(npermx) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "permeability"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jyypermx_hi(npermx) = jyypermx_lo(npermx)   !  Assume jxxpermx_hi=jxxpermx_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No Y location given for permeability zone'
              WRITE(*,*) ' Permeability zone number ',npermx
              WRITE(*,*)
              STOP
            END IF
              
            WRITE(*,*)
            WRITE(*,*) ' Jyypermx_lo = ',jyypermx_lo(npermx)
            WRITE(*,*) ' Jyypermx_hi = ',jyypermx_hi(npermx)
            WRITE(*,*)    
            
            id = ids + ls
            CALL sschaine_hyph(zone,id,iff,ssch_a,ssch_b,ids,ls_a,ls_b,ls)
            IF(ls /= 0) THEN
              lzs=ls_a
              CALL convan(ssch_a,lzs,res)
              IF (res == 'n') THEN
                jzzpermx_lo(npermx) = JNUM(ssch_a)
              ELSE                !  An ascii string--so bag it.
                WRITE(*,*)
                WRITE(*,*) ' No Z location for permeability '
                WRITE(*,*)
                STOP
              END IF
              IF (ls_b /= 0) THEN
                lzs=ls_b
                CALL convan(ssch_b,lzs,res)
                IF (res == 'n') THEN
                  jzzpermx_hi(npermx) = JNUM(ssch_b)
                ELSE                !  An ascii string--so bag it.
                  WRITE(*,*)
                  WRITE(*,*) ' A grid location should follow zone specification'
                  WRITE(*,*) ' Dont know what to do with this string after "permeability"'
                  WRITE(*,*)
                  STOP
                END IF
              ELSE
                jzzpermx_hi(npermx) = jzzpermx_lo(npermx)   !  Assume jxxpermx_hi=jxxpermx_lo
              END IF
            ELSE                  ! Zero length trailing string
              WRITE(*,*)
              WRITE(*,*) ' No Z location given for permeability zone'
              WRITE(*,*) ' Permeability zone number ',npermx
              WRITE(*,*)
              STOP
            END IF
              
            WRITE(*,*)
            WRITE(*,*) ' Jzzpermx_lo = ',jzzpermx_lo(npermx)
            WRITE(*,*) ' Jzzpermx_hi = ',jzzpermx_hi(npermx)
            WRITE(*,*)    

          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Dont understand string following permeability value'
            WRITE(*,*) ssch(1:ls)
            WRITE(*,*)
            STOP
          END IF
          
        ELSE                !  A number--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' Cant interpret string following permeability value'
          WRITE(*,*) ' Looking for ASCII string'
          WRITE(*,*)
          STOP
        END IF
      ELSE   ! Assume this is default if nothing else given
        permzonex(0) = permx_tmp
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value given for permeability'
      WRITE(*,*) ' Permeability specification ignored'
      WRITE(*,*)
    END IF
  ELSE
    GO TO 10
  END IF
  
END IF

GO TO 10

500 DO l = 1,npermx
  IF (jxxpermx_hi(l) > nx+1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an X permeability at JX > NX+1'
    WRITE(*,*)
    STOP
  END IF
  IF (jyypermx_hi(l) > ny) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an X permeability at JY > NY'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzpermx_hi(l) > nz) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an X permeability at JZ > NZ'
    WRITE(*,*)
    STOP
  END IF
  IF (jxxpermx_lo(l) < 0) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an X permeability at JX < 0'
    WRITE(*,*)
    STOP
  END IF
  IF (jyypermx_lo(l) < 1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an X permeability at JY < 1'
    WRITE(*,*)
    STOP
  END IF
  IF (jzzpermx_lo(l) < 1) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified an X permeability at JZ < 1'
    STOP
  END IF
END DO

RETURN
END SUBROUTINE read_permx
