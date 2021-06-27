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


SUBROUTINE read_gaspump(nout,nx,ny,nz,nchem,ngaspump)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE transport
USE flow
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(IN)                                    :: nchem
INTEGER(I4B), INTENT(OUT)                                   :: ngaspump

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: nco
INTEGER(I4B)                                                :: intbnd_tmp
INTEGER(I4B)                                                :: jxxtemp
INTEGER(I4B)                                                :: jyytemp
INTEGER(I4B)                                                :: jzztemp

REAL(DP)                                                    :: qtemp

nxyz = nx*ny*nz

REWIND nout

ngaspump = 0
10 READ(nout,'(a)',END=500) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  IF (ssch == 'gaspump') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qtemp = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "gaspump"'
        WRITE(*,*) ' Looking for numerical value'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
!  Now, look for geochemical condition following pumping rate (only used if rate is positive)
      
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        
!  Check to see that heterogeneity label matches one of the labels
!  for geochemical conditions (condlabel)
        
        DO nco = 1,nchem
          IF (ssch == condlabel(nco)) THEN
            GO TO 50
          END IF
        END DO
        WRITE(*,*)
        WRITE(*,*) ' Geochemical condition for pumping well not found'
        WRITE(*,*) ' Label = ',ssch
        WRITE(*,*)
        READ(*,*)
        STOP
        50         ngaspump = ngaspump+ 1
        intbnd_tmp = nco
      ELSE         !  Blank string
        WRITE(*,*)
        WRITE(*,*) ' No geochemical condition for gas pumping well provided'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
! Now look for pumping well
      
      id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
      IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
          jxxtemp = JNUM(ssch)
        ELSE                !  An ascii string--so bag it.
          WRITE(*,*)
          WRITE(*,*) ' A grid location should follow gas pumping well specification'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      ELSE                  ! Zero length trailing string
        WRITE(*,*)
        WRITE(*,*) ' No grid location given for gas pumping zone'
        WRITE(*,*) ' Gas pumping zone ',ngaspump
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
      IF (jxxtemp > nx) THEN
        WRITE(*,*)
        WRITE(*,*) ' You have specified a gas pumping zone at JX > NX'
        WRITE(*,*) ' Gas pumping zone number ',ngaspump
        READ(*,*)
        STOP
      END IF
      IF (jxxtemp < 1) THEN
        WRITE(*,*)
        WRITE(*,*) ' You have specified a gas pumping zone at JX < 1'
        WRITE(*,*) ' Gas pumping zone number ',ngaspump
        READ(*,*)
        STOP
      END IF
      
      WRITE(*,*)
      WRITE(*,*) ' Gas pumping zone number ',ngaspump
      WRITE(*,*) ' Jxx location = ', jxxtemp
      WRITE(*,*)
      
!!      IF (ny > 1) THEN
        id = ids + ls
        CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF(ls /= 0) THEN
          lzs=ls
          CALL convan(ssch,lzs,res)
          IF (res == 'n') THEN
            jyytemp = JNUM(ssch)
          ELSE                !  An ascii string--so bag it.
            WRITE(*,*)
            WRITE(*,*) ' No Y location for gas pumping zone'
            WRITE(*,*) ' Gas pumping zone ',ngaspump
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
        ELSE                  ! Zero length trailing string
          WRITE(*,*)
          WRITE(*,*) ' No Y location for gas pumping zone'
          WRITE(*,*) ' Gas pumping zone ',ngaspump
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        
        IF (jyytemp > ny) THEN
          WRITE(*,*)
          WRITE(*,*) ' You have specified a gas pumping zone at JY > NY'
          WRITE(*,*) ' Gas pumping zone number ',ngaspump
          READ(*,*)
          STOP
        END IF
        IF (jyytemp < 1) THEN
          WRITE(*,*)
          WRITE(*,*) ' You have specified a gas pumping zone at JY < 1'
          WRITE(*,*) ' Gas pumping zone number ',ngaspump
          READ(*,*)
          STOP
        END IF
        
        WRITE(*,*)
        WRITE(*,*) ' Gas pumping zone number ',ngaspump
        WRITE(*,*) ' Jyy location = ', jyytemp
        WRITE(*,*)
        
!!      ELSE
!!        jyytemp = 1
!!      END IF


        id = ids + ls
        CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF(ls /= 0) THEN
          lzs=ls
          CALL convan(ssch,lzs,res)
          IF (res == 'n') THEN
            jzztemp = JNUM(ssch)
          ELSE                !  An ascii string--so bag it.
            WRITE(*,*)
            WRITE(*,*) ' No Z location for gas pumping zone'
            WRITE(*,*) ' Gas pumping zone ',ngaspump
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
        ELSE                  ! Zero length trailing string
          WRITE(*,*)
          WRITE(*,*) ' No Z location for gas pumping zone'
          WRITE(*,*) ' Gas pumping zone ',ngaspump
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        
        IF (jzztemp > nz) THEN
          WRITE(*,*)
          WRITE(*,*) ' You have specified a gas pumping zone at JZ > NZ'
          WRITE(*,*) ' Gas pumping zone number ',ngaspump
          READ(*,*)
          STOP
        END IF
        IF (jzztemp < 1) THEN
          WRITE(*,*)
          WRITE(*,*) ' You have specified a gas pumping zone at JZ < 1'
          WRITE(*,*) ' Gas pumping zone number ',ngaspump
          READ(*,*)
          STOP
        END IF
        
        WRITE(*,*)
        WRITE(*,*) ' Gas pumping zone number ',ngaspump
        WRITE(*,*) ' Jzz location = ', jzztemp
        WRITE(*,*)

      gaspump(1,jxxtemp,jyytemp,jzztemp) = qtemp
      intbndgas(1,jxxtemp,jyytemp,jzztemp) = intbnd_tmp
      
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No gas pumping rate given'
      WRITE(*,*) ' Gas pumping zone ignored'
      WRITE(*,*)
      GO TO 10
    END IF
  ELSE
    GO TO 10
  END IF
ELSE
  GO TO 10
END IF

GO TO 10

500 RETURN
END SUBROUTINE read_gaspump
