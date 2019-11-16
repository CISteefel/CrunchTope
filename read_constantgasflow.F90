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

SUBROUTINE read_constantgasflow(nout,nx,ny,nz,constant_gasflow,  &
    qxgasinit,qygasinit,qzgasinit)
USE crunchtype
USE CrunchFunctions
USE params
USE concentration
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
LOGICAL(LGT), INTENT(OUT)                                   :: constant_gasflow
REAL(DP), INTENT(OUT)                                       :: qxgasinit
REAL(DP), INTENT(OUT)                                       :: qygasinit
REAL(DP), INTENT(OUT)                                       :: qzgasinit

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1

nxyz = nx*ny*nz
constant_gasflow = .false.
qxgasinit = 0.0
qygasinit = 0.0
qzgasinit = 0.0

REWIND nout

10  READ(nout,'(a)',END=1000) zone
nlen1 = LEN(zone)
CALL majuscules(zone,nlen1)
id = 1
iff = mls
CALL sschaine(zone,id,iff,ssch,ids,ls)
IF(ls /= 0) THEN
  lzs=ls
  CALL convan(ssch,lzs,res)
  
!  Check to see if initial substring is "constant_gasflow"
  
  IF (ssch == 'constant_gasflow' .OR. ssch == 'constant_gas_flow') THEN
    IF (nxyz == 1) THEN    ! No flow in case of NXYZ = 1
      qxgasinit = 0.0
      qygasinit = 0.0
      qzgasinit = 0.0
      constant_gasflow = .true.
      RETURN
    END IF
    
!  First, looking for X velocity
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qxgasinit = DNUM(ssch)
        constant_gasflow = .true.
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following "constant_gasflow" '
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value following "constant_gasflow" '
      WRITE(*,*) ' Value(s) must be given if constant_gasflow specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
!  Now, looking for Y velocity if NY not equal to 1
    
    IF (ny == 1) THEN
      qygasinit = 0.0
      qzgasinit = 0.0
      GO TO 100
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qygasinit = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following X gas velocity'
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value for Y gas velocity given'
      WRITE(*,*) ' Y velocity must be given if constant_gasflow specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
! Now, looking for Z velocity if NZ not equal to 1
    
    IF (nz == 1) THEN
      qzgasinit = 0.0
      GO TO 100
    END IF
    
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
      lzs=ls
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
        qzgasinit = DNUM(ssch)
      ELSE                !  An ascii string--so bag it.
        WRITE(*,*)
        WRITE(*,*) ' Cant interpret string following Y gas velocity'
        WRITE(*,*) ' A numerical value should follow'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' No value for Z gas velocity given'
      WRITE(*,*) ' Z velocity must be given if constant_gasflow specified'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    
    100     CONTINUE
    
  ELSE
    GO TO 10
  END IF
  GO TO 10
  
ELSE         ! No string found
  GO TO 10
END IF

1000 RETURN
END SUBROUTINE read_constantgasflow
