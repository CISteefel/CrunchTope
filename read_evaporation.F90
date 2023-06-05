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


SUBROUTINE read_evaporation(nout,nx,ny,nz,evapofile,rate,lfile,boolfix,booltimeseries,tslength)
USE crunchtype
USE CrunchFunctions
USE params
USE flow
USE strings
USE medium

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(OUT)                                   :: lfile
INTEGER(I4B), INTENT(OUT)                                   :: tslength
CHARACTER (LEN=mls), INTENT(OUT)                            :: evapofile
REAL(DP), INTENT(OUT)                                       :: rate
LOGICAL(LGT), INTENT(IN OUT)                                :: boolfix
LOGICAL(LGT), INTENT(IN OUT)                                :: booltimeseries

!  Internal variables and arrays

INTEGER(I4B)                                                :: id
INTEGER(I4B)                                                :: iff
INTEGER(I4B)                                                :: ids
INTEGER(I4B)                                                :: ls
INTEGER(I4B)                                                :: lzs
INTEGER(I4B)                                                :: nxyz
INTEGER(I4B)                                                :: nlen1
INTEGER(I4B)                                                :: lenformat

nxyz = nx*ny*nz
evapofile = ' '
rate = 0.0d0

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
  
  IF (ssch == 'evaporation') THEN
    id = ids + ls
    CALL sschaine(zone,id,iff,ssch,ids,ls)
    IF(ls /= 0) THEN
        lzs=ls
      !!  *****************************************************************     
      !!check for numerical values
      CALL convan(ssch,lzs,res)
      IF (res == 'n') THEN
      boolfix = .true.
      rate = DNUM(ssch)
      rate = (rate/1000)*dxx(nx)*dzz(nx,ny,nz) ! convert in mm/year in m3/year
      !!  *****************************************************************
      !!check for the name of timeseries file
      ELSEIF (res /= 'n') THEN
        CALL stringtype(ssch,lzs,res)
        evapofile = ssch
        lfile = ls
        booltimeseries = .true.
        id = ids + ls
      CALL sschaine(zone,id,iff,ssch,ids,ls)
        IF(ls /= 0) THEN
        lzs=ls
        CALL convan(ssch,lzs,res)
        IF (res == 'n') THEN
        tslength = int(DNUM(ssch))
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Must provide a numerical value for length time series evapo '
          WRITE(*,*)
          READ(*,*)
          STOP  
        ENDIF
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Must provide a length for time series evapo '
        WRITE(*,*)
        READ(*,*)
        STOP  
        ENDIF

      ENDIF
  !!  *****************************************************************

      ELSE
        WRITE(*,*)
        WRITE(*,*) ' No info provided on evapo '
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      
  ELSE
    GO TO 10
  END IF
  GO TO 10
  
ELSE         ! No string found
  GO TO 10
END IF

1000 RETURN
END SUBROUTINE read_evaporation
