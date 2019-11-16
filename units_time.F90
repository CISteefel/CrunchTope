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

SUBROUTINE units_time(nout,section,time_scale)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
CHARACTER (LEN=mls), INTENT(IN)                                 :: section
REAL(DP), INTENT(OUT)                                           :: time_scale

!  Internal variables and arrays

CHARACTER (LEN=mls)                                             :: parchar
CHARACTER (LEN=mls)                                             :: parfind
CHARACTER (LEN=mls)                                             :: dumstring
CHARACTER (LEN=mls)                                             :: time_units

INTEGER(I4B)                                                    :: lchar

time_scale = 1.0d0
parchar = 'time_units'
parfind = ' '
time_units = ' '
CALL read_string(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN  ! Parameter "time_units" not found
  time_units = 'years'             ! Use default
ELSE
  time_units = dumstring
  
!  Check to see that the units are recognized
  IF (time_units == 'year' .OR. time_units == 'years') THEN
    time_units = 'years'
  END IF
  IF (time_units == 'day' .OR. time_units == 'days') THEN
    time_units = 'days'
  END IF
  IF (time_units == 'hour' .OR. time_units == 'hours') THEN
    time_units = 'hours'
  END IF
  IF (time_units == 'minute' .OR. time_units == 'minutes') THEN
    time_units = 'minutes'
  END IF
  IF (time_units == 'second' .OR. time_units == 'seconds') THEN
    time_units = 'seconds'
  END IF
  
  IF (time_units == 'years') THEN
    time_scale = 1.0d0
  ELSE IF (time_units == 'days') THEN
    time_scale = 1.0d0/(365.0d0)
  ELSE IF (time_units == 'hours') THEN
    time_scale = 1.0d0/(365.0d0*24.0d0)
  ELSE IF (time_units == 'minutes') THEN
    time_scale = 1.0d0/(365.0d0*24.0d0*60.0d0)
  ELSE IF (time_units == 'seconds') THEN
    time_scale = 1.0d0/(365.0d0*24.0d0*60.0d0*60.0d0)
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Time units not recognized'
    WRITE(*,*) ' Section = ',section
    WRITE(*,*) ' Using "years" as time unit'
    WRITE(*,*)
    time_scale = 1.0d0
    READ(*,*)
    STOP
  END IF
END IF

RETURN
END SUBROUTINE units_time
