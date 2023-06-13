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

SUBROUTINE units_distance(nout,section,dist_scale)
USE crunchtype
USE params
USE strings

IMPLICIT NONE

!  External variables and arrays

INTEGER(I4B), INTENT(IN)                                        :: nout
CHARACTER (LEN=mls), INTENT(IN)                                 :: section
REAL(DP), INTENT(OUT)                                           :: dist_scale

!  Internal variables and arrays

CHARACTER (LEN=mls)                                             :: parchar
CHARACTER (LEN=mls)                                             :: parfind
CHARACTER (LEN=mls)                                             :: dumstring
CHARACTER (LEN=mls)                                             :: distance_units

INTEGER(I4B)                                                    :: lchar

!  Search for time units (default is years)

dist_scale = 1.0
parchar = 'distance_units'
parfind = ' '
distance_units = ' '
CALL read_string(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == 'distance_units' .OR. parfind == 'distance_unit') THEN  
  distance_units = dumstring
  
!  Check to see that the units are recognized
  IF (distance_units == 'meter' .OR. distance_units == 'meters' .OR. distance_units == 'm') THEN
    distance_units = 'meters'
  END IF
  IF (distance_units == 'centimeter' .OR. distance_units == 'centimeters'.OR. distance_units == 'cm') THEN
    distance_units = 'centimeters'
  END IF
  IF (distance_units == 'millimeter' .OR. distance_units == 'millimeters'.OR. distance_units == 'mm') THEN
    distance_units = 'millimeters'
  END IF
  IF (distance_units == 'kilometer' .OR. distance_units == 'kilometers'.OR. distance_units == 'km') THEN
    distance_units = 'kilometers'
  END IF
  IF (distance_units == 'micrometer' .OR. distance_units == 'micrometers'.OR.  distance_units == 'um') THEN
    distance_units = 'micrometers'
  END IF
  IF (distance_units == 'micron' .OR. distance_units == 'microns') THEN
    distance_units = 'micrometers'
  END IF
  IF (distance_units == 'nanometer' .OR. distance_units == 'nanometers'.OR. distance_units == 'nm') THEN
    distance_units = 'nanometers'
  END IF
  IF (distance_units == 'm') THEN
    distance_units = 'meters'
  END IF
  IF (distance_units == 'cm') THEN
    distance_units = 'centimeters'
  END IF
  IF (distance_units == 'mm') THEN
    distance_units = 'millimeters'
  END IF
  IF (distance_units == 'km') THEN
    distance_units = 'kilometers'
  END IF
  IF (distance_units == 'um') THEN
    distance_units = 'micrometers'
  END IF
  IF (distance_units == 'nm') THEN
    distance_units = 'nanometers'
  END IF
  
  IF (distance_units == 'meters') THEN
    dist_scale = 1.0
  ELSE IF (distance_units == 'centimeters') THEN
    dist_scale = 100.0
  ELSE IF (distance_units == 'millimeters') THEN
    dist_scale = 1000.0
  ELSE IF (distance_units == 'kilometers') THEN
    dist_scale = 0.001
  ELSE IF (distance_units == 'micrometers') THEN
    dist_scale = 1.0E06
  ELSE IF (distance_units == 'nanometers') THEN
    dist_scale = 1.0E09
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Distance units not recognized'
    WRITE(*,*) ' Section = ',section
    WRITE(*,*) ' Using "meters" as distance unit'
    WRITE(*,*)
    dist_scale = 1.0
    READ(*,*)
    STOP
  END IF

ELSE
  distance_units = 'meters'             ! Use default
  dist_scale = 1.0
END IF

RETURN
END SUBROUTINE units_distance
