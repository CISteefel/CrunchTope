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

MODULE CrunchFunctions

CONTAINS

#ifdef __GNUC__

!!  *******************************************
  PURE FUNCTION jnum (string)
  USE CrunchType
  USE params
  IMPLICIT NONE
  CHARACTER (LEN=mls), INTENT(IN)          :: string
  INTEGER(I4B)                           :: jnum
  READ(string, *) jnum
  END FUNCTION jnum
!! **********************************************
  PURE FUNCTION dnum (string)
  USE CrunchType
  USE params
  IMPLICIT NONE
  CHARACTER (LEN=mls), INTENT(IN)          :: string
  REAL(DP)                               :: dnum
  READ(string, *) dnum
  END FUNCTION dnum
!! **********************************************
  PURE FUNCTION cosd (degrees)
  USE CrunchType
  IMPLICIT NONE
  REAL(DP),INTENT(IN)                  :: degrees
  REAL(DP)                             :: cosd
  real(dp)                             :: pi
  real(dp)                             :: rad
  pi = dacos(-1.0d0)
  rad = degrees * pi / 180.0d0
  cosd = dcos(rad)
  END FUNCTION cosd
!! **********************************************
#endif

  PURE FUNCTION IntegerToCharacter (i)
  USE CrunchType
  USE params
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)               :: i
  CHARACTER (LEN=mls)                    :: IntegerToCharacter
  CHARACTER (LEN=mls)                    :: dumstring
  WRITE(dumstring,*) i
  IntegerToCharacter = ADJUSTL(dumstring)
  END FUNCTION IntegerToCharacter

!! **********************************************
  PURE FUNCTION RealToCharacter (x)
  USE CrunchType
  USE params
  IMPLICIT NONE
  REAL(DP), INTENT(IN)                   :: X
  CHARACTER (LEN=mls)                    :: RealToCharacter
  CHARACTER (LEN=mls)                    :: dumstring
  WRITE(dumstring,*) x
  RealToCharacter = ADJUSTL(dumstring)
  END FUNCTION RealToCharacter
!! **********************************************
PURE FUNCTION ArithmeticMean (value1,value2)
USE crunchtype
USE params
IMPLICIT NONE
REAL(DP),INTENT(IN)                           :: value1
REAL(DP),INTENT(IN)                           :: value2
REAL(DP)                                      :: ArithmeticMean
ArithmeticMean = 0.5d0*(value1+value2)
END FUNCTION ArithmeticMean
!! **********************************************
PURE FUNCTION HarmonicMean (value1,value2)
USE crunchtype
USE params
IMPLICIT NONE
REAL(DP),INTENT(IN)                           :: value1
REAL(DP),INTENT(IN)                           :: value2
REAL(DP)                                      :: HarmonicMean
REAL(DP)                                      :: denominator
denominator = value1 + value2
IF (denominator /= 0.0d0) THEN
  HarmonicMean = 2.0d0*value1*value2/denominator
ELSE
  HarmonicMean = 0.0d0
END IF
END FUNCTION HarmonicMean
!! **********************************************
PURE FUNCTION GeometricMean (value1,value2)
USE crunchtype
USE params
IMPLICIT NONE
REAL(DP),INTENT(IN)                           :: value1
REAL(DP),INTENT(IN)                           :: value2
REAL(DP)                                      :: GeometricMean
GeometricMean = DSQRT(value1*value2)
END FUNCTION GeometricMean
!! **********************************************
PURE FUNCTION GetPrimarySpeciesNumber (ncomp,dumstring)
USE CrunchType
USE params, ONLY: mls
USE concentration,ONLY: ulab
IMPLICIT NONE
INTEGER(I4B), INTENT(IN)               :: ncomp
CHARACTER (LEN=mls), INTENT(IN)        :: dumstring
INTEGER(I4B)                           :: GetPrimarySpeciesNumber
INTEGER(I4B)                           :: i
INTEGER(I4B)                           :: ls
GetPrimarySpeciesNumber = 0
DO i = 1,ncomp
  IF (dumstring == ulab(i)) THEN
    GetPrimarySpeciesNumber = i
  END IF
END DO
END FUNCTION GetPrimarySpeciesNumber

PURE FUNCTION imaxloc(arr)
USE crunchtype
REAL(DP), DIMENSION(:), INTENT(IN) :: arr
INTEGER(I4B) :: imaxloc
INTEGER(I4B), DIMENSION(1) :: imax

imax=maxloc(arr(:))
imaxloc=imax(1)
END FUNCTION imaxloc

PURE FUNCTION outerprod(a,b)
USE crunchtype
REAL(DP), DIMENSION(:), intent(in)                         :: a,b
REAL(DP), DIMENSION(size(a),size(b))                       :: outerprod

outerprod = spread(a,dim=2,ncopies=size(b)) * spread(b,dim=1,ncopies=size(a))
END FUNCTION outerprod

END MODULE CrunchFunctions