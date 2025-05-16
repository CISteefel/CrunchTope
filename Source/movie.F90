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

SUBROUTINE movie(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,  &
    ndecay,ikin,nx,ny,nz,realtime,nn,nint,ikmast,ikph,delt,jpor,FirstCall)
USE crunchtype
USE CrunchFunctions
USE params
USE runtime
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE flow
USE temperature
USE strings

IMPLICIT NONE

!  External variables and arrays

REAL(DP), INTENT(IN)                               :: realtime
REAL(DP), INTENT(IN)                               :: delt

INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nrct
INTEGER(I4B), INTENT(IN)                           :: nspec
INTEGER(I4B), INTENT(IN)                           :: ndecay
INTEGER(I4B), INTENT(IN)                           :: nsurf
INTEGER(I4B), INTENT(IN)                           :: nsurf_sec
INTEGER(I4B), INTENT(IN)                           :: ikin
INTEGER(I4B), INTENT(IN)                           :: nkin
INTEGER(I4B), INTENT(IN)                           :: nexchange
INTEGER(I4B), INTENT(IN)                           :: nexch_sec
INTEGER(I4B), INTENT(IN)                           :: nx
INTEGER(I4B), INTENT(IN)                           :: ny
INTEGER(I4B), INTENT(IN)                           :: nz
INTEGER(I4B), INTENT(IN)                           :: nn
INTEGER(I4B), INTENT(IN)                           :: nint
INTEGER(I4B), INTENT(IN)                           :: ikmast
INTEGER(I4B), INTENT(IN)                           :: ikph
INTEGER(I4B), INTENT(IN)                           :: jpor
LOGICAL(LGT), INTENT(IN)                           :: FirstCall

!  Internal variables and arrays

CHARACTER (LEN=mls)                                 :: fn
CHARACTER (LEN=mls)                                  :: suf
CHARACTER (LEN=mls)                                  :: suf1
CHARACTER (LEN=mls)                                 :: fnv
CHARACTER (LEN=1)                                  :: tab
CHARACTER (LEN=mls)                                 :: char_time
CHARACTER (LEN=40)                                 :: prtspecies
 
INTEGER(I4B)                                       :: lspecies
 
INTEGER(I4B), DIMENSION(nrct)                      :: len_min
INTEGER(I4B)                                       :: j
INTEGER(I4B)                                       :: jx
INTEGER(I4B)                                       :: jy
INTEGER(I4B)                                       :: jz
INTEGER(I4B)                                       :: ilength
INTEGER(I4B)                                       :: ik
INTEGER(I4B)                                       :: k
INTEGER(I4B)                                       :: ks
INTEGER(I4B)                                       :: ns
INTEGER(I4B)                                       :: i
INTEGER(I4B)                                       :: nex
INTEGER(I4B)                                       :: ir
INTEGER(I4B)                                       :: lsjx
INTEGER(I4B)                                       :: ls
INTEGER(I4B)                                       :: nlen

CHARACTER (LEN=mls)                                        :: namtemp

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: is

REAL(DP)                                                   :: phwrite
REAL(DP)                                                   :: phprt


jz = 1

suf='.out'
suf1 ='.out'
tab = CHAR(9)

!  Write out master variable

OPEN(UNIT=15,STATUS='scratch')
WRITE(15,*) realtime
REWIND 15
READ(15,'(a)') char_time
CLOSE(UNIT=15)

CALL stringlen(char_time,ls)

200 FORMAT(1PE9.2)

IF (FirstCall) THEN
  IF (ikTracer /= 0) THEN
    OPEN(UNIT=89,FILE='TracerMovie.out', ACCESS='sequential',STATUS='unknown')
    WRITE(89,1013)
    WRITE(89,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(89,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,s(ikTracer,jx,jy,jz)
      END DO
    END DO
  END IF
ELSE
  IF (ikTracer /= 0) THEN
    WRITE(89,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(89,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,s(ikTracer,jx,jy,jz)
      END DO
    END DO
  END IF
END IF

IF (FirstCall) THEN
  IF (ikpH /= 0) THEN
    OPEN(UNIT=90,FILE='pHMovie.out', ACCESS='sequential',STATUS='unknown')
    WRITE(90,1014)
    WRITE(90,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        phprt =  -(sp(ikph,jx,jy,jz)+lngamma(ikph,jx,jy,jz))/clg
        WRITE(90,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,phprt
      END DO
    END DO
  END IF
ELSE
  IF (ikph /= 0) THEN
    WRITE(90,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        phprt =  -(sp(ikph,jx,jy,jz)+lngamma(ikph,jx,jy,jz))/clg
        WRITE(90,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,phprt
      END DO
    END DO
  END IF
END IF

!!!    IF (FirstCall) THEN
!!!      
!!!      OPEN(UNIT=91,FILE='HaliteMovie.out', ACCESS='sequential',STATUS='unknown')
!!!      WRITE(91,1019)
!!!      WRITE(91,*) 'ZONE F=POINT,I=', nx,  ', J=',ny    
!!!
!!!        DO jy = 1,ny
!!!          DO jx = 1,nx
!!!            WRITE(91,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,volfx(1,jx,jy,1)
!!!          END DO
!!!        END DO
!!!
!!!    ELSE
!!!      WRITE(91,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
!!!        DO jy = 1,ny
!!!          DO jx = 1,nx
!!!            WRITE(91,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,volfx(1,jx,jy,1)
!!!          END DO
!!!        END DO
!!!
!!!    END IF
    
!!!    IF (FirstCall) THEN
!!!      
!!!      OPEN(UNIT=92,FILE='HaliteRateMovie.out', ACCESS='sequential',STATUS='unknown')
!!!      WRITE(92,1019)
!!!      WRITE(92,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
!!!    
!!!        DO jy = 1,ny
!!!          DO jx = 1,nx
!!!            WRITE(92,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,dppt(1,jx,jy,1)
!!!          END DO
!!!        END DO
!!!
!!!    ELSE
!!!      WRITE(92,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
!!!        DO jy = 1,ny
!!!          DO jx = 1,nx
!!!            WRITE(92,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,dppt(1,jx,jy,1)
!!!          END DO
!!!        END DO
!!!
!!!    END IF
    
    IF (FirstCall) THEN
      
      OPEN(UNIT=93,FILE='PorosityMovie.out', ACCESS='sequential',STATUS='unknown')
      WRITE(93,1020)
      WRITE(93,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    

        DO jy = 1,ny
          DO jx = 1,nx
            WRITE(93,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,por(jx,jy,1)
          END DO
        END DO

    ELSE
      
      WRITE(93,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
        DO jy = 1,ny
          DO jx = 1,nx
            WRITE(93,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,por(jx,jy,1)
          END DO
        END DO

    END IF
  
  IF (FirstCall) THEN
    
    OPEN(UNIT=94,FILE='H2Evolve.out', ACCESS='sequential',STATUS='unknown')
    WRITE(94,1021)
    WRITE(94,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
      
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(94,184)  x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, spgas10(1,jx,jy,1)
      END DO
    END DO

  ELSE
    
    WRITE(94,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(94,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, spgas10(1,jx,jy,1)
      END DO
    END DO

  END IF


185 FORMAT(1PE12.5,12x,100(1X,1PE13.5))
1010 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "Porosity"')
1011 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "Permeability"')
1013 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "Tracer"')
1014 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "pH"')
1015 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "Fe(II)"')
1016 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "Acetate"')
1017 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "Goethite"')
1018 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "Calcite"')
1019 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "Halite"')
1020 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "Porosity"')
1021 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "H2(g)"')

1009 FORMAT('VARIABLES = " X (meters)", "  Y (meters)  "',100(', "',A10,'"'))
2009 FORMAT('VARIABLES = " X (meters)", "  Z (meters)  "',100(', "',A10,'"'))
2001 FORMAT('VARIABLES = "X (meters)"',                   100(', "',A10,'"'))

1012 FORMAT('VARIABLES = " X (meters)", " Y (meters)", "X Velocity", "Y Velocity"')
2012 FORMAT('VARIABLES = " X (meters)", " Z (meters)", "X Velocity", "Z Velocity"')
!!1013 FORMAT('VARIABLES = " X (meters)", " Y (meters)", "Pressure"')
2014 FORMAT('VARIABLES = " X (meters)", " Y (meters)", "Tortuosity"')

182 FORMAT(100(1X,1PE12.4))
183 FORMAT(1PE12.4,2X,1PE12.4,2X,1PE12.4)
184 FORMAT(100(1X,1PE15.7))

2283 FORMAT('# Time (yrs) ',2X,1PE12.4)
2284 FORMAT('      X        ','     Y        ',a18)
2282 FORMAT('   X        ','     Y        ','        pH')
2281 FORMAT('   X        ','     Y        ',4X,a18)
2285 FORMAT('    X        ','     Y        ',3X,30(1X,a13))

RETURN
END SUBROUTINE movie
!  *******************************************************
