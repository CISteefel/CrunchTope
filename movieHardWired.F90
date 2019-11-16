!************** (C) COPYRIGHT 1993 Carl I. Steefel *******************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:43:15
 
!                      All Rights Reserved

!  GIMRT IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR
!  IMPLIED. THE USER ASSUMES ALL RISKS OF USING 1DREACT. THERE  IS
!  NO CLAIM OF THE MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO USERS AT
!  ANY SITES OTHER THAN YOUR OWN.
!**********************************************************************

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

200 FORMAT(1PE8.2)

IF (FirstCall) THEN
  OPEN(UNIT=88,FILE='PorosityMovie.out', ACCESS='sequential',STATUS='unknown')
  WRITE(88,1010)
  WRITE(88,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jy = 1,ny
      DO jx = 1,nx
        WRITE(88,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,por(jx,jy,jz)
      END DO
  END DO
ELSE
  WRITE(88,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jy = 1,ny
      DO jx = 1,nx
        WRITE(88,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,por(jx,jy,jz)
      END DO
  END DO
END IF

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
    OPEN(UNIT=91,FILE='pHMovie.out', ACCESS='sequential',STATUS='unknown')
    WRITE(91,1014)
    WRITE(91,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        phwrite = -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz) )/clg
        WRITE(91,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,phwrite
      END DO
    END DO
  END IF
ELSE
  IF (ikpH /= 0) THEN
    WRITE(91,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        phwrite = -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz) )/clg
        WRITE(91,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,phwrite
      END DO
    END DO
  END IF
END IF



IF (FirstCall) THEN
  OPEN(UNIT=90,FILE='PermMovie.out', ACCESS='sequential',STATUS='unknown')
  WRITE(90,1011)
  WRITE(90,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jy = 1,ny
      DO jx = 1,nx
        WRITE(90,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,permx(jx,jy,jz)
      END DO
  END DO
ELSE
  WRITE(90,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jy = 1,ny
      DO jx = 1,nx
        WRITE(90,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,permx(jx,jy,jz)
      END DO
  END DO
END IF

IF (FirstCall) THEN
  OPEN(UNIT=92,FILE='VelocityMovie.out', ACCESS='sequential',STATUS='unknown')
  WRITE(92,1012)
  WRITE(92,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jy = 1,ny
      DO jx = 1,nx
        WRITE(92,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,qx(jx,jy,jz),qy(jx,jy,jz)
      END DO
  END DO
ELSE
  WRITE(92,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jy = 1,ny
      DO jx = 1,nx
        WRITE(92,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,qx(jx,jy,jz),qy(jx,jy,jz)
      END DO
  END DO
END IF

!!  *************  Hardwired numbers  ************************************

!!  For Fe(II)

IF (FirstCall) THEN
  IF (ncomp >= 2) THEN
    OPEN(UNIT=93,FILE='FeIIMovie.out', ACCESS='sequential',STATUS='unknown')
    WRITE(93,1015)
    WRITE(93,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(93,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,s(2,jx,jy,jz)
      END DO
    END DO
  END IF
ELSE
  IF (ncomp >= 2) THEN
    WRITE(93,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(93,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,s(2,jx,jy,jz)
      END DO
    END DO
  END IF
END IF

!!  For Acetate

IF (FirstCall) THEN
  IF (ncomp >= 8) THEN
    OPEN(UNIT=94,FILE='AcetateMovie.out', ACCESS='sequential',STATUS='unknown')
    WRITE(94,1016)
    WRITE(94,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(94,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,s(8,jx,jy,jz)
      END DO
    END DO
  END IF
ELSE
  IF (ncomp >= 8) THEN
    WRITE(94,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(94,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,s(8,jx,jy,jz)
      END DO
    END DO
  END IF
END IF

!!  For Goethite

IF (FirstCall) THEN
  IF (nrct >= 1) THEN
    OPEN(UNIT=95,FILE='GoethiteMovie.out', ACCESS='sequential',STATUS='unknown')
    WRITE(95,1017)
    WRITE(95,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(95,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,volfx(1,jx,jy,jz)
      END DO
    END DO
  END IF
ELSE
  IF (nrct >= 1) THEN
    WRITE(95,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(95,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,volfx(1,jx,jy,jz)
      END DO
    END DO
  END IF
END IF

!!  For Calcite

IF (FirstCall) THEN
  IF (nrct >= 3) THEN
    OPEN(UNIT=96,FILE='CalciteMovie.out', ACCESS='sequential',STATUS='unknown')
    WRITE(96,1018)
    WRITE(96,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(96,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,volfx(3,jx,jy,jz)
      END DO
    END DO
  END IF
ELSE
  IF (nrct >= 3) THEN
    WRITE(96,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(96,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,volfx(3,jx,jy,jz)
      END DO
    END DO
  END IF
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
