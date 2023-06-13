!---------------------------------------------------------------------------
! Subroutine writing the header of output files in NVIEW format
!
! original:      Yunwei Sun
! last modified: Olivier Bildstein 11/13/00
!---------------------------------------------------------------------------  
SUBROUTINE xtoolInit(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,ltitle,  &
    nsurf_sec,ndecay,ikin,nx,ny,nz,realtime,nn,nint,ikmast,ikph,ikO2,master,delt)
USE crunchtype
USE params
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

CHARACTER (LEN=2*mls)                              :: ltitle

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
INTEGER(I4B), INTENT(IN)                           :: ikO2

CHARACTER (LEN=mls), INTENT(IN)                    :: master

!  Internal variables and arrays

CHARACTER (LEN=13), DIMENSION(nrct)                :: uminprnt
CHARACTER (LEN=13), DIMENSION(ncomp+nspec)         :: ulabprnt
CHARACTER (LEN=20)                                 :: fn
CHARACTER (LEN=4)                                  :: suf
CHARACTER (LEN=4)                                  :: suf1
CHARACTER (LEN=20)                                 :: fnv
CHARACTER (LEN=1)                                  :: tab
CHARACTER (LEN=mls), DIMENSION(nsurf+nsurf_sec)    :: prtsurf
CHARACTER (LEN=mls)                                :: tempstring
CHARACTER (LEN=mls)                                :: tempstring2
 
INTEGER(I4B), DIMENSION(ncomp+nspec)               :: len_sp
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
INTEGER(I4B)                                       :: nlen
INTEGER(I4B)                                       :: noutput
INTEGER(I4B)                                       :: nexcomp
INTEGER(I4B)                                       :: nsurfcomp

REAL(DP), DIMENSION(ncomp)                         :: totex_bas
REAL(DP), DIMENSION(nrct)                          :: dptprt
REAL(DP), DIMENSION(nrct)                          :: dsat
REAL(DP), DIMENSION(nrct)                          :: dvolpr
REAL(DP)                                           :: sum
REAL(DP)                                           :: porprt
REAL(DP)                                           :: phprt
REAL(DP)                                           :: porcalc
REAL(DP)                                           :: dzzz
REAL(DP)                                           :: xcumulative


ALLOCATE(exflag(ncomp))
DO i= 1, ncomp
   exflag(i)=.FALSE.
   DO nex=1,nexch_sec
      IF (muexc(nex,i) /= 0.) THEN
         exflag(i)=.TRUE.
      ENDIF
   END DO 
END DO
nexcomp = 0
DO i= 1, ncomp
   IF (exflag(i)) THEN
      nexcomp = nexcomp + 1
   ENDIF
END DO

nsurfcomp = 0
ALLOCATE(surflag(ncomp))
surflag = .FALSE.
IF (nsurf > 0) THEN
   DO i= 1, ncomp
     DO ns=1,nsurf_sec
       IF (musurf(ns,i) /= 0.) THEN
         surflag(i)=.TRUE.
       ENDIF
     END DO 
   END DO
   DO i= 1, ncomp
     IF (surflag(i)) THEN
       nsurfcomp = nsurfcomp + 1
     ENDIF
   END DO
ENDIF

DO ks = 1,nsurf
  prtsurf(ks) = namsurf(ks)
END DO
DO ns = 1,nsurf_sec
  prtsurf(ns+nsurf) = namsurf_sec(ns)
END DO

432   FORMAT(e14.7)
!!432   FORMAT(1pe13.7)
465   FORMAT('      ',I5)

!         files header (Common)
!
OPEN(32,FILE='crunch-solution.ext',STATUS='UNKNOWN')
WRITE(32,401) 
401   FORMAT('         1')
WRITE(32,406) ltitle
406   FORMAT(a132)
WRITE(32,402) 
402   FORMAT('         10')
WRITE(32,400)
400   FORMAT('internal')
WRITE(32,403) 
403   FORMAT('Output for Aqueous phase Properties')
WRITE(32,404) 
404   FORMAT('0')
WRITE(32,405)  
405   FORMAT('VARIABLE')

WRITE(32,407) 
407   FORMAT('$gdef')
WRITE(32,408) 
408   FORMAT('$type rect')

WRITE(32,409) nx
WRITE(32,410) ny
WRITE(32,411) nz
WRITE(32,412) 


409   FORMAT('$nx ',I5)
410   FORMAT('$ny ',I5)
411   FORMAT('$nz ',I5)
412   FORMAT('$order yzx')

WRITE(32,421) 
xcumulative = 0.0
DO jx=1, nx
   WRITE(32,432) dxx(jx)
   xcumulative = xcumulative + dxx(jx)
END DO

IF (ny == 1 .AND. nz == 1) THEN
  WRITE(32,422) 
  dzzz=1.0
  WRITE(32,432) dzzz
  WRITE(32,423) 
  xcumulative = xcumulative/3.0    !  Use an aspect ratio of 3 to 1 for 1D cases
  WRITE(32,432) xcumulative
ELSE
  WRITE(32,422) 
  DO jy=1,ny
    WRITE(32,432) dyy(jy)
  END DO
  WRITE(32,423) 
  DO jz=1,nz
    WRITE(32,432) dzz(1,1,jz)
  END DO
END IF

421 FORMAT('$dx ')
422   FORMAT('$dy ')
423   FORMAT('$dz ')


WRITE(32,424) 
424   FORMAT('$end_internal_grid')
WRITE(32,425) 
425   FORMAT('$OperatingSystem ')
WRITE(32,426) 
426   FORMAT('$C-Compiler  ')
WRITE(32,427) 
427   FORMAT('$FortranCompiler  ')
WRITE(32,428) 
428   FORMAT('$RunID  0')
WRITE(32,429) 
429   FORMAT('$RunDate ')

! number of output variables 
! pH + O2(aq) + total_conc + exchange + totexchange +
! surface + totsurface
!!noutput = 2*ncomp+nspec+nexch_sec+nexcomp+nsurf+nsurf_sec+nsurfcomp
noutput = ncomp+nexch_sec+nexcomp+nsurf+nsurf_sec+nsurfcomp
IF (ikph /= 0) THEN
   noutput = noutput + 1
ENDIF ! for pH
IF (ikO2 /= 0) THEN
   noutput = noutput + 1
ENDIF ! for O2(aq)
WRITE(32,465) noutput 

230 FORMAT(a14)
IF (ikph /= 0) THEN
   WRITE(32,230) "pH            "
ENDIF ! for pH
IF (ikO2 /= 0) THEN
   WRITE(32,230) "O2(aq)        "
ENDIF ! for O2(aq)
240 FORMAT("TOT-",a14)
DO i=1,ncomp
   WRITE(32,240) ulab(i)
END DO
!!DO i=1,ncomp + nspec
!!   WRITE(32,230) ulab(i)
!!END DO
DO nex=1,nexch_sec
   WRITE(32,230) nam_exchsec(nex)
END DO
241 FORMAT("TOTEXCH-",a14)
DO i=1,ncomp
   IF (exflag(i)) then
      WRITE(32,241) ulab(i)
   ENDIF
END DO
DO ns=1,nsurf+nsurf_sec
   WRITE(32,230) prtsurf(ns)
END DO
242 FORMAT("TOTSURF-",a14)
DO i=1,ncomp
   IF (surflag(i)) then
     WRITE(32,242) ulab(i)
   ENDIF
END DO

WRITE(32,430) 
430   FORMAT('        0')
WRITE(32,465) nx*ny*nz 

DO jx=1,nx
   DO jy=1,ny
      DO jz=1,nz
!!         WRITE(32, *) "rock#",jx,":",jy,":",jz
         WRITE(32, 1001) jx,jz,jy
      END DO
   END DO
END DO

1001 FORMAT("rock#",i3,":",i3,":",i3)
231 FORMAT("THIST")
232 FORMAT("EVAR ",a14)
233 FORMAT("T *")
234 FORMAT("E *")
235 FORMAT(" ")
WRITE(32,430) 
IF (ikph /= 0) THEN
   WRITE(32,231)
   WRITE(32,232) "pH            "
   WRITE(32,233)
   WRITE(32,234)
   WRITE(32,235)
ENDIF ! for pH
IF (ikO2 /= 0) THEN
   WRITE(32,231)
   WRITE(32,232) "O2(aq)        "
   WRITE(32,233)
   WRITE(32,234)
   WRITE(32,235)
ENDIF ! for O2(aq)
250 FORMAT("EVAR TOT-",a14)
DO i=1,ncomp
   WRITE(32,231)
   WRITE(32,250) ulab(i)
   WRITE(32,233)
   WRITE(32,234)
   WRITE(32,235)
END DO
!!DO i=1,ncomp+nspec
!!   WRITE(32,231)
!!   WRITE(32,232) ulab(i)
!!   WRITE(32,233)
!!   WRITE(32,234)
!!   WRITE(32,235)
!!END DO
DO nex=1,nexch_sec
   WRITE(32,231)
   WRITE(32,232) nam_exchsec(nex)
   WRITE(32,233)
   WRITE(32,234)
   WRITE(32,235)
END DO
251 FORMAT("EVAR TOTEXCH-",a14)
DO i=1,ncomp
   IF (exflag(i)) then
      WRITE(32,231)
      WRITE(32,251) ulab(i)
      WRITE(32,233)
      WRITE(32,234)
      WRITE(32,235)
   ENDIF
END DO
DO ns=1,nsurf+nsurf_sec
   WRITE(32,231)
   WRITE(32,232) prtsurf(ns)
   WRITE(32,233)
   WRITE(32,234)
   WRITE(32,235)
END DO
252 FORMAT("EVAR TOTSURF-",a14)
DO i=1,ncomp
   IF (surflag(i)) then
      WRITE(32,231)
      WRITE(32,252) ulab(i)
      WRITE(32,233)
      WRITE(32,234)
      WRITE(32,235)
   ENDIF
END DO

!--------------------------------------------------!
OPEN(33,FILE='crunch-mineral.ext',STATUS='UNKNOWN')
WRITE(33,401) 
WRITE(33,406) ltitle
WRITE(33,402) 
WRITE(33,400)
WRITE(33,431) 
431   FORMAT('Output for Mineral Properties')
WRITE(33,404) 
WRITE(33,405)  

WRITE(33,407) 
WRITE(33,408) 
WRITE(33,409) nx
WRITE(33,410) nz
WRITE(33,411) ny
WRITE(33,412) 

WRITE(33,421) 
DO jx=1, nx
   WRITE(33,432) dxx(jx)
END DO

IF (ny == 1 .AND. nz == 1) THEN
  WRITE(33,422) 
  dzzz=1.0
  WRITE(33,432) dzzz
  WRITE(33,423) 
  WRITE(33,432) xcumulative
ELSE
  WRITE(33,422) 
  dzzz=1
  DO jz=1, nz
    WRITE(33,432) dzzz
  END DO
  WRITE(33,423) 
  DO jy=1, ny
    WRITE(33,432) dyy(jy)
  END DO
END IF

WRITE(33,424) 
WRITE(33,425) 
WRITE(33,426) 
WRITE(33,427) 
WRITE(33,428) 
WRITE(33,429) 

! number of output variables (volume+saturation+rate+porosity)
noutput = 3*nrct 

IF (nrct > 1) THEN
   noutput = noutput + 1
ENDIF ! for porosity
WRITE(33,465) noutput 

IF (nrct > 1) THEN
   WRITE(33,230) "Porosity      "
ENDIF ! for porosity
DO k=1,nrct
   WRITE(33,230) umin(k)
END DO ! for volume
260 FORMAT("RATE-",a14)
DO k=1,nrct
   WRITE(33,260) umin(k)
END DO ! for rates
261 FORMAT("SAT-",a14)
DO k=1,nrct
   WRITE(33,261) umin(k)
END DO ! for saturation

WRITE(33,430) 
WRITE(33,465) nx*ny*nz 

DO jx=1,nx
   DO jy=1,ny
      DO jz=1,nz
         WRITE(33, 1001) jx,jz,jy
      END DO
   END DO
END DO

WRITE(33,430) 

IF (nrct > 1) THEN
   WRITE(33,231)
   WRITE(33,232) "Porosity      "
   WRITE(33,233)
   WRITE(33,234)
   WRITE(33,235)
ENDIF ! for porosity
!!270 FORMAT("EVAR VOL-",a14)
DO k=1,nrct
   WRITE(33,231)
   WRITE(33,232) umin(k)
   WRITE(33,233)
   WRITE(33,234)
   WRITE(33,235)
END DO
271 FORMAT("EVAR RATE-",a14)
DO k=1,nrct
   WRITE(33,231)
   WRITE(33,271) umin(k)
   WRITE(33,233)
   WRITE(33,234)
   WRITE(33,235)
END DO
272 FORMAT("EVAR SAT-",a14)
DO k=1,nrct
   WRITE(33,231)
   WRITE(33,272) umin(k)
   WRITE(33,233)
   WRITE(33,234)
   WRITE(33,235)
END DO

RETURN
END SUBROUTINE xtoolInit

