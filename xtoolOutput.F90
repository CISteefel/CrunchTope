!---------------------------------------------------------------------------
! Subroutine writing graphics output files in XPLOT format
!
! original:      Yunwei Sun
! last modified: Olivier Bildstein 10/30/00
!---------------------------------------------------------------------------  
SUBROUTINE xtoolOutput(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,  &
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
INTEGER(I4B)                                       :: nxyz
INTEGER(I4B)                                       :: noutput

REAL(DP)                                           :: dptprt
REAL(DP), DIMENSION(nrct)                          :: dsat
REAL(DP), DIMENSION(nrct)                          :: dvolpr
REAL(DP)                                           :: sum
REAL(DP)                                           :: porprt
REAL(DP)                                           :: phprt
REAL(DP)                                           :: porcalc
REAL(DP)                                           :: times
REAL(DP)                                           :: totex_bas
       
nxyz = nx*ny*nz
times = realtime*31557600.

432   FORMAT(e14.7)
465   FORMAT('       ',I5)

noutput = 0 
! pH
IF (ikph /= 0) THEN
   noutput = noutput + 1
   WRITE(32,432) times
   WRITE(32,465) noutput
   WRITE(32,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            WRITE(32,432) -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
         END DO
      END DO
   END DO
ENDIF   ! END OF pH

! O2
IF (ikO2 /= 0) THEN
   noutput = noutput + 1
   WRITE(32,432) times
   WRITE(32,465) noutput
   WRITE(32,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            WRITE(32,432) (sp(ikO2,jx,jy,jz)+gam(ikO2,jx,jy,jz))/clg
         END DO
      END DO
   END DO
ENDIF   ! END OF O2

! Total concentration
DO i=1,ncomp
   noutput = noutput + 1
   WRITE(32,432) times
   WRITE(32,465) noutput
   WRITE(32,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            WRITE(32,432) s(i,jx,jy,jz)
         END DO
      END DO
   END DO
END DO   ! END OF ncomp

! Log species concentration
!!DO i=1,ncomp + nspec
!!   noutput = noutput + 1
!!   WRITE(32,432) times
!!   WRITE(32,465) noutput
!!   WRITE(32,465) nxyz
!!   DO jx=1,nx
!!      DO jy=1,ny
!!         DO jz=1,nz
!!            WRITE(32,432) sp(i,jx,jy,jz)/clg
!!         END DO
!!      END DO
!!   END DO
!!END DO   ! END OF ncomp

! Adsorbates
DO nex=1,nexch_sec
   noutput = noutput + 1
   WRITE(32,432) times
   WRITE(32,465) noutput
   WRITE(32,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            WRITE(32,432) spex10(nex+nexchange,jx,jy,jz)
         END DO
      END DO
   END DO
END DO   ! END OF nexch_sec

! Total adsorbates
DO i=1,ncomp
   IF (exflag(i)) then
      noutput = noutput + 1
      WRITE(32,432) times
      WRITE(32,465) noutput
      WRITE(32,465) nxyz
      DO jx=1,nx
        DO jy=1,ny
          DO jz=1,nz
            totex_bas = 0.0
            DO nex = 1,nexch_sec
            totex_bas = totex_bas + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
            END DO
            WRITE(32,432) totex_bas
          END DO
        END DO
      END DO
    ENDIF  
END DO   ! END OF ncomp

! Surface complexes
DO ns=1,nsurf+nsurf_sec
   noutput = noutput + 1
   WRITE(32,432) times
   WRITE(32,465) noutput
   WRITE(32,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            WRITE(32,432) spsurf10(ns,jx,jy,jz)
         END DO
      END DO
   END DO
END DO   ! END OF nexch_sec

! Total surface complexes
IF (nsurf > 0) THEN
  DO i=1,ncomp
    IF (surflag(i)) then
      noutput = noutput + 1
      WRITE(32,432) times
      WRITE(32,465) noutput
      WRITE(32,465) nxyz
      DO jx=1,nx
        DO jy=1,ny
          DO jz=1,nz
             totex_bas = 0.0
             DO ns = 1,nsurf_sec
                totex_bas = totex_bas + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
             END DO            
             WRITE(32,432) totex_bas
          END DO
        END DO
      END DO
    ENDIF
  END DO   ! END OF ncomp
ENDIF

!------------------------!
noutput = 0

! Porosity
IF (nrct > 1) THEN
   noutput = noutput + 1
   WRITE(33,432) times
   WRITE(33,*) noutput
   WRITE(33,*) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            sum = 0.0
            DO k = 1,nrct
              sum = sum + volfx(k,jx,jy,jz)
            END DO
            WRITE(33,432) (1.0-sum)*100.0
         END DO
      END DO
   END DO
ENDIF

! Mineral volume fraction
DO k=1,nrct
   noutput = noutput + 1
   WRITE(33,432) times
   WRITE(33,465) noutput
   WRITE(33,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            WRITE(33,432) volfx(k,jx,jy,jz)*100.0
         END DO
      END DO
   END DO
END DO   ! END OF nrct

! Mineral rate
DO k=1,nrct
   noutput = noutput + 1
   WRITE(33,432) times
   WRITE(33,465) noutput
   WRITE(33,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
!************************
!  For units of volume %/year, uncomment the following line and
!  recompile
!            dptprt = dppt(k,jx,jy,jz)*volmol(k)*100.0  ! volume %/yr
!***********************
!************************
!  For units of mol/L(BV)/sec:
            dptprt = dppt(k,jx,jy,jz)/(secyr*1000.0)    ! mol/L(BV)/sec
!*************************************
            WRITE(33,432) dptprt
         END DO
      END DO
   END DO
END DO   ! END OF nrct

! Mineral saturation
DO k=1,nrct
   noutput = noutput + 1
   WRITE(33,432) times
   WRITE(33,465) noutput
   WRITE(33,465) nxyz
   DO jx=1,nx
      DO jy=1,ny
         DO jz=1,nz
            CALL reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt,realtime)
            CALL satcalc(ncomp,nrct,jx,jy,jz)
            WRITE(33,432) silog(1,k)
         END DO
      END DO
   END DO
END DO   ! END OF nrct



RETURN
END SUBROUTINE xtoolOutput

