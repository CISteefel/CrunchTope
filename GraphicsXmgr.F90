!************** (C) COPYRIGHT 1993 Carl I. Steefel *******************
 
! Code converted using TO_F90 by Alan Miller
! Date: 2000-07-27  Time: 10:43:21
 
!                      All Rights Reserved

!  GIMRT IS PROVIDED "AS IS" AND WITHOUT ANY WARRANTY EXPRESS OR
!  IMPLIED. THE USER ASSUMES ALL RISKS OF USING 1DREACT. THERE  IS
!  NO CLAIM OF THE MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

!  YOU MAY MODIFY THE SOURCE CODE FOR YOUR OWN USE, BUT YOU MAY NOT
!  DISTRIBUTE EITHER THE ORIGINAL OR THE MODIFIED CODE TO USERS AT
!  ANY SITES OTHER THAN YOUR OWN.
!**********************************************************************

SUBROUTINE GraphicsXmgr(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,  &
    ndecay,ikin,nx,ny,nz,realtime,nn,nint,ikmast,ikph,delt)
USE crunchtype
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
!fp! auto_par_loops=0;

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

!  Internal variables and arrays

CHARACTER (LEN=13), DIMENSION(nrct)                :: uminprnt
CHARACTER (LEN=13), DIMENSION(ncomp+nspec)         :: ulabprnt
CHARACTER (LEN=mls)                                 :: fn
CHARACTER (LEN=mls)                                  :: suf
CHARACTER (LEN=mls)                                  :: suf1
CHARACTER (LEN=mls)                                 :: fnv
CHARACTER (LEN=1)                                  :: tab
CHARACTER (LEN=mls), DIMENSION(nsurf+nsurf_sec)    :: prtsurf
 
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

REAL(DP), DIMENSION(ncomp)                         :: totex_bas
REAL(DP), DIMENSION(nrct)                          :: dptprt
REAL(DP), DIMENSION(nrct)                          :: dsat
REAL(DP), DIMENSION(nrct)                          :: dvolpr
REAL(DP)                                           :: sum
REAL(DP)                                           :: porprt
REAL(DP)                                           :: phprt
REAL(DP)                                           :: porcalc

REAL(DP)                                                   :: sumiap
REAL(DP)                                                   :: pHprint
REAL(DP)                                                   :: peprint
REAL(DP)                                                   :: Ehprint
REAL(DP)                                                   :: spprint
REAL(DP)                                                   :: totcharge
REAL(DP)                                                   :: siprnt
REAL(DP)                                                   :: actprint
REAL(DP)                                                   :: actprint10
REAL(DP)                                                   :: spbase
REAL(DP)                                                   :: rone
REAL(DP)                                                   :: PrintTime
REAL(DP)                                                   :: alk
REAL(DP)                                                   :: tflux_top
REAL(DP)                                                   :: tflux_bot
REAL(DP)                                                   :: top_norm
REAL(DP)                                                   :: bot_norm
REAL(DP)                                                   :: aflux_net
REAL(DP)                                                   :: ad_net_bot

CHARACTER (LEN=mls)                                        :: namtemp

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: is

PrintTime = realtime*OutputTimeScale
rone = 1.0d0

suf='.out'
suf1 ='.out'
tab = CHAR(9)

DO k = 1,nrct
  uminprnt(k) = umin(k)
END DO
DO ik = 1,ncomp+nspec
  ulabprnt(ik) = ulab(ik)
END DO
DO ks = 1,nsurf
  prtsurf(ks) = namsurf(ks)
END DO
DO ns = 1,nsurf_sec
  prtsurf(ns+nsurf) = namsurf_sec(ns)
END DO

!  Write out master variable

IF (ikph /= 0) THEN
  fn='pH'
  ilength = 2
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,2282)
  jy = 1
  jz = 1
  DO jx = 1,nx
!fp! if_onproc({#expr# sp(ikph,jx,jy,jz) #});
    phprt =  -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
    WRITE(8,183) x(jx)*OutputDistanceScale,phprt
!fp! end_onproc();
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

fn='conc'
ilength = 4
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
IF (ikph /= 0) THEN
  WRITE(8,2288) (ulabprnt(ik),ik=1,ncomp+nspec)
  jy = 1
  jz = 1
  DO jx = 1,nx
!fp! if_onproc({#expr# sp(ik,jx,jy,jz) #});
    phprt =  -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
    WRITE(8,184) x(jx)*OutputDistanceScale, phprt,(sp(ik,jx,jy,jz)/clg,ik = 1,ncomp+nspec)
!        write(8,184) x(jx)*OutputDistanceScale,
!     1     (sp10(ik,jx,jy,jz),IK = 1,ncomp+nspec)
!        write(8,184) x(jx)*OutputDistanceScale,(sp10(ik,jx,jy,jz),IK = 1,ncomp)
!fp! end_onproc();
  END DO
ELSE
  WRITE(8,2285) (ulabprnt(ik),ik=1,ncomp+nspec)
  jy = 1
  jz = 1
  DO jx = 1,nx
!fp! if_onproc({#expr# sp(ik,jx,jy,jz) #});
    WRITE(8,184) x(jx)*OutputDistanceScale, (sp(ik,jx,jy,jz)/clg,ik = 1,ncomp+nspec)
!        write(8,184) x(jx)*OutputDistanceScale,
!     1     (sp10(ik,jx,jy,jz),IK = 1,ncomp+nspec)
!        write(8,184) x(jx)*OutputDistanceScale,(sp10(ik,jx,jy,jz),IK = 1,ncomp)
!fp! end_onproc();
  END DO
END IF
CLOSE(UNIT=8,STATUS='keep')

fn='totcon'
ilength = 6
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,2285) (ulabprnt(ik),ik=1,ncomp)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# s(i,jx,jy,jz) #});
  WRITE(8,184) x(jx)*OutputDistanceScale,(s(i,jx,jy,jz),i = 1,ncomp)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

fn='exchange'
ilength = 8
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,2285) (nam_exchsec(nex),nex=1,nexch_sec)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# spex10(nex,jx,jy,jz) #});
  WRITE(8,184) x(jx)*OutputDistanceScale,(spex10(nex+nexchange,jx,jy,jz),nex = 1,nexch_sec)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

fn='totexchange'
ilength = 11
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# spex10(nex,jx,jy,jz) #});
  totex_bas = 0.0
  DO i = 1,ncomp  
    DO nex = 1,nexch_sec
      totex_bas(i) = totex_bas(i) + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
    END DO
  END DO
  WRITE(8,184) x(jx)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

fn='surface'
ilength = 7
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,2285) (prtsurf(ks),ks=1,nsurf+nsurf_sec)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# spsurf10(ns,jx,jy,jz) #});
  WRITE(8,184) x(jx)*OutputDistanceScale,(spsurf10(ns,jx,jy,jz),ns = 1,nsurf+nsurf_sec)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

fn='totsurface'
ilength = 10
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# spsurf10(ns,jx,jy,jz) #});
  totex_bas = 0.0
  DO i = 1,ncomp  
    DO ns = 1,nsurf_sec
      totex_bas(i) = totex_bas(i) + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
    END DO
  END DO
  WRITE(8,184) x(jx)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

fn='totmineral'
ilength = 10
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# volfx(k,jx,jy,jz) #});
  totex_bas = 0.0
  DO i = 1,ncomp  
    DO k = 1,nrct
      IF (volmol(k) /= 0.0) THEN
        IF (nradmax > 0) THEN
          totex_bas(i) = totex_bas(i) + 0.001*mumin_decay(1,k,i,jx,1,1)*volfx(k,jx,jy,jz)/volmol(k)
        ELSE 
          totex_bas(i) = totex_bas(i) + 0.001*mumin(1,k,i)*volfx(k,jx,jy,jz)/volmol(k)
        END IF
      ENDIF
    END DO
  END DO
  WRITE(8,184) x(jx)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

!  Write out the reaction rates in units of mol/L(bulk vol.)/sec

fn='rate'

ilength = 4
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,2285)  (uminprnt(k),k=1,nrct)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# dppt(k,jx,jy,jz) #});
  sum = 0.0
  DO k = 1,nrct
!************************
!  For units of volume %/year, uncomment the following line and
!  recompile
!          dptprt(k) = dppt(k,jx,jy,jz)*volmol(k)*100.0  ! volume %/yr
!***********************
!************************
!  For units of mol/L(BV)/sec, uncomment the following line and
    dptprt(k) = dppt(k,jx,jy,jz)/(secyr*1000.0)    ! mol/L(BV)/sec
!*************************************
    sum = sum + dptprt(k)
  END DO
  porcalc = sum
  WRITE(8,184) x(jx)*OutputDistanceScale,(dptprt(k),k=1,nrct)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

!   Write out the reaction rates in units of mol/L(bulk vol.)/sec

fn='AqRate'

ilength = 6
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,2285)  (namkin(ir),ir=1,ikin)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# raq_tot(ir,jx,jy,jz) #});
  sum = 0.0
  WRITE(8,184) x(jx)*OutputDistanceScale,(raq_tot(ir,jx,jy,jz),ir=1,ikin)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

fn='volume'

!  Volumes in %

ilength = 6
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,2286)  (uminprnt(k),k=1,nrct)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# volfx(k,jx,jy,jz) #});
  sum = 0.0
  DO k = 1,nrct
    dvolpr(k) = 100*volfx(k,jx,jy,jz)
    sum = sum + volfx(k,jx,jy,jz)
  END DO
  porprt = (1.0-sum)*100.0
  WRITE(8,184) x(jx)*OutputDistanceScale,(dvolpr(k),k=1,nrct)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

fn = 'porosity'

ilength = 8
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# volfx(k,jx,jy,jz) #});
  sum = 0.0
  DO k = 1,nrct
    sum = sum + volfx(k,jx,jy,jz)
  END DO
  porprt = (1.0-sum)*100.0
  WRITE(8,184) x(jx)*OutputDistanceScale,porprt
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

!  Write out the saturation indices of the minerals (log Q/K).

fn='saturation'
ilength = 10
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,2285)  (uminprnt(k),k=1,nrct)
jy = 1
jz = 1
DO jx = 1,nx
!fp! if_onproc({#expr# sp(1,jx,jy,jz) #});
!  CALL reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt)
  CALL satcalc(ncomp,nrct,jx,jy,jz)
  DO k = 1,nrct
    dsat(k) = silog(1,k)
  END DO
  WRITE(8,184) x(jx)*OutputDistanceScale,(dsat(k),k=1,nrct)
!fp! end_onproc();
END DO
CLOSE(UNIT=8,STATUS='keep')

!  Write out the alkalinity

IF (giambalvo) THEN

fn='alkalinity'
ilength = 10
CALL NEWFILE(fn,suf1,fnv,nint,ilength)
OPEN(unit=8,file=fnv,access='sequential',status='unknown')
WRITE(8,2283) realtime
WRITE(8,2287)
DO jx = 1,nx
  alk = -s(ikph,jx,jy,jz)
  WRITE(8,184) x(jx)*OutputDistanceScale,alk
END DO
CLOSE(unit=8,status='keep')

!  Write out the advective and diffusive fluxes across top boundary (jx=1)
!  and bottom boundary (jx=nx).
!  The calculation (in fx.f) does not work when dispersion .ne. 0)
fn='flux'
ilength = 4
CALL newfile (fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv,ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) realtime
WRITE(8,2296) netflowx(0,1,1)
WRITE(8,2290)'adv top','dif top','tot top','adv bot',  &
    'dif bot','tot bot', 'Tb / Tt', 'nor top', 'nor bot'
WRITE(8,2291)
jy = 1
jz = 1
DO ik = 1, ncomp
  tflux_top = advflux_x(ik,jy,1)+dflux_x(ik,jy,1)
  tflux_bot = advflux_x(ik,jy,2)+dflux_x(ik,jy,2)
  top_norm = tflux_top / ABS(advflux_x(ik,jy,2))
  bot_norm = tflux_bot / ABS(advflux_x(ik,jy,2))
  WRITE(8,2292) ulabprnt(ik), advflux_x(ik,jy,1), dflux_x(ik,jy,1),  &
      tflux_top, advflux_x(ik,jy,2), dflux_x(ik,jy,2),  &
      tflux_bot, tflux_top/tflux_bot, top_norm, bot_norm
END DO
CLOSE(UNIT=8,STATUS='keep')

!  Write out the advective and diffusive fluxes across top boundary (jx=1)
!  with advective relative to seawater concentrations
!  The calculation (in fx.f) does not work when dispersion .ne. 0)
!  Because upflow is negative, >0 is SINK to ocean, <0 is SOURCE to ocean
!  Also write normalized flux (net out sediment/basement direct to ocean)

fn='netflux'
ilength = 7
CALL newfile (fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv,ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) realtime
WRITE(8,2297) netflowx(0,1,1)
WRITE(8,2298) netflowx(nx,1,1)
WRITE(8,2290)'NetAdvT','Dif top','NetTotT', 'NetAdvB','NTT/NAB','ConcSW'
WRITE(8,2293)
jy = 1
jz = 1
DO ik = 1, ncomp
  aflux_net = advflux_x(ik,jy,1) - sbnd(ik,1)*netflowx(0,jy,jz) *ro(1,jy,jz)
  tflux_top = aflux_net + dflux_x(ik,jy,1)
  
!  normalized to advective flux out of basement directly to ocean
  ad_net_bot = advflux_x(ik,jy,2) - sbnd(ik,1) * netflowx(nx,jy,jz)*ro(nx,jy,jz)
  top_norm = tflux_top / ad_net_bot
    WRITE(8,2294) ulabprnt(ik), aflux_net, dflux_x(ik,jy,1),  &
        tflux_top, ad_net_bot, top_norm, sbnd(ik,1)
END DO
CLOSE(UNIT=8,STATUS='keep')

END IF

502 FORMAT('temperature    ' ,f8.2)
503 FORMAT(a20,4X,1PE12.4)
504 FORMAT('END')
182 FORMAT(80(1X,1PE12.4))
183 FORMAT(1PE12.4,2X,1PE12.4)
184 FORMAT(100(1X,1PE14.6))

!2283 FORMAT('# Time (yrs) ',2X,1PE12.4)
2283 FORMAT('# Time      ',2X,1PE12.4)
2284 FORMAT('#   Distance ',a18)
2282 FORMAT('#   Distance ','        pH')
2281 FORMAT('#   Distance ',4X,a18)
2285 FORMAT('#   Distance    ',100(1X,a14))
2286 FORMAT('#   Distance    ',100(1X,a14))


600 FORMAT(2X,f10.2,2X,a15)
201 FORMAT(2X,a18,2X,f8.2)
202 FORMAT(2X,a18,2X,f8.3,3X,f8.3,2X,1PE12.3,2X,1PE12.3,2X,1PE12.3,2x,a8)
211 FORMAT(2X,a18,2X,f8.3,3X,f8.3,2X,1PE12.3,2X,1PE12.3,2X,'            ',2x,a8)
203 FORMAT(2X,a18)

2288 FORMAT('#   Distance    ','  pH           ',100(1X,a14))
2287 FORMAT('#   Distance ',' alkalinity (eq/kg)')
2290 FORMAT('#  Component',5X,12(a7,7X))
2291 FORMAT('#           ',5X,6('mol/m2/yr',5X))
2292 FORMAT(a15, 9(1X,1PE13.6), 2(1X,1I2), 2(1X,1PE13.6))
2293 FORMAT('#           ',5X,3('mol/m2/yr',5X))
2294 FORMAT(a15, 9(1X,1PE13.6))
2296 FORMAT('# Net flow at top: ',1X,1PE13.6)
2297 FORMAT('# Net flow at top:    ',1X,1PE13.6)
2298 FORMAT('# Net flow at bottom: ',1X,1PE13.6)



RETURN
END SUBROUTINE GraphicsXmgr
!  *******************************************************
