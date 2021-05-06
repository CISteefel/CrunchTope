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

SUBROUTINE GraphicsKaleidagraph(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,  &
    ndecay,ikin,nx,ny,nz,realtime,nn,nint,ikmast,ikph,delt,jpor)
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
USE NanoCrystal
USE isotope

IMPLICIT NONE
!  *********************  INTERFACE BLOCKS  *****************************
INTERFACE
  SUBROUTINE GasPartialPressure(ncomp,ngas,gastmp10,jx,jy,jz)
    USE crunchtype
    INTEGER(I4B), INTENT(IN)                                   :: ncomp
    INTEGER(I4B), INTENT(IN)                                   :: ngas
    REAL(DP), DIMENSION(:)                                     :: gastmp10
    INTEGER(I4B), INTENT(IN)                                   :: jx
    INTEGER(I4B), INTENT(IN)                                   :: jy
    INTEGER(I4B), INTENT(IN)                                   :: jz
  END SUBROUTINE GasPartialPressure
END INTERFACE
!  **********************************************************************

!  External variables and arrays

REAL(DP), INTENT(IN)                               :: realtime
REAL(DP), INTENT(IN)                               :: delt

INTEGER(I4B), INTENT(IN)                           :: ncomp
INTEGER(I4B), INTENT(IN)                           :: nrct
INTEGER(I4B), INTENT(IN)                           :: nspec
INTEGER(I4B), INTENT(IN)                           :: ngas
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
INTEGER(I4B)                                       :: l

REAL(DP), DIMENSION(ncomp)                         :: totex_bas
REAL(DP), DIMENSION(ncomp)                         :: WeightPercent
REAL(DP), DIMENSION(nrct)                          :: MineralPercent
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
REAL(DP)                                                   :: SolidRatio
REAL(DP)                                                   :: FluidRatio
REAL(DP)                                                   :: Solid234U
REAL(DP)                                                   :: Solid238U
REAL(DP)                                                   :: SolidCa
REAL(DP)                                                   :: SolidRatioMarine1
REAL(DP)                                                   :: SolidRatioMarine2
REAL(DP)                                                   :: SolidSolutionRatioTemp

REAL(DP)                                                   :: AreaWrite3
REAL(DP)                                                   :: AreaWrite4

REAL(DP)                                                   :: Del34S_sulfate
REAL(DP)                                                   :: Del34S_sulfide
REAL(DP)                                                   :: DelCa44

REAL(DP)                                                   :: MineralVolumeSum
REAL(DP)                                                   :: PorosityFromSum

CHARACTER (LEN=mls)                                        :: namtemp

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: is
INTEGER(I4B)                                               :: kk

REAL(DP), DIMENSION(ngas)                                  :: gastmp10
REAL(DP)                                                   :: denmol
REAL(DP)                                                   :: tk

REAL(DP)                                                   :: MoleFraction40Mineral, MoleFraction44Mineral, MoleFraction40, MoleFraction44

REAL(DP)                                                        :: pi
REAL(DP)                                                        :: checkRatio

REAL(DP), DIMENSION(ikin)                                       :: satlogWrite
REAL(DP), DIMENSION(ikin)                                       :: satkinWrite

!! Isotopes
CHARACTER (LEN=mls),DIMENSION(ncomp+nspec+nrct)                 :: WriteString
CHARACTER (LEN=mls)                                        :: StringProper
CHARACTER (LEN=mls)                                        :: StringTemp
INTEGER(I4B)                                               :: id
INTEGER(I4B)                                               :: kIsotopologue
INTEGER(I4B)                                               :: isotopologue
INTEGER(I4B)                                               :: kMineralCommon
INTEGER(I4B)                                               :: kMineralRare
REAL(DP)                                                   :: totRare
REAL(DP)                                                   :: totCommon
REAL(DP), DIMENSION(ncomp)                                 :: IsotopeRatio
INTEGER(I4B)                                       :: ls

pi = DACOS(-1.0d0)

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
    phprt =  -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
    WRITE(8,183) x(jx)*OutputDistanceScale,phprt
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

fn='conc'
ilength = 4
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,101)
101 FORMAT('# Units: Log10 mol/kgw')

IF (ikph /= 0) THEN
  WRITE(8,2288) (ulabprnt(ik),ik=1,ncomp+nspec)
  jy = 1
  jz = 1
  DO jx = 1,nx
    phprt =  -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
    WRITE(8,184) x(jx)*OutputDistanceScale, phprt,(sp(ik,jx,jy,jz)/clg,ik = 1,ncomp+nspec)
  END DO
ELSE
  WRITE(8,2285) (ulabprnt(ik),ik=1,ncomp+nspec)
  jy = 1
  jz = 1
  DO jx = 1,nx
    WRITE(8,184) x(jx)*OutputDistanceScale, (sp(ik,jx,jy,jz)/clg,ik = 1,ncomp+nspec)
  END DO
END IF
CLOSE(UNIT=8,STATUS='keep')

IF (nIsotopePrimary > 0) THEN
  fn='toperatio_aq'
  ilength = 12
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Isotope ratios" '
  DO id= 1, nIsotopePrimary
    StringTemp = nameIsotopeCommon(id)
    CALL stringlen(StringTemp,ls)
    IF (ls > 16) THEN
      ls = 16
    END IF
    WriteString(id) = StringTemp(1:ls)
  END DO
  WRITE(8,2285) (WriteString(id),id=1,nIsotopePrimary)   
  jz = 1
  jy = 1
  DO jx = 1,nx
    DO id = 1,nIsotopePrimary
      totCommon = s(isotopeCommon(id),jx,jy,jz)
      totRare = s(isotopeRare(id),jx,jy,jz)
    checkRatio = totrare/totcommon
      IsotopeRatio(id) = ( (totRare/totCommon)/IsotopeReference(id) - 1.0d0 ) *1000.0d0
    END DO
    WRITE(8,184) x(jx)*OutputDistanceScale,(IsotopeRatio(id),id = 1,nIsotopePrimary)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

IF (nIsotopeMineral > 0) THEN

  fn='toperatio_min'
  ilength = 13
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Isotope ratios" '
  DO id= 1, nIsotopeMineral
    StringTemp = nameIsotopeMineralCommon(id)
    CALL stringlen(StringTemp,ls)
    IF (ls > 16) THEN
      ls = 16
    END IF
    WriteString(id) = StringTemp(1:ls)
  END DO
  WRITE(8,2285) (WriteString(id),id=1,nIsotopeMineral)   


  DO jx = 1,nx
    DO kIsotopologue = 1,nIsotopeMineral

      kMineralRare = kIsotopeRare(kIsotopologue)
      KMineralCommon = kIsotopeCommon(kIsotopologue)
      isotopologue = PointerToPrimaryIsotope(kIsotopologue)

      totCommon = volfx(kMineralCommon,jx,jy,jz)
      totRare   = volfx(kMineralRare,jx,jy,jz)
      IF (totCommon == 0.0) THEN
        IsotopeRatio(kIsotopologue) = 0.0
      ELSE
        IsotopeRatio(kIsotopologue) = ( (totRare/totCommon)/IsotopeReference(isotopologue) - 1.0d0 ) *1000.0d0
      END IF
    END DO
    WRITE(8,184) x(jx)*OutputDistanceScale,(IsotopeRatio(kIsotopologue),kIsotopologue = 1,nIsotopeMineral)
  END DO

  CLOSE(UNIT=8,STATUS='keep')

END IF

fn='totcon'
ilength = 6
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,102)
102 FORMAT('# Units: Mol/kgw')
WRITE(8,2285) (ulabprnt(ik),ik=1,ncomp)
jy = 1
jz = 1
DO jx = 1,nx
  WRITE(8,184) x(jx)*OutputDistanceScale,(s(i,jx,jy,jz),i = 1,ncomp)
END DO
CLOSE(UNIT=8,STATUS='keep')


fn='gases'
ilength = 5
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,103)
103 FORMAT('# Units: Bars')
WRITE(8,2285) (namg(kk),kk=1,ngas)
jy = 1
jz = 1
DO jx = 1,nx
  tk = 273.15d0 + t(jx,jy,jz)
  denmol = 1.e05/(8.314*tk)                      ! P/RT = n/V, with pressure converted from bars to Pascals
  CALL GasPartialPressure(ncomp,ngas,gastmp10,jx,jy,jz)
  WRITE(8,184) x(jx)*OutputDistanceScale,(gastmp10(kk),kk = 1,ngas)
END DO
CLOSE(UNIT=8,STATUS='keep')

IF (nexchange > 0) THEN
  fn='exchange'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,104)
!!  104 FORMAT('# Units: Mol/m^3 porous medium')
  104 FORMAT('# Units: Mol/g solid')
  WRITE(8,2285) (nam_exchsec(nex),nex=1,nexch_sec)
  jy = 1
  jz = 1
  DO jx = 1,nx
    SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
    WRITE(8,184) x(jx)*OutputDistanceScale,(spex10(nex+nexchange,jx,jy,jz)/SolidSolutionRatioTemp,nex = 1,nexch_sec)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='totexchange'
  ilength = 11
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,104)
  WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
  jy = 1
  jz = 1
  DO jx = 1,nx
    SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
    totex_bas = 0.0
    DO i = 1,ncomp  
      DO nex = 1,nexch_sec
        totex_bas(i) = totex_bas(i) + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
      END DO
    END DO
    WRITE(8,184) x(jx)*OutputDistanceScale,(totex_bas(i)/SolidSolutionRatioTemp,i = 1,ncomp)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF            !!  End of exchange block

IF (nsurf > 0) THEN
  fn='surface'
  ilength = 7
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,104)
  WRITE(8,2285) (prtsurf(ks),ks=1,nsurf+nsurf_sec)
  jy = 1
  jz = 1
  DO jx = 1,nx
    SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
    WRITE(8,184) x(jx)*OutputDistanceScale,(spsurf10(ns,jx,jy,jz)/SolidSolutionRatioTemp,ns = 1,nsurf+nsurf_sec)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='totsurface'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,104)
  WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
  jy = 1
  jz = 1
  DO jx = 1,nx
    totex_bas = 0.0
    DO i = 1,ncomp  
      DO ns = 1,nsurf_sec
        totex_bas(i) = totex_bas(i) + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
      END DO
    END DO
    WRITE(8,184) x(jx)*OutputDistanceScale,(totex_bas(i)/SolidSolutionRatioTemp,i = 1,ncomp)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

IF (nrct > 0) THEN
  fn='TotMineral'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,105)
  105 FORMAT('# Units: Mole component in mineral/m^3 porous medium')
  WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
  jy = 1
  jz = 1
  DO jx = 1,nx
    totex_bas = 0.0
    DO i = 1,ncomp  
      DO k = 1,nrct
        IF (volmol(k) /= 0.0) THEN
          IF (nradmax > 0) THEN
            totex_bas(i) = totex_bas(i) + mumin_decay(1,k,i,jx,1,1)*volfx(k,jx,jy,jz)/volmol(k)
          ELSE 
            totex_bas(i) = totex_bas(i) + mumin(1,k,i)*volfx(k,jx,jy,jz)/volmol(k)
          END IF
        ENDIF
      END DO
    END DO
    WRITE(8,184) x(jx)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='WeightPercent'
  ilength = 13
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,106)
106 FORMAT('# Units: Weight % Component')
  WRITE(8,2285) (ulabprnt(i),i=1,ncomp)
  jy = 1
  jz = 1
  DO jx = 1,nx
    totex_bas = 0.0
    sum = 0.0
    DO k = 1,nrct
      sum = sum + wtmin(k)*volfx(k,jx,jy,jz)/volmol(k)
    END DO
    DO i = 1,ncomp  
      DO k = 1,nrct
        IF (volmol(k) /= 0.0) THEN
          IF (nradmax > 0) THEN
            totex_bas(i) = totex_bas(i) + mumin_decay(1,k,i,jx,1,1)*volfx(k,jx,jy,jz)/volmol(k)
          ELSE 
            totex_bas(i) = totex_bas(i) + mumin(1,k,i)*volfx(k,jx,jy,jz)/volmol(k)
          END IF
        ENDIF
      END DO
    END DO
    DO i = 1,ncomp

      IF (sum > 0.0) THEN
        WeightPercent(i) = 100.0*totex_bas(i)*wtcomp(i)/sum
      ELSE
        WeightPercent(i) = 0.0
      END IF

    END DO
   
    WRITE(8,184) x(jx)*OutputDistanceScale,(WeightPercent(i),i = 1,ncomp)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='MineralPercent'
  ilength = 14
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,107)
107 FORMAT('# Units: Weight % Mineral')
  WRITE(8,2285)  (uminprnt(k),k=1,nrct)
  jy = 1
  jz = 1
  DO jx = 1,nx
    sum = 0.0
    DO k = 1,nrct
      sum = sum + wtmin(k)*volfx(k,jx,jy,jz)/volmol(k)
    END DO 
    DO k = 1,nrct
      IF (volmol(k) /= 0.0 .AND. sum /= 0.0) THEN
        MineralPercent(k) = 100.0*wtmin(k)*volfx(k,jx,jy,jz)/volmol(k)/sum
        IF (MineralPercent(k) < 1.0E-30) THEN
          MineralPercent(k) = 1.0E-30
        END IF
      END IF
    END DO
    WRITE(8,184) x(jx)*OutputDistanceScale,(MineralPercent(k),k = 1,nrct)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!  Write out the reaction rates in units of mol/m**3(bulk vol.)/sec

  fn='rate'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,108)
108 FORMAT('# Units: Mol/m**3 Porous Medium/s')
  WRITE(8,2285)  (uminprnt(k),k=1,nrct)
  jy = 1
  jz = 1
  DO jx = 1,nx
    sum = 0.0
    DO k = 1,nrct
!************************
!  For units of volume %/year, uncomment the following line and
!  recompile
!          dptprt(k) = dppt(k,jx,jy,jz)*volmol(k)*100.0  ! volume %/yr
!***********************
!************************
!  For units of mol/L(BV)/sec, uncomment the following line and
!!        dptprt(k) = dppt(k,jx,jy,jz)/(secyr*1000.0)    ! mol/L(BV)/sec
!!        dptprt(k) = dppt(k,jx,jy,jz)/(por(jx,jy,jz))    ! mol/L fluid/yr
        dptprt(k) = dppt(k,jx,jy,jz)/(secyr)    ! mol/m**3 porous medium/s
!*************************************
      sum = sum + dptprt(k)
    END DO
    porcalc = sum
    WRITE(8,184) x(jx)*OutputDistanceScale,(dptprt(k),k=1,nrct)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

ENDIF

!   Write out the reaction rates in units of mol/kgw/sec

IF (ikin > 0) THEN

  fn='AqRate'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,109)
109 FORMAT('# Units: Mol/kgw/yr')
  WRITE(8,2285)  (namkin(ir),ir=1,ikin)
  jy = 1
  jz = 1
  DO jx = 1,nx
    sum = 0.0
    WRITE(8,184) x(jx)*OutputDistanceScale,(raq_tot(ir,jx,jy,jz),ir=1,ikin)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

IF (ikin > 0) THEN
  fn='AqSat'
  ilength = 5
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,119)
119 FORMAT('# Units: Log Q/Keq (dimensionless)')
  WRITE(8,2285)  (namkin(ir),ir=1,ikin)
  jy = 1
  jz = 1
  DO jx = 1,nx
    DO ir = 1,ikin
      sum = 0.0d0
      DO i = 1,ncomp
        sum = sum + mukin(ir,i)*sp(i,jx,jy,jz)
      END DO
      satlogWrite(ir) = sum - clg*keqkin(ir)
      satkinWrite(ir) = satlog(ir,jx,jy,jz)/clg
    END DO
    WRITE(8,184) x(jx)*OutputDistanceScale,(satkinWrite(ir),ir=1,ikin)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

IF (ikin > 0) THEN
  fn='DelGBiomass'
  ilength = 11
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,109)
  WRITE(8,2285)  (namkin(ir),ir=1,ikin)
  jy = 1
  jz = 1
  DO jx = 1,nx
    Tk = 273.15d0 + t(jx,jy,jz)
    sum = 0.0
    WRITE(8,184) x(jx)*OutputDistanceScale,(rgas*Tk*satlog(ir,jx,jy,jz),ir=1,ikin)
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

IF (ikin > 0) THEN
  fn='fTBiomass'
  ilength = 9
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,109)
  WRITE(8,2285)  (namkin(ir),ir=1,ikin)
  jy = 1
  jz = 1
  DO jx = 1,nx
    Tk = 273.15d0 + t(jx,jy,jz)
    sum = 0.0
    WRITE(8,184) x(jx)*OutputDistanceScale,(1.0-DEXP(satlog(ir,jx,jy,jz)),ir=1,ikin)
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nrct > 0) THEN

!!  Volumes in %
  fn='volume'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,110)
  110 FORMAT('# Units: Volume % (m^3 mineral/m^3 porous medium)')
  WRITE(8,2286)  (uminprnt(k),k=1,nrct)
  jy = 1
  jz = 1
  DO jx = 1,nx
    sum = 0.0
    DO k = 1,nrct
      if (mintype(k) == 0) then
        dvolpr(k) = 100*volfx(k,jx,jy,jz)
        sum = sum + volfx(k,jx,jy,jz)
        if (dvolpr(k) < 1.0E-30 ) THEN
          dvolpr(k) = 0.0
        end if
      else if (mintype(k) == 1) then
        dvolpr(k) = volfx(k,jx,jy,jz)
      end if
    END DO
    porprt = (1.0-sum)*100.0
    WRITE(8,184) x(jx)*OutputDistanceScale,(dvolpr(k),k=1,nrct)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!!  Bulk areas in m2/m3 
  fn='area'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,111)
  111 FORMAT('# Units: M^2 mineral/m^3 porous medium')
  WRITE(8,2286)  (uminprnt(k),k=1,nrct)
  jy = 1
  jz = 1
  DO jx = 1,nx
    if (area(k,jx,jy,jz) < 1.0E-35) then
      area(k,jx,jy,jz) = 0.0d0
    end if
    WRITE(8,184) x(jx)*OutputDistanceScale,(area(k,jx,jy,jz),k=1,nrct)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

fn = 'porosity'
ilength = 8
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,112)
112 FORMAT('    Distance       Transport_Porosity  Mineral_Porosity')
jy = 1
jz = 1
DO jx = 1,nx
  MineralVolumeSum = 0.0d0
  DO k = 1,nrct
    MineralVolumeSum = MineralVolumeSum + volfx(k,jx,jy,jz)
  END DO
  porprt = por(jx,jy,jz)
  PorosityFromSum = 1.0d0 - MineralVolumeSum
  WRITE(8,184) x(jx)*OutputDistanceScale,porprt, PorosityFromSum
END DO
CLOSE(UNIT=8,STATUS='keep')

IF (CalculateFlow) THEN
fn = 'permeability'
ilength = 12
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,123)
123 FORMAT('# Units:           X Permeability')
jy = 1
jz = 1
DO jx = 1,nx
  WRITE(8,184) x(jx)*OutputDistanceScale,permx(jx,jy,jz)
END DO
CLOSE(UNIT=8,STATUS='keep')

fn = 'DarcyFlux'
ilength = 9
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,124)
124 FORMAT('# Units:           X Darcy Flux')
jy = 1
jz = 1
DO jx = 0,nx
  WRITE(8,184) x(jx)*OutputDistanceScale,qx(jx,jy,jz)
END DO
CLOSE(UNIT=8,STATUS='keep')

END IF

!  Write out the saturation indices of the minerals (log Q/K).

IF (nrct > 0) THEN
  fn='saturation'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,113)
113 FORMAT('# Units: Dimensionless (Log10(Q/Keq)')
  WRITE(8,2285)  (uminprnt(k),k=1,nrct)
  jy = 1
  jz = 1
  DO jx = 1,nx
!!!    CALL satcalc(ncomp,nrct,jx,jy,jz)
    DO k = 1,nrct
!!!      dsat(k) = silog(1,k)
      dsat(k) = silogGlobal(1,k,jx,jy,jz)
    END DO
    WRITE(8,184) x(jx)*OutputDistanceScale,(dsat(k),k=1,nrct)
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

!  Write out pressure

IF (calculateflow) THEN
  fn='pressure'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,2283) PrintTime
  WRITE(8,115)
  115 FORMAT('# Units: Pascals')
    DO jx = 1,nx
      WRITE(8,184) x(jx)*OutputDistanceScale,pres(jx,jy,jz)
    END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

!  Write out Darcy fluxes

write(*,*) ' Writing velocity file '
fn='velocity'
ilength = 8
CALL newfile(fn,suf1,fnv,nint,ilength)
OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
WRITE(8,2283) PrintTime
WRITE(8,116)
116 FORMAT('# Units: m^3 fluid/m^2 porous medium/yr')
jy = 1
jz = 1
DO jx = 0,nx
  WRITE(8,184) x(jx)*OutputDistanceScale,qx(jx,1,1)
END DO
CLOSE(UNIT=8,STATUS='keep')

502 FORMAT('temperature    ' ,f8.2)
503 FORMAT(a20,4X,1PE12.4)
504 FORMAT('END')
182 FORMAT(80(1X,1PE12.4))
183 FORMAT(1PE12.4,2X,1PE12.4)
184 FORMAT(100(1X,1PE16.4))

!2283 FORMAT('# Time (yrs) ',2X,1PE12.4)
2283 FORMAT('# Time      ',2X,1PE12.4)
2284 FORMAT('    Distance ',a18)
2282 FORMAT('    Distance ','        pH')
2281 FORMAT('    Distance ',4X,a18)
2285 FORMAT('    Distance       ',100(1X,a16))
2299 FORMAT('    Distance    ',100(1X,a16))
2288 FORMAT('    Distance    ','  pH           ',100(1X,a16))
2286 FORMAT('    Distance    ',100(1X,a16))
2289 FORMAT('    Distance    ',1x,' CalciteAR',1x,' PoreWaterAR', 1x, '234U', 1x, '238U', 1x, 'Ca')
514 FORMAT(1X,a10,1X,1PE11.4)
513 FORMAT(1X,a18,1X,1PE12.4,5X,a12)

600 FORMAT(2X,f10.2,2X,a15)
201 FORMAT(2X,a18,2X,f8.2)
202 FORMAT(2X,a18,2X,f8.3,3X,f8.3,2X,1PE12.3,2X,1PE12.3,2X,1PE12.3,2x,a8)
211 FORMAT(2X,a18,2X,f8.3,3X,f8.3,2X,1PE12.3,2X,1PE12.3,2X,'            ',2x,a8)
203 FORMAT(2X,a18)
207 FORMAT('                 ','          Log',1X,'       Log',1x,  &
'              ',1x,'          ', '       Activity')
204 FORMAT(' Species         ','     Molality',1X,'  Activity',1x,  &
'     Molality ',1x,'    Activity ', '  Coefficient','    Type')
205  FORMAT(2X,'Total Charge    = ',1pe10.3)

509 FORMAT(2X,a18,2X,f12.4)
510 FORMAT(2X,'GEOCHEMICAL CONDITION NUMBER',i3)
555 FORMAT(2X,'Temperature (C) = ',f10.3)
411 FORMAT(2X,'Ionic Strength  = ',f10.3)
5022 FORMAT(2X,'Solution pH     = ',f10.3)
5023 FORMAT(2X,'Solution pe     = ',f10.3)
5024 FORMAT(2X,'Solution Eh     = ',f10.3)
512 FORMAT(' Basis species    ','     Molality  ', '  Constraint type')
511 FORMAT(1X,a18,1X,1PE12.4,5X,a13,1X,a18)

542  FORMAT(5X,'Primary species ',a18,1X,'at grid pt. ',i5)
544  FORMAT(5X,'Check geochem. condition # ',i2)
543  FORMAT(5X,'Constraint mineral ',2X,a18)
643  FORMAT(5X,'Constraint gas',2X,a18)
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
END SUBROUTINE GraphicsKaleidagraph
!  *******************************************************
