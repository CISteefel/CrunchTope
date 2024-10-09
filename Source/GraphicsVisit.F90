
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

SUBROUTINE GraphicsVisit(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,  &
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
LOGICAL(LGT), INTENT(IN)                           :: FirstCall

!  Internal variables and arrays

CHARACTER (LEN=13), DIMENSION(nrct)                :: uminprnt
CHARACTER (LEN=13), DIMENSION(ncomp+nspec)         :: ulabprnt
CHARACTER (LEN=mls)                                 :: fn
CHARACTER (LEN=mls)                                  :: suf
CHARACTER (LEN=mls)                                  :: suf1
CHARACTER (LEN=mls)                                 :: fnv
CHARACTER (LEN=1)                                  :: tab
CHARACTER (LEN=mls), DIMENSION(nsurf+nsurf_sec)    :: prtsurf
CHARACTER (LEN=mls)                                 :: char_time
CHARACTER (LEN=40)                                 :: prtspecies

INTEGER(I4B), DIMENSION(ncomp+nspec)               :: len_sp
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
INTEGER(I4B)                                       :: kk

REAL(DP), DIMENSION(ncomp)                         :: totex_bas
REAL(DP), DIMENSION(nrct)                          :: dptprt
REAL(DP), DIMENSION(nrct)                          :: dsat
REAL(DP), DIMENSION(nrct)                          :: dvolpr
REAL(DP)                                           :: sum
REAL(DP)                                           :: porprt
REAL(DP)                                           :: StressPrt
REAL(DP)                                           :: phprt
REAL(DP)                                           :: porcalc

REAL(DP)                                                   :: sumiap
REAL(DP)                                                   :: pHprint
REAL(DP)                                                   :: peprint
REAL(DP)                                                   :: fe2print
REAL(DP)                                                   :: fe3print
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
REAL(DP)                                                   :: AqueousToBulk
REAL(DP)                                                   :: SolidSolutionRatioTemp

REAL(DP), DIMENSION(ngas)                                  :: gastmp10
REAL(DP)                                                   :: denmol
REAL(DP)                                                   :: tk


REAL(DP), PARAMETER                                        :: zero = 0.0d0

CHARACTER (LEN=mls)                                        :: namtemp

INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: is

REAL(DP)                                                   :: CellVolume
REAL(DP)                                                   :: SumPorosity
REAL(DP)                                                   :: SumVolumeAllMinerals
REAL(DP), DIMENSION(nrct)                                  :: sumMineralRate
REAL(DP), DIMENSION(nrct)                                  :: sumVolumeMineral
REAL(DP), DIMENSION(nrct)                                  :: sumMoleMineral

CHARACTER (LEN=mls),DIMENSION(ncomp+nspec+nrct)                 :: WriteString
CHARACTER (LEN=mls)                                        :: StringProper
CHARACTER (LEN=mls)                                        :: StringTemp

REAL(DP), DIMENSION(ncomp)                                 :: sprint

REAL(DP)                                                   :: WritePermx
REAL(DP)                                                   :: WritePermy
REAL(DP)                                                   :: WritePermz

INTEGER(I4B)                                               :: id
INTEGER(I4B)                                               :: kIsotopologue
INTEGER(I4B)                                               :: isotopologue
INTEGER(I4B)                                               :: kMineralCommon
INTEGER(I4B)                                               :: kMineralRare
REAL(DP)                                                   :: totRare
REAL(DP)                                                   :: totCommon
REAL(DP), DIMENSION(ncomp)                                 :: IsotopeRatio

REAL(DP), DIMENSION(ncomp)                         :: gflux_ver
REAL(DP), DIMENSION(ncomp)                         :: gflux_hor
REAL(DP), DIMENSION(nrct)                          :: MineralPercent
REAL(DP), DIMENSION(ikin)                          :: dummy_raq_tot

!!!jz = 1
PrintTime = realtime*OutputTimeScale
rone = 1.0d0

suf='.tec'
suf1 ='.tec'
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

OPEN(UNIT=15,STATUS='scratch')
WRITE(15,*) realtime
REWIND 15
READ(15,'(a)') char_time
CLOSE(UNIT=15)

CALL stringlen(char_time,ls)

200 FORMAT(1PE9.2)


!****************
!! Begin default
!****************

fn='Aq_totconc'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Total Concentrations (mol/kgw)" '
  DO ik= 1, ncomp
    StringTemp = ulab(ik)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(ik) = StringProper(1:ls+3)
  END DO
    WRITE(8,2009) (WriteString(ik),ik=1,ncomp)
!!!  WRITE(8,2009) (ulab(ik),ik=1,ncomp)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          if (activecell(jx,jy,jz) == 0) THEN
            do i = 1,ncomp
              sprint(i) = -0.001 
            end do
          ELSE
            do i = 1,ncomp
              if (s(i,jx,jy,jz) < 1.0E-30) THEN
                sprint(i) = 1.0E-30
                ELSE
                sprint(i) = s(i,jx,jy,jz)
                END IF
            end do
          END IF
        
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,(sprint(i),i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='Aq_conc'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Total Concentrations (log mol/kgw)" '
  DO ik= 1, ncomp+nspec
    StringTemp = ulab(ik)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(ik) = StringProper(1:ls+3)
  END DO
    WRITE(8,2009) (WriteString(ik),ik=1,ncomp+nspec)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
                z(jz)*OutputDistanceScale,(sp(IK,jx,jy,jz)/clg,IK = 1,ncomp+nspec)
        END DO
      END DO
    END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='velocity'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Velocity (m/yr)" '
  WRITE(8,2012)
  2012 FORMAT('VARIABLES = "X"          "Y"              "Z"           "X Velocity"     "Y-Velocity"     "Z-Velocity" ')
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
                z(jz)*OutputDistanceScale,qx(jx,jy,jz),qy(jx,jy,jz),qz(jx,jy,jz)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn = 'temperature'
  ilength = 11
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Temperature (C)" '
  WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"     "Wcres" '
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,   &
      t(jx,jy,jz)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
  
  fn = 'porosity'
    ilength = 8
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,*) 'TITLE = "Porosity" '
    WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"          "Porosity" '
    WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          porprt = por(jx,jy,jz)*1.0
          WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,porprt
        END DO
      END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')
    
    fn = 'porositychange'
    ilength = 14
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,*) 'TITLE = "Porosity Change" '
    WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"          "PorosityChange" '
    WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          porprt = (porin(jx,jy,jz) - por(jx,jy,jz))
          WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,porprt
        END DO
      END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')

    fn='tortuosity'
  ilength = 12
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Tortuosity" '
  WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"     "Tortuosity" '
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,   &
                z(jz)*OutputDistanceScale,tortuosity(jx,jy,jz)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!****************
!! End default
!****************

!****************
!! Begin pH
!****************

  IF (ikph /= 0) THEN

    fn='pH'
    ilength = 2
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,*) 'TITLE = "Solution pH" '
    WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"             "pH" '
    WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          phprt =  -(sp(ikph,jx,jy,jz)+lngamma(ikph,jx,jy,jz))/clg
          WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,phprt
        END DO
      END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')

  END IF

!****************
!! End pH
!****************

!****************
!! Begin isotopes
!****************


IF (nIsotopePrimary > 0) THEN
  fn='toperatio_aq'
  ilength = 12
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Isotope ratios (per mil)" '
  DO id= 1, nIsotopePrimary
    StringTemp = nameIsotopeCommon(id)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(id) = StringProper(1:ls+3)
  END DO
  WRITE(8,2009) (WriteString(id),id=1,nIsotopePrimary)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        DO id = 1,nIsotopePrimary
          totCommon = s(isotopeCommon(id),jx,jy,jz)
          totRare = s(isotopeRare(id),jx,jy,jz)
          IsotopeRatio(id) = ( (totRare/totCommon)/IsotopeReference(id) - 1.0d0 ) *1000.0d0
        END DO
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,   &
                z(jz)*OutputDistanceScale,(IsotopeRatio(id),id = 1,nIsotopePrimary)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

IF (nIsotopeMineral > 0) THEN

  fn='toperatio_min'
  ilength = 13
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Isotope ratios (per mil)" '
  DO id= 1, nIsotopeMineral
    StringTemp = nameIsotopeMineralCommon(id)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(id) = StringProper(1:ls+3)
  END DO
  WRITE(8,2009) (WriteString(id),id=1,nIsotopeMineral)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
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
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,      &
                       (IsotopeRatio(kIsotopologue),kIsotopologue = 1,nIsotopeMineral)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

!****************
!! End isotopes
!****************


!****************
!! Begin minerals
!****************

IF (nrct > 0) THEN

  fn='Mineralrate'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Mineral Rate (mol/m^3/sec)" '
  DO k= 1, nrct
    StringTemp = umin(k)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(k) = StringProper(1:ls+3)
  END DO
    WRITE(8,2009) (WriteString(k),k=1,nrct)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        DO k = 1,nrct
          dptprt(k) = dppt(k,jx,jy,jz)/(secyr)    ! mol/m^3/s
          if (abs(dptprt(k)) < 1e-30) then
            dptprt(k) = 1e-30
          endif
        END DO
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,(dptprt(k),k=1,nrct)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn = 'Mineralvolumefraction'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Mineral Volumes (m^3 mineral/m^3 porous medium)" '
  DO k= 1, nrct
    StringTemp = umin(k)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(k) = StringProper(1:ls+3)
  END DO
    WRITE(8,2009) (WriteString(k),k=1,nrct)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        DO k = 1,nrct
          dvolpr(k) = volfx(k,jx,jy,jz)*1.0
        END DO
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,(dvolpr(k),k=1,nrct)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  !   fn='MineralVolfraction'
  !   ilength = 18
  !   CALL newfile(fn,suf1,fnv,nint,ilength)
  !   OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  !   WRITE(8,2283) PrintTime
  !   WRITE(8,131)
  !   131 FORMAT('# Units: % of vol occupied by mineral')
  !   WRITE(8,2285)  (uminprnt(k),k=1,nrct)
  !   jz = 1
  !   DO jy = 1,ny
  !   DO jx = 1,nx
  !     WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(volfx(k,jx,jy,jz),k = 1,nrct)
  !   END DO
  !  END DO
  !   CLOSE(UNIT=8,STATUS='keep')

  fn='Mineralsaturation'
  ilength = 5
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Mineral Saturation State (log Q/Keq)" '
  DO k= 1, nrct
    StringTemp = umin(k)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(k) = StringProper(1:ls+3)
  END DO
  WRITE(8,2009) (WriteString(k),k=1,nrct)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        CALL satcalc(ncomp,nrct,jx,jy,jz)
        DO k = 1,nrct
          dsat(k) = silog(1,k)
        END DO
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,(dsat(k),k=1,nrct)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='Mineralconc'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Total Component Concentration in Minerals (mol/m^3 PM)" '
  DO ik= 1, ncomp
    StringTemp = ulab(ik)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(ik) = StringProper(1:ls+3)
  END DO
    WRITE(8,2009) (WriteString(ik),ik=1,ncomp)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
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
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='Mineralarea'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Mineral Area (m^2/m^3 PM)" '
  DO k= 1, nrct
    StringTemp = umin(k)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(k) = StringProper(1:ls+3)
  END DO
  WRITE(8,2009) (WriteString(k),k=1,nrct)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,(area(k,jx,jy,jz),k=1,nrct)
      END DO
    END DO
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
    jz = 1
    DO jy = 1,ny
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
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale, &
      (MineralPercent(k),k = 1,nrct)
    END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')

  

END IF   

!****************
!! End of minerals
!****************


!****************
!! Begin aqueous reactions
!****************

IF (ikin > 0) THEN

  fn='AqRate'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Aqueous Rate (mol/L/yr) " '
  DO ir=1,ikin
    StringTemp = namkin(ir)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(ir) = StringProper(1:ls+3)
  END DO
    WRITE(8,2009) (WriteString(ir),ir=1,ikin)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        DO ir =1,ikin
          if (abs(raq_tot(ir,jx,jy,jz)) < 1d-70) THEN
            dummy_raq_tot(ir) = 1d-70
          else
            dummy_raq_tot(ir)=raq_tot(ir,jx,jy,jz)
          end if
        END DO
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,(dummy_raq_tot(ir),ir=1,ikin)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

!****************
!! End aqueous reactions
!****************


!****************
!!Begin Velocity field for velocity_read:
!****************

  IF (generate_velocity_vector) THEN
  fn='velocityx'
  ilength = 9
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Velocity (m/yr)" '
  WRITE(8,2012)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 0,nx
            WRITE(8,192) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
            z(jz)*OutputDistanceScale,qx(jx,jy,jz)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
  
  fn='divergence'
  ilength = 9
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Velocity (m/yr)" '
  WRITE(8,2012)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
            WRITE(8,192) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
            z(jz)*OutputDistanceScale,qx(jx,jy,jz) - qx(jx - 1,jy,jz)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='velocityy'
  ilength = 9
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Velocity (m/yr)" '
  WRITE(8,2012)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 0,ny
        DO jx = 1,nx
            WRITE(8,192) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
            z(jz)*OutputDistanceScale,qy(jx,jy,jz)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='velocityz'
  ilength = 9
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Velocity (m/yr)" '
  WRITE(8,2012)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 0,nz
      DO jy = 1,ny
        DO jx = 1,nx
            WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
            z(jz)*OutputDistanceScale,qz(jx,jy,jz)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

ENDIF

!****************
!!End Velocity field for velocity_read:
!****************


!****************
!! Begin Calculateflow
!****************

  IF (calculateflow) THEN

    fn='pressure'
    ilength = 8
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    115 FORMAT('# Units: Pascals')
    WRITE(8,*) 'TITLE = "Pressure" '
    WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"          "Pressure" '
    WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,pres(jx,jy,jz)
      END DO
    END DO
    ENDDO
    CLOSE(UNIT=8,STATUS='keep')

    fn = 'permeability'
    ilength = 12
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,*) 'TITLE = "Log Permeability (m^2)" '
    WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"     "X-Perm" "Y-Perm" "Z-Perm"'
    WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          IF (permx(jx,jy,jz) == 0.0) THEN
            WritePermx = -30.00
          ELSE
            WritePermx = Log10(permx(jx,jy,jz))
          END IF
          IF (permy(jx,jy,jz) == 0.0) THEN
            WritePermy = -30.00
          ELSE
            WritePermy = Log10(permy(jx,jy,jz))
          END IF

          IF (ny==1)  THEN
            WritePermy = 0.0
          ENDIF
          IF (nz==1)  THEN
            WritePermz = 0.0
          ENDIF
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,   &
             WritePermx, WritePermy, WritePermz
        END DO
      END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')

      IF (Richards) THEN
        fn='head'
          ilength = 8
          CALL newfile(fn,suf1,fnv,nint,ilength)
          OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
          WRITE(8,*) 'TITLE = "Pressure head (m)" '
          WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"          "Head" '
          WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
            DO jz = 1,nz
              DO jy = 1,ny
                DO jx = 1,nx
                  WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
                        z(jz)*OutputDistanceScale,head(jx,jy,jz)
              END DO
            END DO
          END DO
          CLOSE(UNIT=8,STATUS='keep')
          
          fn='water_potential'
          ilength = 8
          CALL newfile(fn,suf1,fnv,nint,ilength)
          OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
          WRITE(8,*) 'TITLE = "Water_potential (m)" '
          WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"          "water_potential" '
          WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
            DO jz = 1,nz
              DO jy = 1,ny
                DO jx = 1,nx
                  WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
                        z(jz)*OutputDistanceScale,psi(jx,jy,jz)
              END DO
            END DO
          END DO
          CLOSE(UNIT=8,STATUS='keep')

          fn = 'VG_alpha'
          ilength = 12
          CALL newfile(fn,suf1,fnv,nint,ilength)
          OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
          WRITE(8,*) 'TITLE = "Van Genuchten alpha (m-1)" '
          WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"     "VG_a"'
          WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx = 1,nx
              WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,   &
              VG_alpha(jx,jy,jz)
              END DO
            END DO
          END DO
          CLOSE(UNIT=8,STATUS='keep')
        
          fn = 'VG_n'
          ilength = 12
          CALL newfile(fn,suf1,fnv,nint,ilength)
          OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
          WRITE(8,*) 'TITLE = "Van Genuchten n (-)" '
          WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"     "VG_n"'
          WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx = 1,nx
              WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,   &
              VG_n(jx,jy,jz)
              END DO
            END DO
          END DO
          CLOSE(UNIT=8,STATUS='keep')
        
          fn = 'VG_wcres'
          ilength = 12
          CALL newfile(fn,suf1,fnv,nint,ilength)
          OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
          WRITE(8,*) 'TITLE = "Van Genuchten wat cont residual (-)" '
          WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"     "Wcres" '
          WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx = 1,nx
              WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,   &
              theta_r(jx,jy,jz)
              END DO
            END DO
          END DO
          CLOSE(UNIT=8,STATUS='keep')
        
          IF (evapofix .or. evapotimeseries) THEN
          fn = 'evapo_rate'
          ilength = 12
          CALL newfile(fn,suf1,fnv,nint,ilength)
          OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
          WRITE(8,*) 'TITLE = "evaporate [mm/yr]" '
          WRITE(8,*) 'VARIABLES = "evapo"'
          WRITE(8,*) ' F=POINT'
          WRITE(8,184) evaporate*1000+1e-9
          CLOSE(UNIT=8,STATUS='keep')
      END IF

      IF (transpifix .or. transpitimeseries) THEN
        fn = 'transpi_rate'
        ilength = 12
        CALL newfile(fn,suf1,fnv,nint,ilength)
        OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
        WRITE(8,*) 'TITLE = "transpirate [mm/yr]" '
        WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"     "VG_a"'
        WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
        DO jz = 1,nz
          DO jy = 1,ny
            DO jx = 1,nx
            WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,   &
            transpirate_cell(jx)*1000+1e-9
            END DO
          END DO
        END DO
        CLOSE(UNIT=8,STATUS='keep')
      ENDIF

    ENDIF

  ENDIF

  !****************
  !! End of Calculateflow
  !****************

!****************
!! Begin unsaturated
!****************

  IF (isaturate == 1) THEN

    fn='watercontent'
        ilength = 8
        CALL newfile(fn,suf1,fnv,nint,ilength)
        OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
        WRITE(8,*) 'TITLE = "Water content" '
        WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"          "Water Content" '
        WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx = 1,nx
                WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
                      z(jz)*OutputDistanceScale,satliq(jx,jy,jz)*por(jx,jy,jz)
            END DO
          END DO
        END DO
        CLOSE(UNIT=8,STATUS='keep')
        
        fn='liquid_saturation'
        ilength = 8
        CALL newfile(fn,suf1,fnv,nint,ilength)
        OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
        WRITE(8,*) 'TITLE = "Liquid saturation" '
        WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"          "Liquid saturation" '
        WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx = 1,nx
                WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
                      z(jz)*OutputDistanceScale,satliq(jx,jy,jz)
            END DO
          END DO
        END DO
        CLOSE(UNIT=8,STATUS='keep')
    
  IF (ngas > 0) THEN
      fn='gases_pp'
      ilength = 5
      CALL newfile(fn,suf1,fnv,nint,ilength)
      OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
      WRITE(8,*) 'TITLE = "Gas concentration (bars)" '
      DO kk=1,ngas
        StringTemp = namg(kk)
        CALL stringlen(StringTemp,ls)
        IF (ls > 14) THEN
          ls = 14
        END IF
        StringProper(1:1) = '"'
        StringProper(2:ls+1) = StringTemp(1:ls)
        StringProper(ls+2:ls+3) = '"'
        WriteString(kk) = StringProper(1:ls+3)
      END DO
        WRITE(8,2009) (WriteString(kk),kk=1,ngas)
      WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    
      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            tk = 273.15d0 + t(jx,jy,jz)
            denmol = 1.e05/(8.314*tk)                      ! P/RT = n/V, with pressure converted from bars to Pascals
            CALL GasPartialPressure(ncomp,ngas,gastmp10,jx,jy,jz)
            DO kk = 1,ngas
              if (gastmp10(kk)<1E-30) THEN
                gastmp10(kk)=1.0E-30
              END IF
            END DO
            WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,   &
                    z(jz)*OutputDistanceScale,(gastmp10(kk),kk = 1,ngas)
          END DO
        END DO
      END DO
      CLOSE(UNIT=8,STATUS='keep')
    
      fn='gases_conc'
        ilength = 10
        CALL newfile(fn,suf1,fnv,nint,ilength)
        OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
        WRITE(8,2283) PrintTime
        WRITE(8,130)
        130 FORMAT('# Units: mol/m3')
        WRITE(8,2285) (namg(kk),kk=1,ngas)
        jz = 1
        DO jy = 1,ny
        DO jx = 1,nx
          !!IF (activecellPressure(jx,jy,jz) == 0) THEN
          !!WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(spcondgas10(kk,jinit(jx,jy,1)),kk = 1,ngas)
          !!ELSE
          WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(spgas10(kk,jx,jy,1),kk = 1,ngas)
          !!END IF
        END DO
      END DO
        CLOSE(UNIT=8,STATUS='keep')

    IF (ny > 1) THEN
    fn='gasdifffluxY'
    ilength = 12
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,2283) PrintTime
    117 FORMAT('# Units: mol gas/m2/year')
    WRITE(8,117)
    WRITE(8,2285) (namg(kk),kk=1,ngas)
    jz = 1
    DO jy = 0,ny
    DO jx = 1,nx
      DO i = 1,ngas
        IF (jy == 0 .AND. activecellPressure(jx,jy+1,jz) == 0) THEN
        gflux_ver(i)=0
        ELSEIF (activecellPressure(jx,jy,jz) == 0 .AND. activecellPressure(jx,jy+1,jz) == 0) THEN
        gflux_ver(i)=0
        ELSEIF (jy == 0 .AND. activecellPressure(jx,jy+1,jz) == 1) THEN
        gflux_ver(i)=(fg(jx,jy+1,1))*(spgas10(i,jx,jy+1,1)-spcondgas10(i,jinit(jx,jy,1)))
        ELSEIF (jy<ny .AND. activecellPressure(jx,jy,jz) == 1 .AND. activecellPressure(jx,jy+1,jz) == 0) THEN
        gflux_ver(i)=(fg(jx,jy+1,1))*(spgas10(i,jx,jy+1,1)-spcondgas10(i,jinit(jx,jy,1)))
        ELSE
        gflux_ver(i)=(fg(jx,jy+1,1))*(spgas10(i,jx,jy+1,1)-spgas10(i,jx,jy,1))
        END IF
        if (abs(gflux_ver(i))<1.0E-30) THEN
          gflux_ver(i)=1.0E-30
        END IF
      END DO
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale+dyy(jy)/2,z(jz)*OutputDistanceScale, &
      (gflux_ver(i),i=1,ngas)
    END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')
    ENDIF

    IF (ny == 1 .AND. nz == 1) THEN
    
      fn='gasdifffluxX'
      ilength = 12
      CALL newfile(fn,suf1,fnv,nint,ilength)
      OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
      WRITE(8,2283) PrintTime
      WRITE(8,117)
      WRITE(8,2285) (namg(kk),kk=1,ngas)
      jz = 1
      jy = 1
      DO jx = 0,nx
        DO i = 1,ngas
          if (jx==0) THEN
            gflux_hor(i)=(ag(jx+1,1,1))*(spgas10(i,jx+1,1,1)-spcondgas10(i,jinit(jx,jy,jz)))
          elseif (jx==nx) THEN
            gflux_hor(i)= (cg(jx,1,1))*(spcondgas10(i,jinit(jx+1,jy,jz))-spgas10(i,jx,1,1))
          else
            gflux_hor(i)=(cg(jx,1,1))*(spgas10(i,jx+1,1,1)-spgas10(i,jx,1,1))
          END IF
          if (abs(gflux_hor(i))<1.0E-30) THEN
            gflux_hor(i)=1.0E-30
          END IF      
        END DO
        if (jx==0) THEN
        WRITE(8,184) x(jx)*OutputDistanceScale,(gflux_hor(i),i=1,ngas)
        else
        WRITE(8,184) x(jx)*OutputDistanceScale+dxx(jx)/2,(gflux_hor(i),i=1,ngas)
        END IF
      END DO
      CLOSE(UNIT=8,STATUS='keep')
    ENDIF

  ENDIF

  END IF

!****************
!! End unsaturated
!****************


!****************
!! Begin Makemovie
!****************


  IF (MakeMovie) THEN
    
    IF (FirstCall) THEN
      fn='VelocityEvolve'
      ilength = 14
      CALL newfile(fn,suf1,fnv,nint,ilength)
      OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
      WRITE(8,*) 'TITLE = "Velocity (m/yr)" '
      WRITE(8,2012)
      WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale    &
                    ,z(jz)*OutputDistanceScale,qx(jx,jy,jz),qy(jx,jy,jz),qz(jx,jy,jz)
          END DO
        END DO
      END DO
    ELSE
      WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
                z(jz)*OutputDistanceScale,qx(jx,jy,jz),qy(jx,jy,jz),qz(jx,jy,jz)
          END DO
        END DO
      END DO
    END IF
    
  END IF

!****************
!! End Makemovie
!****************
  
!****************
!! Begin Fracturenetwork
!****************
  
  if (nmmLogical .and. .not. FractureNetwork) THEN
    fn = 'stress'
    ilength = 6
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,*) 'TITLE = "Porosity" '
    WRITE(8,*) 'VARIABLES = "X"          "Y"              "Z"          "Stress (mPa)" '
    WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          StressPrt = stress(jx,jy,jz)*1.0E-06    !! in mPa
          WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,StressPrt
        END DO
      END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')
  END IF

!****************
!! End Fracturenetwork
!****************

!****************
!! Begin MontTerri
!****************

    IF (MontTerri) THEN
    
      fn = 'MontTerri-Slice1-'
      ilength = 17
      CALL newfile(fn,suf1,fnv,nint,ilength)
      OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
      WRITE(8,*) 'TITLE = "Tracer" '
      WRITE(8,*) 'VARIABLES = "X"        "Z"     "Tracer"'
      WRITE(8,*) 'ZONE I=', nx,  ',  K=',nz, ' F=POINT'
      
      jy = 30
      DO jz = 1,nz      
  !!!      DO jy = 1,ny
          DO jx = 1,nx
          WRITE(8,184) x(jx)*OutputDistanceScale, z(jz)*OutputDistanceScale, s(1,jx,jy,jz)
          END DO
          
  !!!      END DO
      END DO
      CLOSE(UNIT=8,STATUS='keep')
      
      fn = 'MontTerri-Porosity1-'
      ilength = 20
      CALL newfile(fn,suf1,fnv,nint,ilength)
      OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
      WRITE(8,*) 'TITLE = "Tracer" '
      WRITE(8,*) 'VARIABLES = "X"        "Z"     "Tracer"'
      WRITE(8,*) 'ZONE I=', nx,  ',  K=',nz, ' F=POINT'
      
      jy = 30
      DO jz = 1,nz      
  !!!      DO jy = 1,ny
          DO jx = 1,nx
          WRITE(8,184) x(jx)*OutputDistanceScale, z(jz)*OutputDistanceScale, por(jx,jy,jz)
          END DO
          
  !!!      END DO
      END DO
      CLOSE(UNIT=8,STATUS='keep')
      
    END IF

!****************
!! End MontTerri
!****************

!****************
!! begin Exchange
!****************

IF (nexchange > 0) THEN

  fn='exchange'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Exchanger Concentration (mol/g solid)" '
  DO nex = 1,nexch_sec
    StringTemp = nam_exchsec(nex)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(nex) = StringProper(1:ls+3)
  END DO
  WRITE(8,2009) (WriteString(nex),nex=1,nexch_sec)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale, &
                (spex10(nex+nexchange,jx,jy,jz)/SolidSolutionRatioTemp,nex = 1,nexch_sec)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='totexchange'
  ilength = 11
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Total Exchange Concentration (mol/g solid)" '
  DO ik= 1, ncomp
    StringTemp = ulab(ik)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(ik) = StringProper(1:ls+3)
  END DO
    WRITE(8,2009) (WriteString(ik),ik=1,ncomp)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
        totex_bas = 0.0
        DO i = 1,ncomp
          DO nex = 1,nexch_sec
            totex_bas(i) = totex_bas(i) + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
          END DO
        END DO
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,   &
                z(jz)*OutputDistanceScale,(totex_bas(i)/SolidSolutionRatioTemp,i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

!****************
!! End Exchange
!****************


!***********
!Begin surface complexation
!***********

IF (nsurf>0) THEN

  fn='surface'
  ilength = 7
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Surface Complex Concentration (mol/g solid)" '
  DO ns=1,nsurf_sec
    StringTemp = namsurf_sec(ns)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(ns) = StringProper(1:ls+3)
  END DO
  WRITE(8,2009) (WriteString(ns),ns=1,nsurf_sec)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale, &
                (spsurf10(ns,jx,jy,jz)/SolidSolutionRatioTemp,ns = 1,nsurf+nsurf_sec)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='totsurface'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "Total Surface Complex Concentration (mol/g solid)" '
  DO ik= 1, ncomp
    StringTemp = ulab(ik)
    CALL stringlen(StringTemp,ls)
    IF (ls > 14) THEN
      ls = 14
    END IF
    StringProper(1:1) = '"'
    StringProper(2:ls+1) = StringTemp(1:ls)
    StringProper(ls+2:ls+3) = '"'
    WriteString(ik) = StringProper(1:ls+3)
  END DO
    WRITE(8,2009) (WriteString(ik),ik=1,ncomp)
  WRITE(8,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
        totex_bas = 0.0
        DO i = 1,ncomp
          DO ns = 1,nsurf_sec
            totex_bas(i) = totex_bas(i) + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
          END DO
        END DO
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,   &
                z(jz)*OutputDistanceScale,(totex_bas(i)/SolidSolutionRatioTemp,i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

END IF

!***********
!End surface complexation
!***********



!!!2009 FORMAT('VARIABLES = "X"          "Y"              "Z"           ', 100(', "',A13,'"') )
!!!2009 FORMAT('VARIABLES = "X"          "Y"              "Z"           ', 100(' "'A13'" ') )
!!!2009 FORMAT('VARIABLES = "X"          "Y"              "Z"           ', 100('"', A13'"') )
!!!2009 FORMAT('VARIABLES = "X"          "Y"              "Z"           ', 100(' "',A13,'"') )
!!!2009 FORMAT('VARIABLES = "X"          "Y"              "Z"           ', 100('"',A13,'"') )
2009 FORMAT('VARIABLES = "X"          "Y"              "Z"           ',100(1X,A16) )

1022 FORMAT('VARIABLES = " X (meters)"," Y (meters)", " Z (meters)","X-Permeability", "Y-Permeability", "Z-Permeability"')
1010 FORMAT('VARIABLES = " X "," Y ", " Z ","Porosity"')
184 FORMAT(100(1X,1PE16.8))
1011 FORMAT('VARIABLES = "X"          "Y"              "Z"           ',100(', "',A13,'"') )






1009 FORMAT('VARIABLES = " X (meters) ", "  Y (meters)  "',100(', "',A13,'"'))

2001 FORMAT('VARIABLES = "X (meters) "',                   100(', "',A13,'"'))

1012 FORMAT('VARIABLES = " X (meters)", " Y (meters)", "X Velocity", "Y Velocity"')
1013 FORMAT('VARIABLES = " X (meters)", " Y (meters)", "Pressure"')
1015 FORMAT('VARIABLES = " X (meters)", "Pressure"')
2014 FORMAT('VARIABLES = " X (meters)", " Y (meters)", "Tortuosity"')
2015 FORMAT('VARIABLES = " X (meters)", "X_Velocity"')

1021 FORMAT('VARIABLES = " X (meters)", "X-Permeability"')

182 FORMAT(100(1X,1PE12.4))
183 FORMAT(1PE12.4,2X,1PE12.4,2X,1PE12.4)

191 FORMAT(100(1X,1PE16.7))
192 FORMAT(100(1X,1PE25.16))
188 FORMAT(100(1X,f15.7))

2283 FORMAT('# Time (yrs) ',2X,1PE12.4)
2284 FORMAT('      X        ','     Y        ',a18)
2282 FORMAT('   X        ','     Y        ','        pH')
2281 FORMAT('   X        ','     Y        ',4X,a18)
2285 FORMAT('    X        ','     Y        ',3X,30(1X,a13))
2286 FORMAT('    X        ','      Y               ',3X,30(1X,a15))

RETURN
END SUBROUTINE GraphicsVisit
!  *******************************************************
