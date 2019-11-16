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

SUBROUTINE GraphicsTecplot(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,  &
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

jz = 1
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

IF (ny > 1 .AND. nz > 1) THEN                 !! 3D case
    
  IF (CalculateFlow) THEN
    fn = 'permeability'
    ilength = 12
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,118)
    WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
    WRITE(8,1022)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny, ', K=',nz  
    DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,   &
             Log10(permx(jx,jy,jz)),Log10(permy(jx,jy,jz)),Log10(permz(jx,jy,jz))
      END DO
    END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')
  END IF

1022 FORMAT('VARIABLES = " X (meters)"," Y (meters)", " z (meters)","X-Permeability", "Y-Permeability", "Z-Permeability"')
     
  fn = 'porosity'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,112)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1010)
  WRITE(8,*) 'ZONE F=POINT,I=', nx 
  do jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      porprt = por(jx,jy,jz)*1.0
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,z(jz)*OutputDistanceScale,porprt
    END DO
  END DO
  end do
  CLOSE(UNIT=8,STATUS='keep')
  
  fn='totcon'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,102)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2009) (ulab(ik),ik=1,ncomp)
  do jy = 1,ny
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
    DO jz = 1,nz
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(s(i,jx,jy,jz),i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
  
  fn='totslice'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,102)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2009) (ulab(ik),ik=1,ncomp)
  jy = ny/2
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
    DO jz = 1,nz
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(s(i,jx,jy,jz),i = 1,ncomp)
      END DO
    END DO

  CLOSE(UNIT=8,STATUS='keep')
  
  fn='totline'
  ilength = 7
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,102)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2009) (ulab(ik),ik=1,ncomp)
  jy = ny/2
  jz = nz/2
  WRITE(8,*) 'ZONE F=POINT,I=', nx
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,(s(i,jx,jy,jz),i = 1,ncomp)
      END DO
  CLOSE(UNIT=8,STATUS='keep')
  
 
  116 FORMAT('# Units: m^3 fluid/m^2 porous medium/yr')

  fn='velocity'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,116)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2012)
  do jy = 1,ny
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
    DO jz = 1,nz
      DO jx = 1,nx
        WRITE(8,191) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,qx(jx,jy,jz),qz(jx,jy,jz)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!  Write out tortuosity

  117 FORMAT('# Units: Dimensionless')

  fn='tortuosity'
  ilength = 12
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,117)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2014)
  do jy = 1,ny
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
    DO jz = 1,nz
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,tortuosity(jx,jy,jz)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')


ELSE IF (ny > 1 .AND. nz == 1) THEN             !!  2D case, X and Y

  IF (FirstCall) THEN
    fn='VelocityEvolve'
    ilength = 14
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,1012)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    jz = 1
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,qx(jx,jy,jz),qy(jx,jy,jz)
      END DO
    END DO
  ELSE
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    jz = 1
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,qx(jx,jy,jz),qy(jx,jy,jz)
      END DO
    END DO
  END IF

  IF (ikph /= 0) THEN
    fn='pH'
    ilength = 2
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
    WRITE(8,*) 'VARIABLES = " X (meters)"," Y (meters)"," pH"'
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    jz = 1
    DO jy = 1,ny
      DO jx = 1,nx
        phprt =  -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
        WRITE(8,183) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,phprt
      END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')
  END IF

!  Write out all of the species concentrations

  DO ik = 1,ncomp+nspec
    CALL stringlen(ulab(ik),len_sp(ik))
  END DO

!!  Individual  species concentrations (log units)


101 FORMAT('# Units: Log10 mol/kgw')

  fn='conc'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,101)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1009) (ulab(ik),ik=1,ncomp+nspec)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(sp(IK,jx,jy,jz)/clg,IK = 1,ncomp+nspec)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!! Total concentrations

  102 FORMAT('# Units: Mol/kgw')

  fn='totcon'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,102)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1009) (ulab(ik),ik=1,ncomp)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(s(i,jx,jy,jz),i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

IF (nexchange>0) THEN
  fn='exchange'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1009) (nam_exchsec(nex),nex=1,nexch_sec)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(spex10(nex+nexchange,jx,jy,jz),nex = 1,nexch_sec)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF


103 FORMAT('# Units: mol component/m^3 gas')

IF (isaturate == 1) THEN
!!  fn='totgas'
  fn='gas'
  ilength = 3
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,123)
123 FORMAT('# Units: Bars')
  WRITE(8,2286) (namg(kk),kk=1,ngas)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
!!        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(sgas(i,jx,jy,jz),i = 1,ncomp)
        tk = 273.15d0 + t(jx,jy,jz)
        denmol = 1.e05/(8.314*tk)                      ! P/RT = n/V, with pressure converted from bars to Pascals
        CALL GasPartialPressure(ncomp,ngas,gastmp10,jx,jy,jz)
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(gastmp10(kk),kk = 1,ngas)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF


!!  104 FORMAT('# Units: Mol/m^3 porous medium')
  104 FORMAT('# Units: Mol/g solid')

IF (nexchange>0) THEN
  fn='totexchange'
  ilength = 11
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1009) (ulab(i),i=1,ncomp)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        totex_bas = 0.0
        SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
        DO i = 1,ncomp  
          DO nex = 1,nexch_sec
            totex_bas(i) = totex_bas(i) + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
          END DO
        END DO
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(totex_bas(i)/SolidSolutionRatioTemp,i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nsurf>0) THEN
  fn='surface'
  ilength = 7
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1009) ( namsurf(is),is=1,nsurf ),(namsurf_sec(ns),ns=1,nsurf_sec)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        SolidSolutionRatioTemp = 1000.d0*SolidDensity(jinit(jx,jy,jz))*(1.0-por(jx,jy,jz))
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(spsurf10(ns,jx,jy,jz)/SolidSolutionRatioTemp,ns = 1,nsurf+nsurf_sec)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nsurf>0) THEN
  fn='totsurface'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1009) (ulab(i),i=1,ncomp)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
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
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(totex_bas(i)/SolidSolutionRatioTemp,i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

105 FORMAT('# Units: Mole component in mineral/m^3 porous medium')

IF (nrct > 0) THEN
  fn='TotMineral'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,105)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1009) (ulab(i),i=1,ncomp)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
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
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!  Write out the reaction rates in units of mol/L(bulk vol.)/sec

!  Temporarily, rates are in units of mol/L fluid/yr

  DO k = 1,nrct
    CALL stringlen(umin(k),len_min(k))
  END DO

108 FORMAT('# Units: Mol/L Porous Medium/s')

IF (nrct > 0) THEN
  fn='rate'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,108)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1011) (umin(k),k=1,nrct)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  sumMineralRate = 0.0d0
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      sum = 0.0
      DO k = 1,nrct
!************************
!  For units of volume %/year, uncomment the following line and
!  recompile
!              dptprt(k) = dppt(k,jx,jy,jz)*volmol(k)*100.0  ! volume %/yr
!***********************
!************************
!  For units of mol/L(BV)/sec, uncomment the following line and
!!        dptprt(k) = dppt(k,jx,jy,jz)/(secyr*1000.0)    ! mol/L(BV)/sec
        dptprt(k) = dppt(k,jx,jy,jz)/(secyr*1000.0*por(jx,jy,jz))    ! mol/L fluid/sec
        CellVolume = dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
        sumMineralRate(k) = sumMineralRate(k) + CellVolume*dppt(k,jx,jy,jz)/(secyr)    !!  mol/m^3 PorMed/sec
!*************************************
        sum = sum + dptprt(k)
      END DO
      porcalc = sum
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(dptprt(k),k=1,nrct)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
  
  fn='BulkRate'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,108)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1011) (umin(k),k=1,nrct)
  WRITE(8,184) (sumMineralRate(k),k=1,nrct) 
  CLOSE(UNIT=8,STATUS='keep')
  
END IF

!   Write out the reaction rates in units of mol/L(bulk vol.)/sec

109 FORMAT('# Units: Mol/kgw/s')

IF (ikin > 0) THEN
  fn='AqRate'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,109)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1011)  (namkin(ir),ir=1,ikin)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
    WRITE(8,185) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(raq_tot(ir,jx,jy,jz),ir=1,ikin)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!  Volumes in %

  110 FORMAT('# Units: Volume % (m^3 mineral/m^3 porous medium)')

IF (nrct > 0) THEN
  fn = 'volume'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,110)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1011) (umin(k),k=1,nrct)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  SumVolumeMineral = 0.0d0
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      CellVolume = dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
      DO k = 1,nrct
        dvolpr(k) = volfx(k,jx,jy,jz)*1.0
        SumVolumeMineral(k) = SumVolumeMineral(k) + volfx(k,jx,jy,jz)*CellVolume
      END DO
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(dvolpr(k),k=1,nrct)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

  fn='BulkMoles'
  ilength = 9
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,108)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  DO k = 1,nrct
    SumMoleMineral(k) = SumVolumeMineral(k)/volmol(k)
  END DO
  WRITE(8,1011) (umin(k),k=1,nrct)
  WRITE(8,184) (SumMoleMineral(k),k=1,nrct) 
  CLOSE(UNIT=8,STATUS='keep')
  
  fn='BulkPorosity'
  ilength = 12
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,108)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  sumVolumeAllMinerals = 0.0d0
  SumPorosity = 0.0d0
  DO k = 1,nrct
    sumVolumeAllMinerals = sumVolumeAllMinerals + SumVolumeMineral(k)
  END DO
  SumPorosity = 1.0d0 - sumVolumeAllMinerals/(x(nx)*y(ny))
  WRITE(8,*) ' Bulk Porosity (m^3/m^3) = ',SumPorosity 
  CLOSE(UNIT=8,STATUS='keep')

!!  Bulk areas in m2/m3 
  fn='area'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,111)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1011) (umin(k),k=1,nrct)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      WRITE(8,184) x(jx)*OutputDistanceScale, y(jy)*OutputdistanceScale,(area(k,jx,jy,jz),k=1,nrct)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  111 FORMAT('# Units: M^2 mineral/m^3 porous medium')

!  Write out the porosity

112 FORMAT('# Units: % Porosity')

IF (ny == 1) THEN
    
  fn = 'porosity'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,112)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1010)
  WRITE(8,*) 'ZONE F=POINT,I=', nx 
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      porprt = por(jx,jy,jz)*1.0
      WRITE(8,184) x(jx)*OutputDistanceScale,porprt
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
  
ELSE
    
  fn = 'porosity'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,112)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1010)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      porprt = por(jx,jy,jz)*100.0
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,porprt
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
  
END IF

  1010 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "Porosity"')

!  Write out the permeability

  118 FORMAT('# Units: Log10 m^2 ')

  IF (CalculateFlow) THEN
    fn = 'permeability'
    ilength = 12
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,118)
    WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
    WRITE(8,1020)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    jz = 1
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,Log10(permx(jx,jy,jz)),Log10(permy(jx,jy,jz))
      END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')
  END IF

  1020 FORMAT('VARIABLES = " X (meters)"," Y (meters)", "X-Permeability", "Y-Permeability"')


!  Write out the saturation indices of the minerals (log Q/K).

113 FORMAT('# Units: Dimensionless (Log10(Q/Keq)')
 
IF (nrct > 0) THEN
  fn='saturation'
  ilength = 5
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,113)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1011) (umin(k),k=1,nrct)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      CALL satcalc(ncomp,nrct,jx,jy,jz)
      DO k = 1,nrct
        dsat(k) = silog(1,k)
      END DO
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(dsat(k),k=1,nrct)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!  Write out Darcy fluxes

  fn='velocity'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,116)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1012)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  jz = 1

!!  WRITE(8,188) x(0)*OutputDistanceScale,y(0)*OutputDistanceScale,zero,zero
!!  DO jx = 1,nx
!!    WRITE(8,188) x(jx)*OutputDistanceScale,y(0)*OutputDistanceScale,zero,qy(jx,0,jz)
!!  END DO
  DO jy = 1,ny
!!    WRITE(8,191) x(0)*OutputDistanceScale,y(jy)*OutputDistanceScale,qx(0,jy,jz),zero
    DO jx = 1,nx
      WRITE(8,191) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,qx(jx,jy,jz),qy(jx,jy,jz)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!  Write out tortuosity

  jz = 1
  fn='tortuosity'
  ilength = 12
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,117)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2014)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,tortuosity(jx,jy,jz)
      END DO
    END DO
  CLOSE(UNIT=8,STATUS='keep')

!  Write out temperatures

  114 FORMAT('# Units: Degrees C')

  fn='temperature'
  ilength = 11
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,114)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1012)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jy = 1,ny
    DO jx = 1,nx
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,t(jx,jy,jz)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!  Write out pressure

  115 FORMAT('# Units: Pascals')

  IF (calculateflow) THEN
  fn='pressure'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,115)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1013)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jy = 1,ny
    DO jx = 1,nx
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,pres(jx,jy,jz)
!      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,pres(jx,jy,jz) &
!         /ro(jx,jy,jz)/9.81 - (COSD(x_angle)*x(jx)+COSD(y_angle)*x(jy)+ & 
!         COSD(x_angle)*x(jx)+COSD(z_angle)*x(jz))*SignGravity

    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')


!  Write out heads

  fn='heads'
  ilength = 5
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1013)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jy = 1,ny
    DO jx = 1,nx
!fp! if_onproc({#expr# qx(jx,jy,jz) #});
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,pres(jx,jy,jz)
!!     WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,pres(jx,jy,jz) &
!!         /ro(jx,jy,jz)/9.81 - (COSD(x_angle)*x(jx)+ COSD(y_angle)*y(jy)+ & 
!!         COSD(x_angle)*x(jx) + COSD(z_angle)*z(jz))*SignGravity

!fp! end_onproc();
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  END IF

  119 FORMAT('# Units: m^3 gas/m^2 porous medium/yr')



IF (isaturate == 1) THEN
  fn='gasflux'
  ilength = 7
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,119)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1012)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,qxgas(jx,jy,jz),qygas(jx,jy,jz)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF


ELSE IF (ny == 1 .AND. nz > 1) THEN

IF (isaturate == 1) THEN
!!  fn='totgas'
  fn='gas'
  ilength = 3
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,123)

  WRITE(8,2286) (namg(kk),kk=1,ngas)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',ny
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
!!        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(sgas(i,jx,jy,jz),i = 1,ncomp)
        tk = 273.15d0 + t(jx,jy,jz)
        denmol = 1.e05/(8.314*tk)                      ! P/RT = n/V, with pressure converted from bars to Pascals
        CALL GasPartialPressure(ncomp,ngas,gastmp10,jx,jy,jz)
        WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,(gastmp10(kk),kk = 1,ngas)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF


  IF (ikph /= 0) THEN
    fn='pH'
    ilength = 2
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
    WRITE(8,*) 'VARIABLES = " X (meters)"," Z (meters)"," pH"'
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
    jy = 1
    DO jz = 1,nz
      DO jx = 1,nx
        phprt =  -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
        WRITE(8,183) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,phprt
      END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')
  END IF

!  Write out all of the species concentrations

  DO ik = 1,ncomp+nspec
    CALL stringlen(ulab(ik),len_sp(ik))
  END DO

!!  Individual  species concentrations (log units)

  fn='conc'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,101)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2009) (ulab(ik),ik=1,ncomp+nspec)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
    jy = 1
    DO jz = 1,nz
      DO jx = 1,nx
!fp! if_onproc({#expr# sp(ik,jx,jy,jz) #});
      WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(sp(IK,jx,jy,jz)/clg,IK = 1,ncomp+nspec)
!fp! end_onproc();
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!! Total concentrations

  fn='totcon'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,102)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2009) (ulab(ik),ik=1,ncomp)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
    DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
!fp! if_onproc({#expr# s(i,jx,jy,jz) #});
        WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(s(i,jx,jy,jz),i = 1,ncomp)
!fp! end_onproc();
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

IF (nexchange > 0) THEN
  fn='exchange'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2009) (nam_exchsec(nex),nex=1,nexch_sec)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
!fp! if_onproc({#expr# spex10(nex,jx,jy,jz) #});
        AqueousToBulk = ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
        WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(spex10(nex+nexchange,jx,jy,jz)/SolidSolutionRatioTemp,nex = 1,nexch_sec)
!fp! end_onproc();
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nexchange > 0) THEN
  fn='totexchange'
  ilength = 11
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2009) (ulab(i),i=1,ncomp)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        AqueousToBulk = ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
        totex_bas = 0.0
        DO i = 1,ncomp  
          DO nex = 1,nexch_sec
            totex_bas(i) = totex_bas(i) + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
          END DO
        END DO
        WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nsurf > 0) THEN
  fn='surface'
  ilength = 7
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2009) ( namsurf(is),is=1,nsurf ),(namsurf_sec(ns),ns=1,nsurf_sec)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
!fp! if_onproc({#expr# spsurf10(ns,jx,jy,jz) #});
        AqueousToBulk = ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
        WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(spsurf10(ns,jx,jy,jz),ns = 1,nsurf+nsurf_sec)
!fp! end_onproc();
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nsurf>0) THEN
  fn='totsurface'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2009) (ulab(i),i=1,ncomp)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
!fp! if_onproc({#expr# spsurf10(ns,jx,jy,jz) #});
        AqueousToBulk = ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
        totex_bas = 0.0
        DO i = 1,ncomp  
          DO ns = 1,nsurf_sec
            totex_bas(i) = totex_bas(i) + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
          END DO
        END DO
        WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
!fp! end_onproc();
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nrct > 0) THEN
  fn='TotMineral'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,105)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2009) (ulab(i),i=1,ncomp)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
!fp! if_onproc({#expr# volfx(k,jx,jy,jz) #});
        AqueousToBulk = ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
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
        WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
!fp! end_onproc();
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!  Write out the reaction rates in units of mol/L(bulk vol.)/sec

  DO k = 1,nrct
    CALL stringlen(umin(k),len_min(k))
  END DO

IF (nrct > 0) THEN
  fn='rate'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,108)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2011) (umin(k),k=1,nrct)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
!fp! if_onproc({#expr# dppt(k,jx,jy,jz) #});
      sum = 0.0
      DO k = 1,nrct
!************************
!  For units of volume %/year, uncomment the following line and
!  recompile
!              dptprt(k) = dppt(k,jx,jy,jz)*volmol(k)*100.0  ! volume %/yr
!***********************
!************************
!  For units of mol/L(BV)/sec, uncomment the following line and
        dptprt(k) = dppt(k,jx,jy,jz)/(secyr*1000.0)    ! mol/L(BV)/sec
!*************************************
        sum = sum + dptprt(k)
      END DO
      porcalc = sum
      WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(dptprt(k),k=1,nrct)
!fp! end_onproc();
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!!! 1011 format('VARIABLES = " X (meters)",
!!!     &      " Y (meters)"',50(', "',a25,'"'))
  1011 FORMAT('VARIABLES = " X (meters)", "  Y (meters)  "',50(', "',A11,'"'))
  2011 FORMAT('VARIABLES = " X (meters)", "  Z (meters)  "',50(', "',A11,'"'))

!   Write out the reaction rates in units of mol/L(bulk vol.)/sec

IF (ikin > 0) THEN
  fn='AqRate'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,109)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1011)  (namkin(ir),ir=1,ikin)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
!fp! if_onproc({#expr# raq_tot(ir,jx,jy,jz) #});
    WRITE(8,185) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(raq_tot(ir,jx,jy,jz),ir=1,ikin)
!fp! end_onproc();
    END DO
  END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!  Mineral volumes

IF (nrct > 0) THEN
  fn = 'volume'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,110)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1011) (umin(k),k=1,nrct)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      DO k = 1,nrct
        dvolpr(k) = volfx(k,jx,jy,jz)*1.0
      END DO
      WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(dvolpr(k),k=1,nrct)
    END DO
  END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!  Write out the porosity

  fn = 'porosity'

  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,112)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1010)
    WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      porprt = por(jx,jy,jz)*100.0
      WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,porprt
    END DO
  END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')    


!  Write out the saturation indices of the minerals (log Q/K).

IF (nrct > 0) THEN
  fn='saturation'
  ilength = 5
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,113)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,1011) (umin(k),k=1,nrct)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
!fp! if_onproc({#expr# sp(1,jx,jy,jz) #});
      CALL satcalc(ncomp,nrct,jx,jy,jz)
      DO k = 1,nrct
        dsat(k) = silog(1,k)
      END DO
      WRITE(8,184) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,(dsat(k),k=1,nrct)
!fp! end_onproc();
    END DO
  END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!  Write out Darcy fluxes

  fn='velocity'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,116)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2012)
  WRITE(8,*) 'ZONE F=POINT,I=', nx,  ', J=',nz
  DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      WRITE(8,191) x(jx)*OutputDistanceScale,z(jz)*OutputDistanceScale,qx(jx,jy,jz),qz(jx,jy,jz)
    END DO
  END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

ELSE   !  XY plot          1D case

!! One-dimensional case

IF (isaturate == 1) THEN
!!  fn='totgas'
  fn='gas'
  ilength = 3
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) (namg(kk),kk=1,ngas)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        tk = 273.15d0 + t(jx,jy,jz)
        denmol = 1.e05/(8.314*tk)                      ! P/RT = n/V, with pressure converted from bars to Pascals
        CALL GasPartialPressure(ncomp,ngas,gastmp10,jx,jy,jz)
        WRITE(8,184) x(jx)*OutputDistanceScale,(gastmp10(kk),kk = 1,ngas)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

  IF (O2Found .and. npointH2gas /= 0) THEN     !  Calculate pe based on O2-H2O couple--O2found refers to O2(aq)

    fn='pe'
    ilength = 2
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
    WRITE(8,*) 'VARIABLES = "X (meters)", "pe "'
    WRITE(8,*) 'ZONE F=POINT,I=', nx
    jz = 1
    DO jy = 1,ny
      DO jx = 1,nx
        pHprint = -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
        IF (ikh2o /= 0) THEN
          peprint = (-0.5d0*keqgas(npointH2gas,jx,jy,jz) +             &
            0.25d0*( sp(npointO2aq,jx,jy,jz)+gam(npointO2aq,jx,jy,jz)) - 0.5*gam(ikh2o,jx,jy,jz) -     &  
            pHprint*clg + 0.25d0*keqgas(npointO2gas,jx,jy,jz) )/clg
!!          Ehprint = peprint*clg*rgas*(t(jx,jy,jz)+273.15)/23.06
        ELSE
          peprint = (-0.5d0*keqgas(npointH2gas,jx,jy,jz) +             &
            0.25d0*( sp(npointO2aq,jx,jy,jz)+gam(npointO2aq,jx,jy,jz)) -     &  
            pHprint*clg + 0.25d0*keqgas(npointO2gas,jx,jy,jz) )/clg
        END IF
        WRITE(8,184) x(jx)*OutputDistanceScale,peprint 

      END DO
    END DO

    CLOSE(UNIT=8,STATUS='keep')

  ELSE IF (ikfe2 /= 0 .and. ikfe3 /= 0) THEN

    fn='pe'
    ilength = 2
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
    WRITE(8,*) 'VARIABLES = "X (meters)", "pe "'
    WRITE(8,*) 'ZONE F=POINT,I=', nx
    jz = 1
    DO jy = 1,ny
      DO jx = 1,nx

        fe2print = ( sp(ikfe2,jx,jy,jz)+gam(ikfe2,jx,jy,jz) )/clg
        fe3print = ( sp(ikfe3,jx,jy,jz)+gam(ikfe3,jx,jy,jz) )/clg
        peprint = 13.02825 + fe3print - fe2print 
        WRITE(8,184) x(jx)*OutputDistanceScale,peprint 

      END DO
    END DO

  ELSE

    CONTINUE

  END IF

!! Total concentrations

  fn='totcon'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,102)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) (ulab(ik),ik=1,ncomp)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(8,185) x(jx)*OutputDistanceScale,(s(i,jx,jy,jz),i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  IF (ikph /= 0) THEN
    fn='pH'
    ilength = 2
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
    WRITE(8,*) 'VARIABLES = "X (meters)", "pH "'
    WRITE(8,*) 'ZONE F=POINT,I=', nx
    jz = 1
    DO jy = 1,ny
      DO jx = 1,nx
        phprt =  -(sp(ikph,jx,jy,jz)+gam(ikph,jx,jy,jz))/clg
        WRITE(8,185) x(jx)*OutputDistanceScale,phprt
      END DO
    END DO
    CLOSE(UNIT=8,STATUS='keep')
  END IF

!  Write out all of the species concentrations

  DO ik = 1,ncomp+nspec
    CALL stringlen(ulab(ik),len_sp(ik))
  END DO

!!  Individual  species concentrations (log units)

  fn='conc'
  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,101)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) (ulab(ik),ik=1,ncomp+nspec)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      WRITE(8,185) x(jx)*OutputDistanceScale,(sp(IK,jx,jy,jz)/clg,IK = 1,ncomp+nspec)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

IF (nexchange > 0) THEN
  fn='exchange'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) (nam_exchsec(nex),nex=1,nexch_sec)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        AqueousToBulk = ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
        WRITE(8,185) x(jx)*OutputDistanceScale,(spex10(nex+nexchange,jx,jy,jz),nex = 1,nexch_sec)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='totexchange'
  ilength = 11
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) (ulab(i),i=1,ncomp)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        AqueousToBulk = ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
        totex_bas = 0.0
        DO i = 1,ncomp  
          DO nex = 1,nexch_sec
            totex_bas(i) = totex_bas(i) + muexc(nex,i)*spex10(nex+nexchange,jx,jy,jz)
          END DO
        END DO
        WRITE(8,185) x(jx)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nsurf > 0) THEN
  fn='surface'
  ilength = 7
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) ( namsurf(is),is=1,nsurf ),(namsurf_sec(ns),ns=1,nsurf_sec)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        AqueousToBulk = ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
        WRITE(8,185) x(jx)*OutputDistanceScale,(spsurf10(ns,jx,jy,jz),ns = 1,nsurf+nsurf_sec)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

  fn='totsurface'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,104)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) (ulab(i),i=1,ncomp)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        AqueousToBulk = ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
        totex_bas = 0.0
        DO i = 1,ncomp  
          DO ns = 1,nsurf_sec
            totex_bas(i) = totex_bas(i) + musurf(ns,i)*spsurf10(ns+nsurf,jx,jy,jz)
          END DO
        END DO
        WRITE(8,185) x(jx)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

IF (nrct > 0) THEN
  fn='TotMineral'
  ilength = 10
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,105)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) (ulab(i),i=1,ncomp)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        AqueousToBulk = ro(jx,jy,jz)*satliq(jx,jy,jz)*por(jx,jy,jz)
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
        WRITE(8,185) x(jx)*OutputDistanceScale,(totex_bas(i),i = 1,ncomp)
      END DO
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')

!  Write out the reaction rates in units of mol/L(bulk vol.)/sec

  DO k = 1,nrct
    CALL stringlen(umin(k),len_min(k))
  END DO

  fn='rate'

  ilength = 4
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,108)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) (umin(k),k=1,nrct)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      sum = 0.0
      DO k = 1,nrct
!************************
!  For units of volume %/year, uncomment the following line and
!  recompile
!              dptprt(k) = dppt(k,jx,jy,jz)*volmol(k)*100.0  ! volume %/yr
!***********************
!************************
!  For units of mol/L(BV)/sec, uncomment the following line and
        dptprt(k) = dppt(k,jx,jy,jz)/(secyr*1000.0)    ! mol/L(BV)/sec
!*************************************
        sum = sum + dptprt(k)
      END DO
      porcalc = sum
      WRITE(8,185) x(jx)*OutputDistanceScale,(dptprt(k),k=1,nrct)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!   Write out the reaction rates in units of mol/L(bulk vol.)/sec

IF (ikin > 0) THEN
  fn='AqRate'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,109)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001)  (namkin(ir),ir=1,ikin)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  jy = 1
  jz = 1
  DO jx = 1,nx
    WRITE(8,185) x(jx)*OutputDistanceScale,(raq_tot(ir,jx,jy,jz),ir=1,ikin)
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!  Volumes %

IF (nrct > 0) THEN
  fn = 'volume'
  ilength = 6
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,110)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) (umin(k),k=1,nrct)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      DO k = 1,nrct
        dvolpr(k) = volfx(k,jx,jy,jz)*1.0
      END DO
      WRITE(8,185) x(jx)*OutputDistanceScale,(dvolpr(k),k=1,nrct)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

!  Write out the porosity

IF (ny == 1) THEN
  fn = 'porosity'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,112)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,*) 'VARIABLES = "X (meters)", "Porosity "'
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
!!      porprt = por(jx,jy,jz)*100.0
      porprt = por(jx,jy,jz)*1.0
      WRITE(8,184) x(jx)*OutputDistanceScale,porprt
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
  
ELSE
  fn = 'porosity'
  ilength = 8
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,112)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,*) 'VARIABLES = "X (meters)", "Porosity "'
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      porprt = por(jx,jy,jz)*1.0
      WRITE(8,184) x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale,porprt
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')    
END IF

!  Write out pressure

  IF (calculateflow) THEN
    fn='pressure'
    ilength = 8
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,115)
    WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
    WRITE(8,1015)
    WRITE(8,*) 'ZONE F=POINT,I=', nx
    jy = 1
    jz = 1
    DO jx = 1,nx
      WRITE(8,184) x(jx)*OutputDistanceScale,pres(jx,jy,jz)
    END DO
    CLOSE(UNIT=8,STATUS='keep')
    
    fn='velocity'
    ilength = 8
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,116)
    WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
    WRITE(8,2015)
    WRITE(8,*) 'ZONE F=POINT,I=', nx
    DO jx = 1,nx
      WRITE(8,191) x(jx)*OutputDistanceScale,qx(jx,jy,jz)
    END DO
    CLOSE(UNIT=8,STATUS='keep')
    
    fn = 'permeability'
    ilength = 12
    CALL newfile(fn,suf1,fnv,nint,ilength)
    OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
    WRITE(8,118)
    WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
    WRITE(8,1021)
    WRITE(8,*) 'ZONE F=POINT,I=', nx
    jz = 1
    jy = 1
      DO jx = 1,nx
        WRITE(8,184) x(jx)*OutputDistanceScale,Log10(permx(jx,jy,jz))
      END DO
    CLOSE(UNIT=8,STATUS='keep')
    
  END IF

!  Write out the saturation indices of the minerals (log Q/K).
 
IF (nrct > 0) THEN
  fn='saturation'
  ilength = 5
  CALL newfile(fn,suf1,fnv,nint,ilength)
  OPEN(UNIT=8,FILE=fnv, ACCESS='sequential',STATUS='unknown')
  WRITE(8,113)
  WRITE(8,*) 'TITLE = "',char_time(1:ls),' Years"'
  WRITE(8,2001) (umin(k),k=1,nrct)
  WRITE(8,*) 'ZONE F=POINT,I=', nx
  jz = 1
  DO jy = 1,ny
    DO jx = 1,nx
      CALL satcalc(ncomp,nrct,jx,jy,jz)
      DO k = 1,nrct
        dsat(k) = silog(1,k)
      END DO
      WRITE(8,185) x(jx)*OutputDistanceScale,(dsat(k),k=1,nrct)
    END DO
  END DO
  CLOSE(UNIT=8,STATUS='keep')
END IF

185 FORMAT(1PE12.5,12x,100(1X,1PE16.8))

END IF

1009 FORMAT('VARIABLES = " X (meters) ", "  Y (meters)  "',100(', "',A13,'"') )
2009 FORMAT('VARIABLES = " X (meters) ", "  Z (meters)  "',100(', "',A13,'"'))
2001 FORMAT('VARIABLES = "X (meters) "',                   100(', "',A13,'"'))

1012 FORMAT('VARIABLES = " X (meters)", " Y (meters)", "X Velocity", "Y Velocity"')
2012 FORMAT('VARIABLES = " X (meters)", " Z (meters)", "X Velocity", "Z Velocity"')
1013 FORMAT('VARIABLES = " X (meters)", " Y (meters)", "Pressure"')
1015 FORMAT('VARIABLES = " X (meters)", "Pressure"')
2014 FORMAT('VARIABLES = " X (meters)", " Y (meters)", "Tortuosity"')
2015 FORMAT('VARIABLES = " X (meters)", "X_Velocity"')
     
1021 FORMAT('VARIABLES = " X (meters)", "X-Permeability"')

182 FORMAT(100(1X,1PE12.4))
183 FORMAT(1PE12.4,2X,1PE12.4,2X,1PE12.4)
184 FORMAT(100(1X,1PE16.8))
191 FORMAT(100(1X,1PE17.7))
188 FORMAT(100(1X,f15.7))

2283 FORMAT('# Time (yrs) ',2X,1PE12.4)
2284 FORMAT('      X        ','     Y        ',a18)
2282 FORMAT('   X        ','     Y        ','        pH')
2281 FORMAT('   X        ','     Y        ',4X,a18)
2285 FORMAT('    X        ','     Y        ',3X,30(1X,a13))
2286 FORMAT('    X        ','      Y               ',3X,30(1X,a15))

RETURN
END SUBROUTINE GraphicsTecplot
!  *******************************************************
