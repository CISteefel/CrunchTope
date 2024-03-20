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


SUBROUTINE StartTope(ncomp,nspec,nkin,nrct,ngas,npot,                   &
  nx,ny,nz,data1,ipath,igamma,ikmast,ikph,iko2,ltitle,    &
  tstep,delt,deltmin,ttol,jpor,ikin,nstop,                          &
  corrmax,nseries,minseries,nexchange,nexch_sec,nsurf,nsurf_sec,ndecay,       &
  str_mon,str_day,str_hr,str_min,str_sec,NumInputFiles,             &
  InputFileCounter,nBoundaryConditionZone)
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
USE io
USE strings
USE ReadFlow
USE modflowModule
!!!USE isotope, ONLY: IsotopeMineralRare, IsotopeMineralCommon,IsotopePrimaryCommon,IsotopePrimaryRare
USE NanoCrystal
USE isotope

IMPLICIT NONE

!  ****************  Beginning of interface blocks  ***************************

INTERFACE
SUBROUTINE readModFlowStress(nout,perlen,nstp,tsmult,ndimdummy)
USE crunchtype
IMPLICIT NONE
INTEGER(I4B), INTENT(IN)                                    :: nout
REAL(DP), DIMENSION(:), INTENT(OUT)                         :: perlen
INTEGER(I4B), DIMENSION(:), INTENT(OUT)                     :: nstp
REAL(DP), DIMENSION(:), INTENT(OUT)                         :: tsmult
INTEGER(I4B), INTENT(IN)                                    :: ndimdummy
END SUBROUTINE readModFlowStress
END INTERFACE

INTERFACE
SUBROUTINE find_condition(nin,nout,found,phfound,ncomp,  &
  nspec,nrct,nkin,ngas,nexchange,nsurf,ndecay,           &
  ph,guessph,constraint,nchem,unitsflag,jpor,            &
  DensityModule,RunningPest)
USE crunchtype
USE params
USE concentration
USE mineral
USE medium
USE temperature
IMPLICIT NONE
INTEGER(I4B), INTENT(IN)                                     :: nin
INTEGER(I4B), INTENT(IN)                                     :: nout
INTEGER(I4B), INTENT(IN)                                     :: ncomp
INTEGER(I4B), INTENT(IN)                                     :: nspec
INTEGER(I4B), INTENT(IN)                                     :: nrct
INTEGER(I4B), INTENT(IN)                                     :: nkin
INTEGER(I4B), INTENT(IN)                                     :: ngas
INTEGER(I4B), INTENT(IN)                                     :: nexchange
INTEGER(I4B), INTENT(IN)                                     :: nsurf
INTEGER(I4B), INTENT(IN)                                     :: ndecay
INTEGER(I4B), INTENT(OUT)                                    :: nchem
INTEGER(I4B), INTENT(IN)                                     :: jpor
LOGICAL(LGT), INTENT(IN OUT)                                 :: found
LOGICAL(LGT), INTENT(IN OUT)                                 :: phfound
REAL(DP), DIMENSION(:), INTENT(OUT)                          :: ph
REAL(DP), DIMENSION(:), INTENT(OUT)                          :: guessph
CHARACTER (LEN=mls), DIMENSION(:,:), INTENT(IN OUT)          :: constraint
INTEGER(I4B), DIMENSION(:), INTENT(INOUT)                    :: unitsflag
CHARACTER (LEN=mls),  INTENT(IN)                             :: DensityModule
LOGICAL(LGT), INTENT(IN)                                     :: RunningPest
END SUBROUTINE find_condition
END INTERFACE

INTERFACE
SUBROUTINE read_multpar(nout,lchar,parchar,parfind,realmult,lenarray,section)
USE crunchtype
USE params
USE strings
IMPLICIT NONE
INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(OUT)                                   :: lchar
CHARACTER (LEN=mls), INTENT(IN)                             :: parchar
CHARACTER (LEN=mls), INTENT(IN OUT)                         :: parfind
INTEGER(I4B), INTENT(OUT)                                   :: lenarray
REAl(DP), DIMENSION(:), INTENT(IN OUT)                      :: realmult
CHARACTER (LEN=mls), INTENT(IN)                             :: section
END SUBROUTINE read_multpar
END INTERFACE

INTERFACE
SUBROUTINE read_snapshot(nout,lchar,parchar,parchar2,parfind,realmult,lenarray,section)
USE crunchtype
USE params
USE strings
IMPLICIT NONE
INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(OUT)                                   :: lchar
CHARACTER (LEN=mls), INTENT(IN)                             :: parchar
CHARACTER (LEN=mls), INTENT(IN)                             :: parchar2
CHARACTER (LEN=mls), INTENT(IN OUT)                         :: parfind
INTEGER(I4B), INTENT(OUT)                                   :: lenarray
REAl(DP), DIMENSION(:), INTENT(IN OUT)                      :: realmult
CHARACTER (LEN=mls), INTENT(IN)                             :: section
END SUBROUTINE read_snapshot
END INTERFACE


INTERFACE
SUBROUTINE read_multstring(nout,lchar,parchar,parfind, stringarray,lenarray,section)
USE crunchtype
USE params
USE strings
IMPLICIT NONE
INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(OUT)                                   :: lchar
CHARACTER (LEN=mls), INTENT(IN)                             :: parchar
CHARACTER (LEN=mls), INTENT(IN OUT)                         :: parfind
CHARACTER (LEN=mls), DIMENSION(:), INTENT(IN OUT)           :: stringarray
INTEGER(I4B), INTENT(OUT)                                   :: lenarray
CHARACTER (LEN=mls), INTENT(IN)                             :: section
END SUBROUTINE read_multstring
END INTERFACE

INTERFACE
SUBROUTINE ModScan(nx,ny,nz,cnhIn,wellIn,riverIn,ndimdummy,mxwell,mxrivr,mxdrn)
USE crunchtype
IMPLICIT NONE
INTEGER(I4B), INTENT(IN)                                 :: nx
INTEGER(I4B), INTENT(IN)                                 :: ny
INTEGER(I4B), INTENT(IN)                                 :: nz
LOGICAL(LGT), DIMENSION(:), INTENT(OUT)                  :: cnhIn
LOGICAL(LGT), DIMENSION(:), INTENT(OUT)                  :: riverIn
LOGICAL(LGT), DIMENSION(:), INTENT(OUT)                  :: wellIn
INTEGER(I4B), INTENT(IN)                                 :: ndimdummy
INTEGER(I4B), INTENT(IN)                                 :: mxwell
INTEGER(I4B), INTENT(IN)                                 :: mxrivr
INTEGER(I4B), INTENT(IN)                                 :: mxdrn
END SUBROUTINE ModScan
END INTERFACE

INTERFACE
SUBROUTINE readModFlowParameters(nout,nchem,nparams,jxTemp,jyTemp,jzTemp,  &
            conditionNum,modflowstring,lenstring)
USE crunchtype
USE strings
IMPLICIT NONE
INTEGER(I4B), INTENT(IN)                                    :: nout
INTEGER(I4B), INTENT(IN)                                    :: nchem
INTEGER(I4B), INTENT(OUT)                                   :: nparams
INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jxTemp
INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jyTemp
INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: jzTemp
INTEGER(I4B), DIMENSION(:),INTENT(OUT)                      :: conditionNum
CHARACTER (LEN=mls), INTENT(IN)                             :: modflowstring
INTEGER(I4B), INTENT(IN)                                    :: lenstring
END SUBROUTINE readModFlowParameters
END INTERFACE

INTERFACE
SUBROUTINE GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,xxx,TEXT)
USE crunchtype
IMPLICIT NONE
INTEGER(I4B), INTENT(IN)                                    :: nx
INTEGER(I4B), INTENT(IN)                                    :: ny
INTEGER(I4B), INTENT(IN)                                    :: nz
INTEGER(I4B), INTENT(IN)                                    :: lowX
INTEGER(I4B), INTENT(IN)                                    :: lowY
INTEGER(I4B), INTENT(IN)                                    :: lowZ
INTEGER(I4B), INTENT(IN)                                    :: highX
INTEGER(I4B), INTENT(IN)                                    :: highY
INTEGER(I4B), INTENT(IN)                                    :: highZ
REAL(DP), DIMENSION(lowX:highX,lowY:highY,lowZ:highZ), INTENT(INOUT)                   :: xxx
CHARACTER (LEN=15), INTENT(IN)                              :: text
END SUBROUTINE GhostCells
END INTERFACE

!  ****************  End of interface blocks  ***************************

!  External variables and arrays

INTEGER(I4B), INTENT(OUT)                                     :: ncomp
INTEGER(I4B), INTENT(OUT)                                     :: nspec
INTEGER(I4B), INTENT(OUT)                                     :: nkin
INTEGER(I4B), INTENT(OUT)                                     :: nrct
INTEGER(I4B), INTENT(OUT)                                     :: ngas
INTEGER(I4B), INTENT(OUT)                                     :: npot
INTEGER(I4B), INTENT(OUT)                                     :: nx
INTEGER(I4B), INTENT(OUT)                                     :: ny
INTEGER(I4B), INTENT(OUT)                                     :: nz
INTEGER(I4B), INTENT(OUT)                                     :: ipath
INTEGER(I4B), INTENT(OUT)                                     :: igamma
INTEGER(I4B), INTENT(OUT)                                     :: ikmast
INTEGER(I4B), INTENT(OUT)                                     :: ikph
INTEGER(I4B), INTENT(OUT)                                     :: ikO2
INTEGER(I4B), INTENT(IN OUT)                                  :: jpor
INTEGER(I4B), INTENT(OUT)                                     :: ikin
INTEGER(I4B), INTENT(OUT)                                     :: nstop
INTEGER(I4B), INTENT(OUT)                                     :: nseries
INTEGER(I4B), INTENT(OUT)                                     :: minseries
INTEGER(I4B), INTENT(OUT)                                      :: nexchange
INTEGER(I4B), INTENT(OUT)                                      :: nexch_sec
INTEGER(I4B), INTENT(OUT)                                      :: nsurf
INTEGER(I4B), INTENT(OUT)                                      :: nsurf_sec
INTEGER(I4B), INTENT(OUT)                                      :: ndecay
INTEGER(I4B), INTENT(IN OUT)                                  :: NumInputFiles
INTEGER(I4B), INTENT(IN OUT)                                  :: InputFileCounter


CHARACTER (LEN=mls), INTENT(IN OUT)                           :: data1
CHARACTER (LEN=2*mls), INTENT(IN OUT)                         :: ltitle

REAL(DP), INTENT(OUT)                                         :: tstep
REAL(DP), INTENT(OUT)                                         :: delt
REAL(DP), INTENT(OUT)                                         :: deltmin
REAL(DP), INTENT(OUT)                                         :: ttol
REAL(DP), INTENT(OUT)                                         :: corrmax

!  Internal variables and arrays

INTEGER(I4B), PARAMETER                                       :: ndimdummy=500

REAL(DP), DIMENSION(:), ALLOCATABLE                           :: work1
REAL(DP), DIMENSION(:,:), ALLOCATABLE                         :: work2
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                       :: work3
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                       :: work3b
REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE                     :: work4

INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: workint1
INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE                     :: workint2
INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE                   :: workint3

CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: workchar1
CHARACTER (LEN=mls), DIMENSION(:,:), ALLOCATABLE              :: workchar2
CHARACTER (LEN=mls), DIMENSION(:,:,:), ALLOCATABLE            :: workchar3

INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: unitsflag

INTEGER(I4B)                                                  :: str_mon
INTEGER(I4B)                                                  :: str_day
INTEGER(I4B)                                                  :: str_hr
INTEGER(I4B)                                                  :: str_min
INTEGER(I4B)                                                  :: str_sec
INTEGER(I4B)                                                  :: str_millisec
INTEGER(I4B)                                                  :: lenInput

CHARACTER (LEN=15)                                            :: text
CHARACTER (LEN=mls)                                           :: filename
CHARACTER (LEN=mls)                                           :: FileOutput
CHARACTER (LEN=mls)                                           :: vxfile
CHARACTER (LEN=mls)                                           :: vyfile
CHARACTER (LEN=mls)                                           :: vzfile
CHARACTER (LEN=mls)                                           :: permxfile
CHARACTER (LEN=mls)                                           :: permyfile
CHARACTER (LEN=mls)                                           :: permzfile
CHARACTER (LEN=mls)                                           :: dummy
CHARACTER (LEN=mls)                                           :: dumstring
CHARACTER (LEN=mls)                                           :: section
CHARACTER (LEN=mls)                                           :: parchar
CHARACTER (LEN=mls)                                           :: parchar2
CHARACTER (LEN=mls)                                           :: parfind
CHARACTER (LEN=mls)                                           :: velocityfile
CHARACTER (LEN=mls)                                           :: TortuosityFile
CHARACTER (LEN=mls)                                           :: gasvelocityfile
CHARACTER (LEN=mls)                                           :: permfile
CHARACTER (LEN=mls)                                           :: vgnfile
CHARACTER (LEN=mls)                                           :: vgafile
CHARACTER (LEN=mls)                                           :: wcrfile
CHARACTER (LEN=mls)                                           :: breakfile
CHARACTER (LEN=mls)                                           :: namtemp
CHARACTER (LEN=12)                                            :: writeph
CHARACTER (LEN=mls)                                           :: PorosityFile
CHARACTER (LEN=mls)                                           :: GridVolumeFile
CHARACTER (LEN=mls)                                           :: SaturationFile
CHARACTER (LEN=mls)                                           :: BurialFile
CHARACTER (LEN=mls)                                           :: ConditionName

INTEGER(I4B)                                                  :: ConditionNumber

LOGICAL(LGT)                                                  :: CheckDuan
LOGICAL(LGT)                                                  :: found
LOGICAL(LGT)                                                  :: activity_dbh
LOGICAL(LGT)                                                  :: lag_activity
LOGICAL(LGT)                                                  :: database_sweep
LOGICAL(LGT)                                                  :: porosity_update
LOGICAL(LGT)                                                  :: phfound,speciesfound
LOGICAL(LGT)                                                  :: reaction_path
LOGICAL(LGT)                                                  :: constant_flow
LOGICAL(LGT)                                                  :: constant_gasflow
LOGICAL(LGT)                                                  :: readvelocity
LOGICAL(LGT)                                                  :: readgasvelocity
LOGICAL(LGT)                                                  :: readperm
LOGICAL(LGT)                                                  :: readvgn
LOGICAL(LGT)                                                  :: readvga
LOGICAL(LGT)                                                  :: readwcr
LOGICAL(LGT)                                                  :: readpermx
LOGICAL(LGT)                                                  :: readpermy
LOGICAL(LGT)                                                  :: onlyspeciate
LOGICAL(LGT)                                                  :: genericrates
LOGICAL(LGT)                                                  :: DaughterFound
LOGICAL(LGT)                                                  :: SolveHindmarsh
LOGICAL(LGT)                                                  :: ext
LOGICAL(LGT)                                                  :: pest
LOGICAL(LGT)                                                  :: readburial
LOGICAL(LGT)                                                  :: NoFluidBury
LOGICAL(LGT)                                                  :: ReadTortuosity

LOGICAL(LGT)                                                  :: streamtube
CHARACTER (LEN=mls), DIMENSION(:,:), ALLOCATABLE              :: constraint
REAL(DP), DIMENSION(:), ALLOCATABLE                           :: realmult
REAL(DP), DIMENSION(:), ALLOCATABLE                           :: pH
REAL(DP), DIMENSION(:), ALLOCATABLE                           :: guessph

INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: list_tmp
INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE                     :: krad_tmp
INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE                     :: nprad_tmp

!  MODFLOW Read parameters

LOGICAL(LGT)                                                  :: GMSsecondary
LOGICAL(LGT)                                                  :: GMSmineral

!!LOGICAL(LGT)                                                  :: ReadGeochemicalConditions
INTEGER(I4B)                                                  :: GridCoordinateX
INTEGER(I4B)                                                  :: GridCoordinateY

INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: jxTemp
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: jyTemp
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: jzTemp
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: conditionNum
INTEGER(I4B)                                                  :: nparams
INTEGER(I4B)                                                  :: lenstring
INTEGER(I4B)                                                  :: mxwell
INTEGER(I4B)                                                  :: mxrivr
INTEGER(I4B)                                                  :: mxdrn
INTEGER(I4B)                                                  :: NeedWellCondition
INTEGER(I4B)                                                  :: NeedRiverCondition
INTEGER(I4B)                                                  :: NeedHeadCondition
INTEGER(I4B)                                                  :: ModFlowTimeUnits

CHARACTER (LEN=mls)                                           :: modflowstring
CHARACTER (LEN=mls)                                           :: SaturationFileFormat
CHARACTER (LEN=mls)                                           :: PorosityFileFormat
CHARACTER (LEN=mls)                                           :: GridVolumeFileFormat
CHARACTER (LEN=mls)                                           :: BurialFileFormat
CHARACTER (LEN=mls)                                           :: TortuosityFileFormat
CHARACTER (LEN=mls)                                           :: GasVelocityFileFormat
CHARACTER (LEN=mls)                                           :: PermFileFormat
CHARACTER (LEN=mls)                                           :: vgnFileFormat
CHARACTER (LEN=mls)                                           :: vgaFileFormat
CHARACTER (LEN=mls)                                           :: wcrFileFormat
CHARACTER (LEN=mls)                                           :: permxFileFormat
CHARACTER (LEN=mls)                                           :: permyFileFormat
CHARACTER (LEN=mls)                                           :: VelocityFileFormat
CHARACTER (LEN=mls)                                           :: TemperatureFileFormat
CHARACTER (LEN=mls)                                           :: FileTemp

INTEGER(I4B)                                                  :: FileNameLength

REAL(DP)                                                      :: dtBase
REAL(DP)                                                      :: geometricSum
REAL(DP)                                                      :: checkSum
INTEGER(I4B)                                                  :: ntt
REAL(DP), DIMENSION(:), ALLOCATABLE                           :: dtTemp
REAL(DP), DIMENSION(:), ALLOCATABLE                           :: perlen
REAL(DP), DIMENSION(:), ALLOCATABLE                           :: tsmult
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                           :: CECconvert
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: nstp

LOGICAL(LGT), DIMENSION(:), ALLOCATABLE                       :: cnhIn
LOGICAL(LGT), DIMENSION(:), ALLOCATABLE                       :: riverIn
LOGICAL(LGT), DIMENSION(:), ALLOCATABLE                       :: wellIn

REAL(DP), PARAMETER                                           :: eps=1.e-12

REAL(DP)                                                      :: xdum
REAL(DP)                                                      :: ydum
REAL(DP)                                                      :: zdum
REAL(DP)                                                      :: permdum
REAL(DP)                                                      :: time_scale
REAL(DP)                                                      :: dist_scale
REAL(DP)                                                      :: realjunk
REAL(DP)                                                      :: time
REAL(DP)                                                      :: sumpor
REAL(DP)                                                      :: portemp
REAL(DP)                                                      :: portemp1
REAL(DP)                                                      :: PressureTemp
REAL(DP)                                                    :: term1
REAL(DP)                                                    :: term2
REAL(DP)                                                    :: termA
REAL(DP)                                                    :: termB
REAL(DP)                                                    :: termPerturb
REAL(DP),PARAMETER                                                    :: smallNumber=1.0E-09
REAL(DP)                                                      :: tempc
REAL(DP)                                                      :: sum
REAL(DP)                                                      :: tgradprt
REAL(DP)                                                      :: qxinit
REAL(DP)                                                      :: qyinit
REAL(DP)                                                      :: qzinit
REAL(DP)                                                      :: qxgasinit
REAL(DP)                                                      :: qygasinit
REAL(DP)                                                      :: qzgasinit
REAL(DP)                                                      :: qxmax
REAL(DP)                                                      :: qymax
REAL(DP)                                                      :: qzmax
REAL(DP)                                                      :: dspxmax
REAL(DP)                                                      :: dspymax
REAL(DP)                                                      :: dspzmax
REAL(DP)                                                      :: qbar
REAL(DP)                                                      :: convert
REAL(DP)                                                      :: MinSaturation
REAL(DP)                                                      :: TortuosityConstant

REAL(DP)                                                      :: LogPermMinX
REAL(DP)                                                      :: LogPermMinY
REAL(DP)                                                      :: LogPermMaxX
REAL(DP)                                                      :: LogPermMaxY
REAL(DP)                                                      :: MidPointX
REAL(DP)                                                      :: MidPointY
REAL(DP)                                                      :: NewPermXRange
REAL(DP)                                                      :: NewPermYRange
REAL(DP)                                                      :: LogPermX
REAL(DP)                                                      :: LogPermY
REAL(DP)                                                      :: NewLogPermX
REAL(DP)                                                      :: NewLogPermY

INTEGER(I4B)                                                  :: nxtemp
INTEGER(I4B)                                                  :: nytemp
INTEGER(I4B)                                                  :: nztemp
INTEGER(I4B)                                                  :: i
INTEGER(I4B)                                                  :: k
INTEGER(I4B)                                                  :: npt
!!INTEGER(I4B)                                                  :: nchem
INTEGER(I4B)                                                  :: nin
INTEGER(I4B)                                                  :: nout
INTEGER(I4B)                                                  :: lchar
INTEGER(I4B)                                                  :: icomplete
INTEGER(I4B)                                                  :: ispeciate
INTEGER(I4B)                                                  :: ncount
INTEGER(I4B)                                                  :: igenericrates
INTEGER(I4B)                                                  :: np
INTEGER(I4B)                                                  :: kk
INTEGER(I4B)                                                  :: ksp
INTEGER(I4B)                                                  :: ik
INTEGER(I4B)                                                  :: nex
INTEGER(I4B)                                                  :: ncnt
INTEGER(I4B)                                                  :: ns
INTEGER(I4B)                                                  :: nco
INTEGER(I4B)                                                  :: ls
INTEGER(I4B)                                                  :: id
INTEGER(I4B)                                                  :: nisotope_max
INTEGER(I4B)                                                  :: nmindecay_max
INTEGER(I4B)                                                  :: kd
INTEGER(I4B)                                                  :: ndim1
INTEGER(I4B)                                                  :: ndim2
INTEGER(I4B)                                                  :: ndim3
INTEGER(I4B)                                                  :: ndim4
INTEGER(I4B)                                                  :: j
INTEGER(I4B)                                                  :: neqn
INTEGER(I4B)                                                  :: l
INTEGER(I4B)                                                  :: is
INTEGER(I4B)                                                  :: iinit
INTEGER(I4B)                                                  :: ix
INTEGER(I4B)                                                  :: ikn
INTEGER(I4B)                                                  :: nretard
INTEGER(I4B)                                                  :: nzonex
INTEGER(I4B)                                                  :: nzoney
INTEGER(I4B)                                                  :: nzonez
INTEGER(I4B)                                                  :: nxyz
INTEGER(I4B)                                                  :: ii
INTEGER(I4B)                                                  :: nsum
INTEGER(I4B)                                                  :: jx
INTEGER(I4B)                                                  :: jy
INTEGER(I4B)                                                  :: jz
INTEGER(I4B)                                                  :: jxx
INTEGER(I4B)                                                  :: jyy
INTEGER(I4B)                                                  :: jzz
INTEGER(I4B)                                                  :: lenarray
INTEGER(I4B)                                                  :: nbig
INTEGER(I4B)                                                  :: ntot
INTEGER(I4B)                                                  :: nhet
INTEGER(I4B)                                                  :: ll
INTEGER(I4B)                                                  :: i2
!!INTEGER(I4B)                                                  :: isotope
INTEGER(I4B)                                                  :: nbnd
INTEGER(I4B)                                                  :: nlength
INTEGER(I4B)                                                  :: intjunk
INTEGER(I4B)                                                  :: lfile
INTEGER(I4B)                                                  :: npermx
INTEGER(I4B)                                                  :: npermy
INTEGER(I4B)                                                  :: npermz
INTEGER(I4B)                                                  :: ngaspump
INTEGER(I4B)                                                  :: intfile
INTEGER(I4B)                                                  :: ndiff
INTEGER(I4B)                                                  :: itmp
INTEGER(I4B)                                                  :: ir
INTEGER(I4B)                                                  :: ir2
INTEGER(I4B)                                                  :: irsave
INTEGER(I4B)                                                  :: idaughter
INTEGER(I4B)                                                  :: ipest
INTEGER(I4B)                                                  :: PestUnit
INTEGER(I4B)                                                  :: npressure
INTEGER(I4B)                                                  :: lowX
INTEGER(I4B)                                                  :: lowY
INTEGER(I4B)                                                  :: lowZ
INTEGER(I4B)                                                  :: highX
INTEGER(I4B)                                                  :: highY
INTEGER(I4B)                                                  :: highZ


INTEGER(I4B)                                                  :: nplotsurface
INTEGER(I4B)                                                  :: nplotexchange
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: ExchangeBasis
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: SurfaceBasis

REAL(DP)                                                      :: MeanSaltConcentration
REAL(DP)                                                      :: MassFraction
REAL(DP)                                                      :: PauseLength
REAL(DP)                                                      :: pi
REAL(DP)                                                      :: erodex
REAL(DP)                                                      :: erodey
REAL(DP)                                                      :: RateGeneric
REAL(DP)                                                      :: denmol
REAL(DP)                                                      :: tk

REAL(DP)                                                      :: SumMineralVolume

REAL(DP)  :: dum1
REAL(DP)  :: dum2
REAL(DP)  :: PorosityRead
REAL(DP)  :: QuartzRead
REAL(DP)  :: ChloriteRead
REAL(DP)  :: IlliteRead
REAL(DP)  :: KaoliniteRead
REAL(DP)  :: SmectiteRead
REAL(DP)  :: FeOxideRead

CHARACTER (LEN=12)                                            :: dumm1
CHARACTER (LEN=12)                                            :: dumm2
CHARACTER (LEN=12)                                            :: dumm3
INTEGER(I4B), DIMENSION(8)                                    :: curr_time

CHARACTER (LEN=mls)                                           :: data2
CHARACTER (LEN=mls)                                           :: data3
CHARACTER (LEN=mls)                                           :: NameMineral
CHARACTER (LEN=mls)                                           :: label
CHARACTER (LEN=mls)                                           :: Surface

CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: NucleationMineral
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: NameNucleationPathway


REAL(DP)                                                      :: A_zero25C
REAL(DP)                                                      :: B_nucleation
REAL(DP)                                                      :: Sigma_mJm2
REAL(DP)                                                      :: SSA_m2g

REAL(DP)                                                      :: ScaleMineralVolumes


INTEGER(I4B)                                                  :: knucl
INTEGER(I4B)                                                  :: ios
INTEGER(I4B)                                                  :: nucleationpaths
INTEGER(I4B)                                                  :: kFlag
INTEGER(I4B)                                                  :: npFlag
INTEGER(I4B)                                                  :: jPoint

REAL(DP)                                                      :: StressMaxVal

INTEGER(I4B)                                                  :: nBoundaryConditionZone

LOGICAL(LGT)                                                  :: ConditionNameFound = .FALSE.

LOGICAL(LGT)                                                  :: NeedNucleationBlock

LOGICAL(LGT)                                                  :: lopen

LOGICAL(LGT)                                                  :: ExportGridLocations

namelist /Nucleation/                                          NameMineral,        &
                                                             label,              &
                                                             A_zero25C,          &
                                                             B_nucleation,       &
                                                             Sigma_mJm2,         &
                                                             SSA_m2g,            &
                                                             Surface

! ************************************
! Edit by Lucien Stolze, June 2023
CHARACTER (LEN=mls)                                           :: evapofile
CHARACTER (LEN=mls)                                           :: transpifile
INTEGER(I4B)                                                  :: tslength
REAL(DP), DIMENSION(:), ALLOCATABLE                           :: realmult_dum
CHARACTER (LEN=mls)                                           :: SnapshotFileFormat
LOGICAL(LGT)                                                  :: boolreg
integer :: IERR = 0
CHARACTER (LEN=mls)                                           :: watertablefile
CHARACTER (LEN=mls)                                           :: WatertableFileFormat
LOGICAL(LGT)                                                  :: readmineral
INTEGER(I4B)                                                  :: mineral_index
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: mineral_id
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: readmin_ssa
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: mineral_name
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: mineral_name_length
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: volfracfile
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: bsafile
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: lfile_volfrac
INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: lfile_bsa
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: FileFormatType_volfrac
CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE                :: FileFormatType_bsa
CHARACTER (LEN=mls)                                           :: vv_file 
INTEGER(I4B)                                                  :: vv_file_l 
CHARACTER (LEN=mls)                                           :: vv_fileformat
CHARACTER (LEN=mls)                                           :: bsa_file 
INTEGER(I4B)                                                  :: bsa_file_l 
CHARACTER (LEN=mls)                                           :: bsa_fileformat
INTEGER(I4B)                                                  :: min_id
CHARACTER (LEN=mls)                                           :: min_name
INTEGER(I4B)                                                  :: min_name_l
REAL(DP)                                                      :: t_default !default temp for zonation case
INTEGER(I4B)                                                  :: ssa_or_bsa
! ************************************
! Edit by Toshiyuki Bandai, 2023 May
INTEGER(I4B)                                 :: VG_error ! error flag for reading van Genuchten parameters
REAL(DP)                                     :: numerator ! used to compute permeability at faces
REAL(DP)                                     :: denominator ! used to compute permeability at faces
CHARACTER (LEN=mls)                          :: Richards_IC_File ! file name for Richards initial condition
CHARACTER (LEN=mls)                          :: Richards_IC_FileFormat ! file name for Richards initial condition (only single column is supported)
INTEGER(I4B)                                 :: BC_location ! ingeger to define the location of the boundary condition (0: lower boundary condition; 1: upper boundary condition)
CHARACTER (LEN=mls)                          :: upper_BC_file ! file name for the upper boundary condition for the Richards equation
CHARACTER (LEN=mls)                          :: lower_BC_file ! file name for the lower boundary condition for the Richards equation
CHARACTER (LEN=mls)                          :: infiltration_file ! file with infiltration data for "enviornmental_forcing" upper boundary condition
! End of edit by Toshiyuki Bandai, 2023 May
! ************************************

#if defined(ALQUIMIA)


include 'mpif.h'
integer :: rank, ierror
character(25) :: fn
#endif

ALLOCATE(realmult_dum(2000))
ALLOCATE(realmult(2000))

pi = DACOS(-1.0d0)

!fp! routine_name="start98";

vxfile = ' '
vyfile = ' '
vzfile = ' '
PorosityFile = ' '
SaturationFile = ' '
FixSaturation = 1.0d0
ModFlowCnv = 1.0
writeph = 'pH           '
ncomp = 0
nspec = 0
nrct = 0
ngas = 0
nchem = 0
ndecay = 0
nkin = 0
irestart = 0
restartfile = ' '
ihindmarsh = 1
gimrt = .TRUE.
petscon = .TRUE.
AppendRestart = .FALSE.
nCSD = 500

time_scale = 1.0d0
dist_scale = 1.0d0
OutputTimeScale = 1.0d0
OutputDistanceScale = 1.0d0
constantpor = 1.0d0
MinimumPorosity = 1.0D-14

nscratch = 9

iunit1 = 2
iunit2 = 3
jz = 1

time = 0.0

IF (NumInputFiles == 1) THEN
INQUIRE(FILE='PestControl.ant',EXIST=ext)
IF (EXT) THEN          !!  Pest Control file exists, so read input filename from it rather than prompting user
  OPEN(iunit1,FILE='PestControl.ant',STATUS='old',ERR=708)
  READ(iunit1,'(a)') filename
  CLOSE(iunit1,STATUS='keep')
  RunningPest = .TRUE.
ELSE                   !!  No Pestcontrol.ant file, so prompt user for the file name

  CALL get_command_argument(1,filename)

  IF (filename == '') THEN
    WRITE(*,*)
    WRITE(*,*) ' Type in your input file name'
    READ(*,'(a)') filename
  END IF

END IF
ELSE
filename = InputFile(InputFileCounter)
END IF

INQUIRE(FILE=filename,EXIST=ext)
IF (.NOT. ext) THEN
CALL stringlen(filename,ls)
WRITE(*,*)
WRITE(*,*) ' Cannot find input file: ', filename(1:ls)
WRITE(*,*)
READ(*,*)
STOP
END IF

OPEN(iunit1,FILE=filename,STATUS='old',ERR=703)
FileOutput = ' '
CALL stringlen(filename,lenInput)
IF (filename(lenInput-1:lenInput) == 'in' .OR. filename(lenInput-1:lenInput) == 'IN') THEN
FileOutput(1:lenInput-2) = filename(1:lenInput-2)
FileOutput(lenInput-1:lenInput+1) = 'out'
OPEN(iunit2,FILE=FileOutput,STATUS='unknown',ERR=704)
SaveInputFileName = FileOutput
ELSE
FileOutput = 'crunch.out'
END IF

CALL date_and_time(dumm1,dumm2,dumm3,curr_time)
str_mon = curr_time(2)
str_day = curr_time(3)
str_hr  = curr_time(5)
str_min = curr_time(6)
str_sec = curr_time(7)
str_millisec = curr_time(8)

nin = iunit1
nout = 4
#if !defined(ALQUIMIA)
inquire(file='CrunchJunk2.out',opened=lopen)
IF (lopen) THEN
CLOSE (unit=nout)
END IF
OPEN(UNIT=nout,FILE='CrunchJunk2.out',STATUS='unknown')
#else
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
write(fn,"(a10,i0,a4)")'CrunchJunk',rank,'.out'
write(*,*)fn
OPEN(UNIT=nout,FILE=fn,STATUS='unknown')
#endif

section = 'title'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
!!  WRITE(*,*)
!!  WRITE(*,*) ' Title block found'
!!  WRITE(*,*)
CALL read_title(nout,ltitle)
ELSE
WRITE(*,*)
WRITE(*,*) ' Failed to find title block'
WRITE(*,*)
END IF

WRITE(iunit2,*)
WRITE(iunit2,*) ' ************************** CrunchFlow ******************************'
WRITE(iunit2,*) '   '
WRITE(iunit2,*) '                  Authors:  C.I. STEEFEL, S. MOLINS '

WRITE(iunit2,*) '                      *** Copyright Notice ***          '
WRITE(iunit2,*) ' �CrunchFlow�, Copyright (c) 2016, The Regents of the University of California, &
                  through Lawrence Berkeley National Laboratory'
WRITE(iunit2,*) ' (subject to receipt of any required approvals from the U.S. Dept. of Energy).� All rights reserved.'
WRITE(iunit2,*)
WRITE(iunit2,*) ' If you have questions about your rights to use or distribute this software, please contact '
WRITE(iunit2,*) ' Berkeley Lab Innovation & Partnerships Office at��IPO@lbl.gov.  '
WRITE(iunit2,*)
WRITE(iunit2,*) ' NOTICE.� This Software was developed under funding from the U.S. Department of Energy and the U.S. Government '
WRITE(iunit2,*) ' consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting '
WRITE(iunit2,*) ' on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, &
                  distribute copies to the public, '
WRITE(iunit2,*) ' prepare derivative works, and perform publicly and display publicly, and to permit other to do so.'
WRITE(iunit2,*) '   '

WRITE(iunit2,*)
WRITE(iunit2,*)
WRITE(iunit2,1010) ltitle
WRITE(iunit2,*)
WRITE(*,*)
CALL stringlen(ltitle,lchar)
WRITE(*,*) ' Title of simulation: ', ltitle(1:72)
WRITE(*,*)

!  *****************RUNTIME BLOCK***********************
section = 'runtime'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
!!  WRITE(*,*)
!!  WRITE(*,*) ' Runtime parameters block found'
!!  WRITE(*,*)

! Now, find the individual parameters within block

!Addition of a walltime [min]] Lucien Stolze
call read_walltime(nin,nx,ny,nz)

parchar = 'gimrt'
parfind = ' '
gimrt = .TRUE.
CALL read_logical(nout,lchar,parchar,parfind,gimrt)
IF (gimrt) THEN
  WRITE(*,*)
  WRITE(*,*) ' --> Running in GIMRT mode (global implicit reaction and transport) '
  WRITE(*,*)
  os3d = .FALSE.
  petscon = .TRUE.
ELSE
  WRITE(*,*)
  WRITE(*,*) ' --> Running in OS3D mode (time splitting of reaction and transport) '
  WRITE(*,*)
  os3d = .TRUE.
  petscon = .FALSE.
END IF

nmmLogical = .FALSE.
parchar = 'nmm'
parfind = ' '
CALL read_logical(nout,lchar,parchar,parfind,nmmLogical)

CriticalZone = .FALSE.
parchar = 'criticalzone'
parfind = ' '
CALL read_logical(nout,lchar,parchar,parfind,CriticalZone)

parchar = 'saltcreep'
parfind = ' '
SaltCreep = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,SaltCreep)

parchar = 'calcitecreep'
parfind = ' '
CalciteCreep = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,CalciteCreep)

parchar = 'montterri'
parfind = ' '
MontTerri = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,MontTerri)

parchar = 'fracturenetwork'
parfind = ' '
FractureNetwork = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,FractureNetwork)

parchar = 'cubiclaw'
parfind = ' '
CubicLaw = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,cubiclaw)

parchar = 'courant_number'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter timestep_max not found
  courfactor = 0.5             ! Use default
ELSE
  courfactor = realjunk
END IF

IF (courfactor <= 0.0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Courant number should be greater than 0'
  WRITE(*,*) ' Courant number specified: ', courfactor
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

IF (os3d) THEN
  IF (courfactor > 1.0) then
    WRITE(*,*)
    WRITE(*,*) ' Using OS3D option, Courant number should be =< 1 '
  WRITE(*,*) ' Explicit transport is unstable at Courant numbers > 1'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

CALL units_time(nout,section,time_scale)

parchar = 'timestep_max'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter timestep_max not found
  tstep = 1.0             ! Use default
ELSE
  tstep = realjunk
  tstep = tstep*time_scale
END IF

parchar = 'timestep_init'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter timestep_max not found
  delt = 1.e-10           ! Use default
ELSE
  delt = realjunk
  ! Zhi Li commented out this line!
  delt = delt*time_scale
END IF
deltmin = delt

parchar = 'database'
parfind = ' '
data1 = ' '
CALL readCaseSensitive(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN  !
  data1 = ' '             ! Use default
ELSE
  data1 = dumstring
END IF

parchar = 'aqueousdatabase'
parfind = ' '
AqueousKineticFile = ' '
CALL readCaseSensitive(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN  !
  AqueousKineticFile = ' '
ELSE
  AqueousKineticFile = dumstring
END IF

parchar = 'catabolicdatabase'
parfind = ' '
CatabolicKineticFile = ' '
CALL readCaseSensitive(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN  !
  CatabolicKineticFile = ' '
ELSE
  CatabolicKineticFile = dumstring
END IF

IF (NumInputFiles > 1) THEN
  CONTINUE
ELSE
  ALLOCATE(stringarray(100))

  parchar = 'later_inputfiles'
  parfind = ' '
  lenarray = 0

  CALL read_multstring(nout,lchar,parchar,parfind,stringarray,lenarray,section)

  IF (parfind == ' ') THEN
    NumInputFiles = 1
    InputFileCounter = 1
  ELSE
    NumInputFiles = lenarray + 1

    IF (ALLOCATED(InputFile)) THEN
      DEALLOCATE(InputFile)
      ALLOCATE(InputFile(NumInputFiles))
    ELSE
      ALLOCATE(InputFile(NumInputFiles))
    END IF

    InputFile(1) = filename
    IF (NumInputFiles > 1) THEN
      DO i = 2,NumInputfiles
        InputFile(i) = stringarray(i-1)
      END DO
    END IF
    InputFileCounter = 1

  END IF

  DEALLOCATE(stringarray)

END IF

time_scale = 1.0d0

parchar = 'density_module'
parfind = ' '
DensityModule = ' '
CALL readCaseSensitive(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN  ! Parameter "time_units" not found
  DensityModule = 'temperature'             ! Use default
ELSE
  DensityModule = dumstring
END IF

! Check to see that the density module chosen is recognized

IF (DensityModule == 'temperature') THEN
  CONTINUE
ELSE IF (DensityModule == 'sodium_nitrate') THEN   !  Check later that sodium and nitrate are in system
  CONTINUE
ELSE IF (DensityModule == 'sodium_chloride') THEN  !  Check later that sodium and chloride are in system
  CONTINUE
ELSE IF (DensityModule == 'potassium_nitrate') THEN  !  Check later that potassium and chloride are in system
  CONTINUE
ELSE IF (DensityModule == 'calcium_nitrate') THEN  !  Check later that potassium and chloride are in system
  CONTINUE
ELSE IF (DensityModule == 'martin') THEN
CONTINUE
ELSE
  CALL stringlen(DensityModule,lchar)
  WRITE(*,*)
  WRITE(*,*) ' Density module not recognized: ', DensityModule(1:lchar)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

parchar = 'fix_saturation'
parfind = ' '
realjunk = 1.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter fix_saturation not found
  FixSaturation = 1.0d0
ELSE
  IF (realjunk > 0.0d0 .AND. realjunk <= 1.0d0) THEN
    FixSaturation = realjunk
    IF (realjunk == 1.0d0) THEN
      isaturate = 0
    ELSE
      isaturate = 1
      WRITE(*,*)
      WRITE(*,*) ' Liquid saturation: ', FixSaturation
      WRITE(*,*) ' Running as an unsaturated problem'
      WRITE(*,*)
    END IF
    WRITE(*,*) ' Constant liquid saturation: ', FixSaturation
    GO TO 5014
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Liquid saturation should be greater than zero and <= 1.0'
    WRITE(*,*) ' Liquid saturation value: ', realjunk
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

parchar = 'read_saturationfile'
parfind = ' '
SaturationFile = ' '
CALL readFileName(nout,lchar,parchar,parfind,dumstring,section,SaturationFileFormat)
IF (parfind == ' ') THEN
 SaturationFile = ' '             ! No default
!!  Check to make sure the user is not using the old designator "read_saturation"
  parchar = 'read_saturation'
  parfind = ' '
  CALL readFileName(nout,lchar,parchar,parfind,dumstring,section,SaturationFileFormat)
  IF (parfind == 'read_saturation') THEN
    WRITE(*,*)
    WRITE(*,*) 'Keyword "read_saturation" now obsolete--use "read_saturationfile"'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  SaturationFile = dumstring
END IF

5014 CONTINUE

parchar = 'CylindricalDivideVolume'
parfind = ' '
CylindricalDivideVolume = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,CylindricalDivideVolume)

parchar = 'Benchmark'
parfind = ' '
Benchmark = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,Benchmark)

parchar = 'DampRateInLowPorosity'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter DampRateInLowPorosity not found
  PorosityDamp = 1.0d0            ! Use default
ELSE
  PorosityDamp = realjunk
END IF

DampRateInLowPorosity = .FALSE.
IF (PorosityDamp < 1.0d0) THEN
  DampRateInLowPorosity = .TRUE.
END IF

!!    **********************  Solver Methods   ***********************

parchar = 'hindmarsh'
parfind = ' '
SolveHindmarsh = .FALSE.
ihindmarsh = 1
CALL read_logical(nout,lchar,parchar,parfind,SolveHindmarsh)
IF (SolveHindmarsh) THEN
  ihindmarsh = 1
  petscon = .FALSE.
ELSE
  ihindmarsh = 0
  IF (gimrt) THEN
    petscon = .TRUE.
  END IF
END IF

!!  Check first for generic solvers and preconditioners (not for GIMRT block solvers)

parchar = 'solver'
parfind = ' '
SolverMethod = ' '
CALL read_string(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN
 SolverMethod = 'bcgs'             ! Use default
ELSE
  SolverMethod = dumstring
END IF

! Check to see that the solver method is recognized, first for generic linear solves (diffusion, etc.)

IF (SolverMethod == 'gmres') THEN
  CONTINUE
ELSE IF (SolverMethod == 'bicg') THEN
  CONTINUE
ELSE IF (SolverMethod == 'bcgs') THEN
  CONTINUE
ELSE IF (SolverMethod == 'direct') THEN
  CONTINUE
ELSE
  CALL stringlen(SolverMethod,lchar)
  WRITE(*,*)
  WRITE(*,*) ' Solver method not recognized: ', SolverMethod(1:lchar)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

parchar = 'pc'
parfind = ' '
PCMethod = ' '
CALL read_string(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN
  PCMethod = 'ilu'             ! Use default
ELSE
  PCMethod = dumstring
END IF

IF (SolverMethod == 'direct') THEN
  PCMethod = 'direct'
END IF

! Check to see that the preconditioner method is recognized

IF (PCMethod == 'ilu') THEN
  CONTINUE
ELSE IF (PCMethod == 'jacobi') THEN
  CONTINUE
ELSE IF (PCMethod == 'direct') THEN
  CONTINUE
ELSE
  CALL stringlen(PCMethod,lchar)
  WRITE(*,*)
  WRITE(*,*) ' Generic preconditioner method not recognized: ', PCMethod(1:lchar)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

parchar = 'pclevel'
parfind = ' '
intjunk = 0
CALL read_integer(nout,lchar,parchar,parfind,intjunk,section)
IF (parfind == ' ') THEN
  Level = 5                 ! Use default
ELSE
  Level = intjunk
  IF (Level < 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Level for ILU preconditioner fill less than 0 not allowed'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
!!    IF (Level > 5) THEN
!!      Level = 5
!!    ELSE IF (Level < 0) THEN
!!      WRITE(*,*)
!!      WRITE(*,*) ' Level for ILU preconditioner fill less than 0 not allowed'
!!      WRITE(*,*)
!!      READ(*,*)
!!      STOP
!!    END IF
END IF
!!    *************************************************

parchar = 'gimrt_solver'
parfind = ' '
GIMRT_SolverMethod = ' '
CALL read_string(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN
 GIMRT_SolverMethod = 'gmres'             ! Use default
ELSE
  GIMRT_SolverMethod = dumstring
END IF

IF (GIMRT_SolverMethod == 'gmres') THEN
  CONTINUE
ELSE IF (GIMRT_SolverMethod == 'bcgs') THEN
  CONTINUE
ELSE
  CALL stringlen(GIMRT_SolverMethod,lchar)
  WRITE(*,*)
  WRITE(*,*) ' GIMRT solver method not recognized: ', GIMRT_SolverMethod(1:lchar)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

parchar = 'gimrt_pc'
parfind = ' '
GIMRT_PCMethod = ' '
CALL read_string(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN
  GIMRT_PCMethod = 'bjacobi'             ! Use default
ELSE
  GIMRT_PCMethod = dumstring
END IF

! Check to see that the preconditioner method is recognized

IF (GIMRT_PCMethod == 'bjacobi') THEN
  CONTINUE
ELSE IF (GIMRT_PCMethod == 'lu') THEN
  CONTINUE
ELSE IF (GIMRT_PCMethod == 'ilu') THEN
  CONTINUE
ELSE
  CALL stringlen(GIMRT_PCMethod,lchar)
  WRITE(*,*)
  WRITE(*,*) ' GIMRT preconditioner method not recognized: ', GIMRT_PCMethod(1:lchar)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

parchar = 'gimrt_pclevel'
parfind = ' '
intjunk = 0
CALL read_integer(nout,lchar,parchar,parfind,intjunk,section)
IF (parfind == ' ') THEN
  GimrtLevel = 1                 ! Use default
ELSE
  GimrtLevel = intjunk
  IF (GimrtLevel > 5) THEN
    GimrtLevel = 5
  ELSE IF (GimrtLevel < 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Level for ILU preconditioner fill less than 0 not allowed'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

parchar = 'gimrt_rtolksp'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter timestep_max not found
  GimrtRTOLKSP = 1.0D-09            ! Use default
ELSE
  GimrtRTOLKSP = realjunk
END IF
!!!  IF (GimrtRTOLKSP > 1.0D-07) THEN
!!!    GimrtRTOLKSP = 1.0D-07
!!!  END IF
IF (GimrtRTOLKSP < 1.0D-10) THEN
  GimrtRTOLKSP = 1.0D-10
END IF

!!   ************************************************

parchar = 'screen_output'
parfind = ' '
intjunk = 1
CALL read_integer(nout,lchar,parchar,parfind,intjunk,section)
IF (parfind == ' ') THEN
  IF (gimrt) THEN
    ScreenInterval = 1                 ! Use default
  ELSE
    ScreenInterval = 10
  END IF
ELSE
  ScreenInterval = intjunk
END IF


parchar = 'time_tolerance'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter timestep_max not found
  ttol = 0.001            ! Use default
ELSE
  ttol = realjunk
END IF

ResidualTolerance = 0.0d0
parchar = 'ResidualTolerance'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter timestep_max not found
  ResidualTolerance = 0.0d0            ! Use default
ELSE
  ResidualTolerance = realjunk
END IF

parchar = 'correction_max'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter timestep_max not found
  corrmax = 2.0           ! Use default
ELSE
  corrmax = realjunk
END IF

parchar = 'dissolution_max'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter timestep_max not found
  vdissmax = 0.001          ! Use default
ELSE
  vdissmax = realjunk
END IF

parchar = 'precipitation_max'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter timestep_max not found
  vpptmax = 0.001          ! Use default
ELSE
  vpptmax = realjunk
END IF

parchar = 'debye-huckel'
parfind = ' '
activity_dbh = .true.
CALL read_logical(nout,lchar,parchar,parfind,activity_dbh)
IF (activity_dbh) THEN
  igamma = 3
  WRITE(*,*) ' Extended Debye-Huckel activity model'
ELSE
  igamma = 0
  WRITE(*,*) ' Unit activity coefficients'
END IF

parchar = 'lag_activity'
parfind = ' '
lag_activity =  .true.
CALL read_logical(nout,lchar,parchar,parfind,lag_activity)
IF (activity_dbh) THEN
  IF (lag_activity) THEN
    igamma = 3
!!      WRITE(*,*) ' Lagging activity coefficients by one timestep'
  ELSE
    igamma = 2
    WRITE(*,*) ' Updating activity coeffs every Newton step'
  END IF
ELSE
  igamma = 0
END IF

parchar = 'database_sweep'
parfind = ' '
database_sweep = .false.
CALL read_logical(nout,lchar,parchar,parfind,database_sweep)
IF (database_sweep) THEN
  icomplete = 1
  WRITE(*,*) ' Sweeping database to find additional species, gases, and minerals'
ELSE
  icomplete = 0
END IF

parchar = 'ReadGeochemicalConditions'
parfind = ' '
ReadGeochemicalConditions = .false.
CALL read_logical(nout,lchar,parchar,parfind,ReadGeochemicalConditions)

parchar = 'ReadGautier'
parfind = ' '
ReadGautier = .false.
CALL read_logical(nout,lchar,parchar,parfind,ReadGautier)

parchar = 'Qingyun'
parfind = ' '
Qingyun = .false.
CALL read_logical(nout,lchar,parchar,parfind,Qingyun)

parchar = 'ForsteriteCapillary'
parfind = ' '
ForsteriteCapillary = .false.
CALL read_logical(nout,lchar,parchar,parfind,ForsteriteCapillary)

parchar = 'HanfordStrontium'
parfind = ' '
HanfordStrontium = .false.
CALL read_logical(nout,lchar,parchar,parfind,HanfordStrontium)

parchar = 'GMSsecondary'
parfind = ' '
GMSsecondary = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,GMSsecondary)

parchar = 'GMSmineral'
parfind = ' '
GMSmineral = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,GMSmineral)

IF (GMSsecondary) THEN
  IF (GMSmineral) THEN
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*) ' Both GMSsecondary and GMSmineral options selected'
    WRITE(*,*) ' ---> If you want to output a sweep of the mineral database, '
    WRITE(*,*) '          Turn GMSsecondary off   '
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

IF (GMSsecondary .OR. GMSmineral) THEN
  icomplete = 1
END IF

parchar = 'reaction_path'
parfind = ' '
reaction_path = .false.
CALL read_logical(nout,lchar,parchar,parfind,reaction_path)
IF (reaction_path) THEN
  ipath = 1
ELSE
  ipath = 0
END IF

!  Logical "speciate_only" instructs code to skip mineral kinetics section

parchar = 'speciate_only'
parfind = ' '
onlyspeciate = .false.
CALL read_logical(nout,lchar,parchar,parfind,onlyspeciate)
IF (onlyspeciate) THEN
  ispeciate = 1
ELSE
  ispeciate= 0
END IF

!  Logical "generic_rates" turns on generic rates (database read skipped)

!!  parchar = 'generic_rates'
!!  parfind = ' '
!!  genericrates = .false.
!!  CALL read_logical(nout,lchar,parchar,parfind,genericrates)
!!  IF (genericrates) THEN
!!    igenericrates = 1
!!  ELSE
!!    igenericrates = 0
!!  END IF

parchar = 'generic_rates'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter "generic_rates" not found
  genericrates = .FALSE.
  igenericrates = 0
ELSE
  genericrates =  .TRUE.
  igenericrates = 1
  RateGeneric = realjunk
END IF

!  If speciation only is specified, turn off "generic rates" option

IF (ispeciate == 1) THEN
  igenericrates = 0
END IF

restartfile = ' '
CALL read_restart(nout)
IF (restartfile == ' ') THEN
  irestart = 0
ELSE
  irestart = 1
END IF

parchar = 'save_restart'
parfind = ' '
RestartOutputFile = ' '
CALL readCaseSensitive(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN  !
  RestartOutputFile = 'crunch.rst'             ! Use default
ELSE
  RestartOutputFile = dumstring
END IF

Rectangular = .TRUE.
Cylindrical = .FALSE.
CALL read_coordinates(nout)

NuftFile = ' '
CALL read_nuft(nout)
IF (NuftFile == ' ') THEN
  ReadNuft = .FALSE.
ELSE
  ReadNuft = .TRUE.
END IF

xtool = .FALSE.
tecplot = .TRUE.
xmgr = .FALSE.
kaleidagraph = .FALSE.
nview = .FALSE.
originlab = .FALSE.
CALL read_graphics(nout)

IF (tecplot .OR. kaleidagraph .OR. originlab) THEN
  CONTINUE
ELSE
  kaleidagraph = .TRUE.
END IF

master = ' '
CALL read_master(nout)

parchar = 'streamtube'
parfind = ' '
streamtube = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,streamtube)
parchar = 'steady_state'
parfind = ' '
RunToSteady = .false.
CALL read_steady(nout,lchar,parchar,parfind,RunToSteady)

parchar = 'giambalvo'
parfind = ' '
giambalvo = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,giambalvo)

parchar = 'KateMaher'
parfind = ' '
KateMaher = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,KateMaher)

parchar = 'JennyDruhan'
parfind = ' '
JennyDruhan = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,JennyDruhan)

parchar = 'JennyFirstOrder'
parfind = ' '
JennyFirstOrder = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,JennyFirstOrder)

parchar = 'Duan'
parfind = ' '
Duan = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,Duan)

parchar = 'Duan2006'
parfind = ' '
Duan2006 = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,Duan2006)
IF (Duan) THEN
  Duan2006 = .FALSE.
END IF

parchar = 'Maggi'
parfind = ' '
Maggi = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,Maggi)

parchar = 'DePaolo'
parfind = ' '
DePaolo = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,DePaolo)


parchar = 'SetSurfaceAreaConstant'
parfind = ' '
SetSurfaceAreaConstant = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,SetSurfaceAreaConstant)

!! Inhibit the consumption of minerals for running model spinup: (Stolze Lucien)
parchar = 'model_spinup'
parfind = ' '
spinup = .false.
CALL read_logical(nout,lchar,parchar,parfind,spinup)

!! Cyclical reading of all time series based on 1 year time series (applied to infiltration, evaporation, transpiration, temperature) Stolze Lucien
TS_1year = .FALSE.
parchar = 'timeseries_cyclic_1year'
parfind = ' '
CALL read_logical(nout,lchar,parchar,parfind,TS_1year)

!! Keep biomass fixed, Stolze Lucien
parchar = 'biomassfixed'
parfind = ' '
biomassfixed = .false.
CALL read_logical(nout,lchar,parchar,parfind,biomassfixed)

!! Generate velocity vector for velocity_read
parchar = 'generate_velocity_vector'
parfind = ' '
generate_velocity_vector = .false.
CALL read_logical(nout,lchar,parchar,parfind,generate_velocity_vector)

parchar = 'Inagaki'
parfind = ' '
inagaki = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,inagaki)

parchar = 'InagakiDensify'
parfind = ' '
inagaki = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,inagaki2)

parchar = 'OelkersRateLaw'
parfind = ' '
OelkersRateLaw = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,OelkersRateLaw)

parchar = 'BurchRateLaw'
parfind = ' '
BurchRateLaw = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,BurchRateLaw)

parchar = 'HellmannRateLaw'
parfind = ' '
HellmannRateLaw = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,HellmannRateLaw)

parchar = 'SilicaRateLaw'
parfind = ' '
SilicaRateLaw = .FALSE.
CALL read_logical(nout,lchar,parchar,parfind,SilicaRateLaw)

IF (JennyDruhan) THEN

  parchar = 'UseBulkMineral'
  parfind = ' '
  UseBulkMineral = .TRUE.
  CALL read_logical(nout,lchar,parchar,parfind,UseBulkMineral)

END IF


IF (KateMaher) THEN
  parchar = 'BurchRateLaw'
  parfind = ' '
  BurchRateLaw = .FALSE.
  CALL read_logical(nout,lchar,parchar,parfind,BurchRateLaw)

  parchar = 'OelkersRateLaw'
  parfind = ' '
  OelkersRateLaw = .FALSE.
  CALL read_logical(nout,lchar,parchar,parfind,OelkersRateLaw)

  parchar = 'DistributionCalcite'
  parfind = ' '
  realjunk = 0.0
  CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
  IF (parfind == ' ') THEN  ! Parameter minimum_porosity not found
    WRITE(*,*)
    WRITE(*,*) ' When running with "KateMaher" option, a value for "DistributionCalcite" must be specified in RUNTIME block'
    WRITE(*,*)
    STOP
  ELSE
    IF (realjunk >= 0.0) THEN
      DistributionCalcite = realjunk
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Distribution coefficient for calcite should be > or = 0'
      WRITE(*,*) ' Distribution coefficient for calcite = ', DistributionCalcite
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  END IF

END IF

parchar = 'OvershootTolerance'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter "OvershootTolerance" not found
  voltol = 1.0D-05             ! Use default
ELSE
  voltol = realjunk
END IF

LagSurface = 0.0

IF (JennyDruhan) THEN
  parchar = 'LagSurface'
  parfind = ' '
  realjunk = 0.0
  CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
  IF (parfind == ' ') THEN  ! Parameter "LagSurface" not found
    LagSurface = 0.00             ! Use default
  ELSE
    LagSurface = realjunk
  END IF
END IF

ELSE
WRITE(*,*) ' Failed to find runtime parameters block'
WRITE(*,*) ' ---> Using default values'
WRITE(*,*)
tstep =  1.0
delt = 1.e-09
ttol = 0.005
corrmax = 2.0
vdissmax = 0.001
vpptmax  = 0.001
master  = ' '  ! If none provided, get one later from species list
igamma  = 3
icomplete = 0
ipath = 1
ispeciate= 0
igenericrates = 0
irestart = 0
gimrt = .TRUE.
os3d = .FALSE.
modflow = .FALSE.
RunToSteady = .FALSE.
giambalvo = .FALSE.
xtool = .FALSE.
tecplot = .FALSE.
xmgr = .FALSE.
nview = .FALSE.
kaleidagraph = .TRUE.
ReadNuft = .FALSE.
Rectangular = .TRUE.
Cylindrical = .FALSE.
Spherical = .FALSE.
courfactor = 0.5
DensityModule = 'temperature'
SolverMethod = 'bicg'
level = 2
GIMRT_PCMethod = 'bjacobi'
GIMRT_SolverMethod = 'gmres'
PCMethod = 'ilu'
GMSsecondary = .FALSE.
GMSmineral = .FALSE.
data1 = ' '
PorosityFile = ' '
SaturationFile = ' '
FixSaturation = 1.0
ihindmarsh = 1
ScreenInterval = 1
RestartOutputFile = 'crunch.rst'
NumInputFiles = 1
InputFileCounter = 1
KateMaher = .FALSE.
BurchRateLaw = .FALSE.
OelkersRateLaw = .FALSE.
HellmannRateLaw = .FALSE.
SilicaRateLaw = .FALSE.
voltol = 1.0D-05
END IF

!  ***********  Database block  *****************

section = 'database'
CALL readblock(nin,nout,section,found,ncount)

IF (data1 == ' ') THEN
IF (found) THEN
  CALL read_dbs(nout,data1)
ELSE
  CONTINUE
END IF
ELSE
CONTINUE
END IF

WRITE(iunit2,*)
WRITE(iunit2,*)
WRITE(iunit2,2101) data1
WRITE(iunit2,*)
WRITE(*,*)
WRITE(*,2101) data1

!  ***********  End Database block  ****************

!  *****  TEMPERATURE BLOCK  ********************

section = 'temperature'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
!!  WRITE(*,*) ' Temperature parameters block found'
!!  WRITE(*,*)

!**************
!Temperature fixed and homogeneous
parchar = 'set_temperature'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter set_temperature not found
  tinit = 25.0
  tgrad = 0.0
  jtemp = 0
ELSE
  tinit = realjunk
  tgrad = 0.0
  jtemp = 0
!!   WRITE(*,5012)  tinit
END IF
!**************

!**************
!Temperature gradient
parchar = 'temperature_gradient'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter temperature_grad not found
  tgrad = 0.0
  jtemp = 0
!!    WRITE(*,*) ' No temperature gradient specified'
ELSE
  tgrad = realjunk
  jtemp = 1
  WRITE(*,5013)  tgrad
END IF
!**************

!**************
!Read temperature distribution from file
parchar = 'read_temperaturefile'
parfind = ' '
TFile = ' '
CALL readFileName(nout,lchar,parchar,parfind,dumstring,section,TemperatureFileFormat)
IF (parfind == ' ') THEN
 TFile = ' '             ! No default
!!  Check to make sure the user is not using the old designator "read_temperature"
  parchar = 'read_temperature'
  parfind = ' '
  CALL readFileName(nout,lchar,parchar,parfind,dumstring,section,TemperatureFileFormat)
  IF (parfind == 'read_temperature') THEN
    WRITE(*,*)
    WRITE(*,*) 'Keyword "read_temperature" now obsolete--use "read_temperaturefile"'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  TFile = dumstring
  jtemp = 2
END IF
!**************

  parchar = 'RunIsothermal'
  parfind = ' '
  RunIsothermal = .false.
  CALL read_logical(nout,lchar,parchar,parfind,RunIsothermal)

  IF (RunIsothermal .AND. tgrad /= 0.0d0) THEN
      WRITE(*,*)
      WRITE(*,*) ' Isothermal run incompatible with temperature gradient'
      WRITE(*,*)
      READ(*,*)
      STOP
  END IF
  IF (RunIsothermal .AND. tfile /= ' ') THEN
      WRITE(*,*)
      WRITE(*,*) ' Isothermal run incompatible with temperature file read'
      WRITE(*,*)
      READ(*,*)
      STOP
  END IF


  !****************
  !Edit by Lucien Stolze, June 2023
  !Temperature time series + zonation
  ! ************************************
  nb_temp_ts = 0
  nb_temp_fix = 0
  boolreg = .false.
  RunTempts = .false.  
  ! ************************************
  CALL read_tempreg(nout,dumstring,TemperatureFileFormat,boolreg)
  IF (boolreg) THEN
    jtemp = 3
    TFile = dumstring
  ENDIF
  t_default = 25
  CALL read_tempts(nout,tslength,t_default)
  !****************


ELSE

WRITE(*,*) ' Temperature parameters not found'
WRITE(*,*) ' Using defaults'
jtemp = 0
tinit = 25.0
tgrad = 0.0
TFile = ' '

END IF

!*************************************************************

!  ****************CALL READ98********************************

CALL FirstAllocation()

dxxt = 0.0d0
dyyt = 0.0d0
dzzt = 0.0d0
dxxt = 0.0d0
dyyt = 0.0d0
dzzt = 0.0d0
muaq = 0.0d0
mugas = 0.0d0
musurf = 0.0d0
muexc = 0.0d0
ndepend = 0
kcrossaff = 0
nmonod = 0
ninhibit = 0
ndepend = 0
bfit = 0.0
jxseries = 1
jyseries = 1
jzseries = 1
iedl = 0

if (.not.allocated(pH)) ALLOCATE(pH(mchem))
if (.not.allocated(guesspH)) ALLOCATE(guesspH(mchem))
if (.not.allocated(constraint)) ALLOCATE(constraint(nc,mchem))

IF (Duan .OR. Duan2006) THEN
if (.not.allocated(vrINitial)) ALLOCATE(vrINitial(mchem))
END IF


AffinityDepend1 = 1.0d0
AffinityDepend2 = 1.0d0
AffinityDepend3 = 1.0d0

! biomass
!!chi = 1
LocalEquilibrium = .FALSE.
!!BQ = 0.0d0
! biomass end

CALL read98(ncomp,nspec,nkin,nrct,ngas,nsurf,nsurf_sec,data1,icomplete,  &
  ispeciate,igenericrates,GMSsecondary,GMSmineral,RateGeneric)

!!  Check to see if a NUCLEATION block is required based on what is in database file (type = nucleation, imintype = 10)

NeedNucleationBlock = .FALSE.
DO k = 1,nrct
DO np = 1,nreactmin(k)
  IF (imintype(np,k) == 10) THEN
     npFlag = np
     kFlag = k
     NeedNucleationBlock = .TRUE.
  END IF
END DO
END DO


!!  Now check for nucleation block

!      *****************NUCLEATION SECTION***********************

section = 'nucleation'

CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN

IF (ALLOCATED(NucleationMineral)) THEN
DEALLOCATE(NucleationMineral)
END IF
ALLOCATE(NucleationMineral(50))
NucleationMineral = ' '

IF (ALLOCATED(NameNucleationPathway)) THEN
DEALLOCATE(NameNucleationPathway)
END IF
ALLOCATE(NameNucleationPathway(50))
NameNucleationPathway = ' '

IF (ALLOCATED(NucleationSurface)) THEN
DEALLOCATE(NucleationSurface)
END IF
ALLOCATE(NucleationSurface(nreactmax,nrct))
NucleationSurface = 0

IF (ALLOCATED(kNucleationPath)) THEN
DEALLOCATE(kNucleationPath)
END IF
ALLOCATE(kNucleationPath(50))
kNucleationPath = 0
IF (ALLOCATED(npNucleationPath)) THEN
DEALLOCATE(npNucleationPath)
END IF
ALLOCATE(npNucleationPath(50))
npNucleationPath = 0

IF (ALLOCATED(Azero25C)) THEN
  DEALLOCATE(Azero25C)
  ALLOCATE(Azero25C(nreactmax,nrct))
ELSE
  ALLOCATE(Azero25C(nreactmax,nrct))
END IF
Azero25C = 0.0d0
IF (ALLOCATED(Bnucleation)) THEN
  DEALLOCATE(Bnucleation)
  ALLOCATE(Bnucleation(nreactmax,nrct))
ELSE
  ALLOCATE(Bnucleation(nreactmax,nrct))
END IF
Bnucleation = 0.0d0
IF (ALLOCATED(sigmaNucleation)) THEN
  DEALLOCATE(sigmaNucleation)
  ALLOCATE(sigmaNucleation(nreactmax,nrct))
ELSE
  ALLOCATE(sigmaNucleation(nreactmax,nrct))
END IF
sigmaNucleation = 0.0d0
IF (ALLOCATED(SurfaceAreaNucleation)) THEN
  DEALLOCATE(SurfaceAreaNucleation)
  ALLOCATE(SurfaceAreaNucleation(nreactmax,nrct))
ELSE
  ALLOCATE(SurfaceAreaNucleation(nreactmax,nrct))
END IF
SurfaceAreaNucleation = 0.0d0
IF (ALLOCATED(SumMineralSurfaceArea)) THEN
  DEALLOCATE(SumMineralSurfaceArea)
  ALLOCATE(SumMineralSurfaceArea(nreactmax,nrct))
ELSE
  ALLOCATE(SumMineralSurfaceArea(nreactmax,nrct))
END IF
IF (ALLOCATED(HomogeneousNucleation)) THEN
  DEALLOCATE(HomogeneousNucleation)
  ALLOCATE(HomogeneousNucleation(nreactmax,nrct))
ELSE
  ALLOCATE(HomogeneousNucleation(nreactmax,nrct))
END IF
SumMineralSurfaceArea = .FALSE.
HomogeneousNucleation = .FALSE.

!!  Based on MINERAL block in input file, identify nucleation pathways

knucl = 0
DO k = 1,nrct
  DO np = 1,nreactmin(k)
    IF (imintype(np,k) == 10) THEN
      knucl = knucl + 1
      NucleationMineral(knucl) = umin(k)
      NameNucleationPathway(knucl) =   rlabel(np,k)
      kNucleationPath(knucl) = k
      npNucleationPath(knucl) = np
    END IF
  END DO
END DO

NucleationPaths = knucl

REWIND nout

!!! Loop over nucleation pathways provided in input file (or database file)
do_input_pathways: DO knucl=1,nucleationpaths

! read all necessary pathways (reactions) from file
  do_nucleationpathways: DO     !!  This is a loop through the multiple? namelist entries in the NUCLEATION block

!!!   initialize namelist variables before reading it in from file
!!!  Mineral        = Calcite
!!!  label          = nucleatecalcite
!!!  Azero25C       = 0.01
!!!  Bnucleation    = 0.009
!!!  Sigma(mJ/m2)   = 97.0
!!!  SSA(m2/g)      = 1.0
!!!  Surface        = all


!     read 'Nucleation' namelist from file
    read(nout,nml=Nucleation,iostat=ios)

!     result from read (IOS)
    IF (ios == 0) THEN

!     successful read, compare to nucleation pathway specified in input file (or database file??)

!! Point from list of names derived from input file (or database) to what is read in the namelist
      IF (NucleationMineral(knucl) == NameMineral .and. NameNucleationPathway(knucl)== label) THEN

        k = kNucleationPath(knucl)
        np = npNucleationPath(knucl)
        Azero25C(np,k) = A_zero25C
        Bnucleation(np,k) = B_nucleation
        SigmaNucleation(np,k) = Sigma_mJm2
        SurfaceAreaNucleation(np,k) = SSA_m2g

        IF (surface == 'none' .OR. surface == 'NONE' .OR. surface == 'None.') THEN
          HomogeneousNucleation(np,k) = .TRUE.
        ELSE IF (surface == 'all' .OR. surface == 'ALL' .OR. surface == 'All') THEN
          HomogeneousNucleation(np,k) = .FALSE.
          SumMineralSurfaceArea(np,k) = .TRUE.
        ELSE
          HomogeneousNucleation(np,k) = .FALSE.
          SumMineralSurfaceArea(np,k) = .FALSE.
          !! Find the mineral number for the nucleation surface
          NucleationSurface(np,k) = 0
          DO kk = 1,nrct
            IF (umin(kk) == Surface) then
              NucleationSurface(np,k) = kk
            END IF
          END DO

          IF (NucleationSurface(np,k) == 0) THEN
            WRITE(*,*)
            WRITE(*,*) ' Mineral surface for nucleation not found in list'
            WRITE(*,*) ' Looking for: ', Surface
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF

        END IF


        rewind(nout)

        EXIT do_nucleationpathways

      ELSE

        CYCLE do_nucleationpathways

      END IF

  else if (ios < 0) then

!     no more nucleation pathways to read
    write(*,*)'End of file'
    exit do_nucleationpathways

  else if (ios > 0) then

!!      write(*,nml=Nucleation)
    write(*,*)' Error reading Nucleation namelist: goodbye'
    stop

  end if

end do do_nucleationpathways

end do do_input_pathways

ELSE

IF (NeedNucleationBlock) THEN
  WRITE(*,*)
  WRITE(*,*) ' Nucleation type rate law found listed in database'
  WRITE(*,*) ' No NUCLEATION block found'
  WRITE(*,*)
  write(*,*) ' Mineral: ',umin(kFlag)
  WRITE(*,*) ' Label:   ',rlabel(npFlag,kFlag)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

END IF          !!!!  End of nucleation read

IF (master == ' ') then
DO ik = 1,ncomp+nspec
  IF (ulab(ik) == 'H+' .OR. ulab(ik) == 'pH' .OR. ulab(ik) == 'ph') THEN
    master = 'H+'
  END IF
END DO
DO ik = 1,ncomp+nspec
  IF (ulab(ik) == 'O2(aq)') THEN
    master = 'O2(aq)'
  END IF
END DO
IF (master == ' ') THEN
  master = ulab(1)
END IF
END IF

!! Check to see if CO2(aq) is a primary species if Duan activity option is selected

CheckDuan = .FALSE.
DO i = 1,ncomp
IF (ulab(i) == 'CO2(aq)') THEN
  CheckDuan = .TRUE.
END IF
END DO

IF (.NOT. CheckDuan) THEN
IF (Duan) THEN
  write(*,*)
  write(*,*) ' Duan option should be used with CO2(aq) as a primary species'
  write(*,*)
  read(*,*)
  stop
ELSE IF (Duan2006) THEN
  write(*,*)
  write(*,*) ' Duan2006 option should be used with CO2(aq) as a primary species'
  write(*,*)
  read(*,*)
  stop
END IF
END IF


IF (OelkersRateLaw) THEN
WRITE(*,*)
WRITE(*,*) ' Oelkers Rate Law no longer hardwired--Use HyperbolicInhibition in "Dependence" '
WRITE(*,*)
READ(*,*)
STOP

!!  kUPlag = 1
!!  ikAl = 0
!!  DO ik = 1,ncomp
!!    IF (ulab(ik) == 'HAlO2(aq)') THEN
!!      ikAl = ik
!!    END IF
!!  END DO

!!  IF (ikAl == 0) THEN
!!    WRITE(*,*)
!!    WRITE(*,*) ' Primary species "HAlO2(aq)" not found when using "OelkersRateLaw" option '
!!    WRITE(*,*)
!!    STOP
!!  END IF

END IF

IF (BurchRateLaw .OR. HellmannRateLaw .OR. SilicaRateLaw) THEN
kUPlag = 1
END IF

IF (KateMaher) THEN

kUCalcite = 0
DO k = 1,nkin
  IF (umin(k) == 'Uranium-Calcite') THEN
    kUCalcite = k
  END IF
END DO
IF (kUCalcite == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Mineral number of "Uranium-Calcite" when using "KateMaher" option not found '
  WRITE(*,*) ' Mineral should be specified as "Uranium-Calcite" in database and input file'
  WRITE(*,*)
  STOP
END IF

kMarineCalcite = 0
DO k = 1,nkin
  IF (umin(k) == 'Marine-Calcite') THEN
    kMarineCalcite = k
  END IF
END DO
IF (kMarineCalcite == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Mineral number of "Marine-Calcite" when using "KateMaher" option not found '
  WRITE(*,*) ' Mineral should be specified as "Marine-Calcite" in database and input file'
  WRITE(*,*)
  STOP
END IF

kUPlag = 0
DO k = 1,nkin
  IF (umin(k) == 'Plag_U') THEN
    kUPlag = k
  END IF
END DO
IF (kUPlag == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Mineral number of "Plag_U" when using "KateMaher" option not found '
  WRITE(*,*) ' Mineral should be specified as "Plag_U" in database and input file'
  WRITE(*,*)
  STOP
END IF

ikCa = 0
ik234U = 0
ik238U = 0
ikCO3 = 0
ikAl = 0
DO ik = 1,ncomp+nspec
  IF (ulab(ik) == 'Ca++') THEN
    ikCa = ik
  END IF
  IF (ulab(ik) == 'CO3--') THEN
    ikCO3 = ik
  END IF
  IF (ulab(ik) == 'U_234O2++') THEN
    ik234U = ik
  END IF
  IF (ulab(ik) == 'U_238O2++') THEN
    ik238U = ik
  END IF
END DO

IF (OelkersRateLaw) THEN
  DO ik = 1,ncomp
    IF (ulab(ik) == 'HAlO2(aq)') THEN
      ikAl = ik
    END IF
  END DO

  IF (ikAl == 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Primary species "HAlO2(aq)" not found when using "OelkersRateLaw" option '
    WRITE(*,*)
    STOP
  END IF
END IF

IF (ikCa == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Species number for "Ca++" not found when using "KateMaher" option '
  WRITE(*,*)
  STOP
END IF
IF (ikCO3 == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Species number for "CO3--" not found when using "KateMaher" option '
  WRITE(*,*)
  STOP
END IF
IF (ik234U == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Species number for "U_234O2++" not found when using "KateMaher" option '
  WRITE(*,*)
  STOP
END IF
IF (ik238U == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Species number for "U_238O2++" not found when using "KateMaher" option '
  WRITE(*,*)
  STOP
END IF
END IF

!!  The following is now disabled, since the AffinityDependence are read in the database

IF (KateMaher .AND. BurchRateLaw) THEN
WRITE(*,*)
WRITE(*,*) ' Keyword "BurchRateLaw" is disabled--set AffinityDepend in kinetic database'
WRITE(*,*)
STOP
!!  AffinityDepend1 = 1.00d0
!!  AffinityDepend2 = 1.00d0
!!  AffinityDepend3 = 1.00d0
!!  AffinityDepend1(1,kUPlag) = 100000.00d0
ELSE IF (OelkersRateLaw) THEN
!!  AffinityDepend1 = 1.00d0
!!  AffinityDepend2 = 1.00d0
!!  AffinityDepend3 = 1.00d0
!!  AffinityDepend2(1,kUPlag) = 0.3333333333333d0
ELSE IF (BurchRateLaw) THEN
WRITE(*,*)
WRITE(*,*) ' Keyword "BurchRateLaw" is disabled--set AffinityDepend in kinetic database'
WRITE(*,*)
STOP
!!  AffinityDepend1 = 1.00d0
!!  AffinityDepend2 = 1.00d0
!!  AffinityDepend3 = 1.00d0
!!  AffinityDepend1(1,kUPlag) = 100000.00d0
!!  AffinityDepend1 = 1.00d0
!!  AffinityDepend2 = 1.00d0
!!  AffinityDepend3 = 1.00d0
!!  AffinityDepend1(2,kUPlag) = 1.45
!!  AffinityDepend2(1,kUPlag) = 8.4D-17
!!  AffinityDepend3(1,kUPlag) = 15.0
ELSE IF (HellmannRateLaw) THEN
WRITE(*,*)
WRITE(*,*) ' Keyword "HellmannRateLaw" is disabled--set AffinityDepend in kinetic database'
WRITE(*,*)
STOP
!!  AffinityDepend1 = 1.00d0
!!  AffinityDepend2 = 1.00d0
!!  AffinityDepend3 = 1.00d0
!!  AffinityDepend1(2,kUPlag) = 1.17d0
!!  AffinityDepend2(1,kUPlag) = 0.0000798
!!  AffinityDepend3(1,kUPlag) = 3.81d0
ELSE IF (SilicaRateLaw) THEN
WRITE(*,*)
WRITE(*,*) ' Keyword "SilicaRateLaw" is disabled--set AffinityDepend in kinetic database'
WRITE(*,*)
STOP
!!  AffinityDepend1 = 1.00d0
!!  AffinityDepend2 = 1.00d0
!!  AffinityDepend3 = 1.00d0
!!  AffinityDepend2(1,kUPlag) = 2.000d0
!!ELSE
!!  AffinityDepend1 = 1.00d0
!!  AffinityDepend2 = 1.00d0
!!  AffinityDepend3 = 1.00d0
ELSE
CONTINUE
END IF

!  Output info about minerals here


WRITE(iunit2,*) ' ***KINETIC INPUTS***'
WRITE(iunit2,*)
WRITE(iunit2,*) '  MINERAL'
WRITE(iunit2,*)

!      do k = 1,nkin
!        write(iunit2,1101) umin(k)
!        write(iunit2,1102) nreact(k)
!        write(iunit2,1104) thresh(np,k)
!        do ll = 1,nreact(k)
!           write(iunit2,1105) ll,npre(ll,k)
!           write(iunit2,1106) rate0(ll,k)/secyr
!           write(iunit2,1108) Ea(ll,k)
!           write(iunit2,1103) sat1(ll,k),sat2(ll,k)
!           if (npre(ll,k).gt.0) then
!            write(iunit2,*) '     Species              Exp. Dependence'
!             do mm = 1,npre(ll,k)
!               write(iunit2,1107)  ulab(ispec(mm,ll,k)),depend(mm,ll,k)
!             end do
!           endif
!        end do
!        write(iunit2,*)
!      end do


DO i = 1,nc
DO nco = 1,mchem
  constraint(i,nco) = ' '
END DO
END DO

!         *********** ION EXCHANGE SECTION **************

section='ion_exchange'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
!!  WRITE(*,*)
!!  WRITE(*,*) ' Ion exchange block found'
!!  WRITE(*,*)

CALL read_exchange(nout,ncomp,nexchange,data1,nexch_sec,nkin)

IF (ALLOCATED(icec)) THEN
  DEALLOCATE(icec)
  ALLOCATE(icec(nexchange))
ELSE
  ALLOCATE(icec(nexchange))
END IF

!  If there is exchange, change to IGAMMA = 2 (no lag of activity coefficients)

!  IF (igamma == 3 .AND. nexchange > 0) THEN
!    igamma = 2
!  END IF

ELSE
WRITE(*,*) ' No ion exchange block found'
END IF

!  Now, check to see that species dependences specified for mineral
!  rate laws involve secondary aqueous species, exchange species, or
!  surface complexes (only primary species checked so far)

IF (ALLOCATED(ndependsurf)) THEN
DEALLOCATE(ndependsurf)
ALLOCATE(ndependsurf(nreactmax,nkin))
ELSE
ALLOCATE(ndependsurf(nreactmax,nkin))
END IF
IF (ALLOCATED(ndependex)) THEN
DEALLOCATE(ndependex)
ALLOCATE(ndependex(nreactmax,nkin))
ELSE
ALLOCATE(ndependex(nreactmax,nkin))
END IF
IF (ALLOCATED(ixdepend)) THEN
DEALLOCATE(ixdepend)
ALLOCATE(ixdepend(nexch_sec,nreactmax,nkin))
ELSE
ALLOCATE(ixdepend(nexch_sec,nreactmax,nkin))
END IF
IF (ALLOCATED(isdepend)) THEN
DEALLOCATE(isdepend)
ALLOCATE(isdepend(nsurf+nsurf_sec,nreactmax,nkin))
ELSE
ALLOCATE(isdepend(nsurf+nsurf_sec,nreactmax,nkin))
END IF
IF (ALLOCATED(dependex)) THEN
DEALLOCATE(dependex)
ALLOCATE(dependex(nexch_sec,nreactmax,nkin))
ELSE
ALLOCATE(dependex(nexch_sec,nreactmax,nkin))
END IF
IF (ALLOCATED(dependsurf)) THEN
DEALLOCATE(dependsurf)
ALLOCATE(dependsurf(nsurf+nsurf_sec,nreactmax,nkin))
ELSE
ALLOCATE(dependsurf(nsurf+nsurf_sec,nreactmax,nkin))
END IF
IF (ALLOCATED(ispot)) THEN
DEALLOCATE(ispot)
ALLOCATE(ispot(nsurf))
ELSE
ALLOCATE(ispot(nsurf))
END IF

IF (ALLOCATED(kpot)) THEN
DEALLOCATE(kpot)
ALLOCATE(kpot(50))
ELSE
ALLOCATE(kpot(50))
END IF

IF (ALLOCATED(kPotential)) THEN
DEALLOCATE(kPotential)
ALLOCATE(kPotential(500))
ELSE
ALLOCATE(kPotential(500))
END IF

ndependex = 0
ndependsurf = 0
ixdepend = 0
isdepend = 0
dependex = 0.0
dependsurf = 0.0d0

DO is = 1,nsurf
ispot(is) = is
END DO

ALLOCATE(stringarray(ncomp+nspec+nrct))
stringarray = ' '

DO k = 1,nkin
DO np = 1,nreactmin(k)
  ndependex(np,k) = 0
  ndependsurf(np,k) = 0
  DO kk = 1,ndepend(np,k)
    IF (namdep_nyf(kk,np,k) /= 'found') THEN
      speciesfound = .false.
!  Search through secondary aqueous species
      DO ksp = 1,nspec
        ik = ncomp + ksp
        IF (ulab(ik) == namdep_nyf(kk,np,k)) THEN
          idepend(kk,np,k) = ik
          speciesfound = .true.
        END IF
      END DO
      IF (speciesfound) THEN
        CYCLE
      END IF
!  SearcH through exchange species
      DO nex = 1,nexch_sec
        IF (nam_exchsec(nex) == namdep_nyf(kk,np,k)) THEN
          ndependex(np,k) = ndependex(np,k) + 1
          ncnt = ndependex(np,k)
          ixdepend(ncnt,np,k) = nex
          dependex(ncnt,np,k) = depend(kk,np,k)
          depend(kk,np,k) = 0.0
          speciesfound = .true.
        END IF
      END DO
      IF (speciesfound) THEN
        CYCLE
      END IF
!  Search through surface complexes
      DO is = 1,nsurf
        IF (namsurf(is) == namdep_nyf(kk,np,k)  ) THEN
          ndependsurf(np,k) = ndependsurf(np,k) + 1
          ncnt = ndependsurf(np,k)
          isdepend(ncnt,np,k) = is
          dependsurf(ncnt,np,k) = depend(kk,np,k)
          depend(kk,np,k) = 0.0
          speciesfound = .true.
        END IF
      END DO
      IF (speciesfound) THEN
        CYCLE
      END IF
      DO ns = 1,nsurf_sec
        IF (namsurf_sec(ns) == namdep_nyf(kk,np,k)  ) THEN
          ndependsurf(np,k) = ndependsurf(np,k) + 1
          ncnt = ndependsurf(np,k)
          isdepend(ncnt,np,k) = ns + nsurf
          dependsurf(ncnt,np,k) = depend(kk,np,k)
          depend(kk,np,k) = 0.0
          speciesfound = .true.
        END IF
      END DO
      IF (.NOT. speciesfound) THEN
        namtemp = namdep_nyf(kk,np,k)
        CALL stringlen(namtemp,ls)
        WRITE(*,*)
        WRITE(*,*) ' Species in mineral reaction not found'
        WRITE(*,*) ' Species: ',namtemp(1:ls)
        WRITE(*,*) ' In parallel reaction ', np
        namtemp = umin(k)
        CALL stringlen(namtemp,ls)
        WRITE(*,*) ' For mineral: ',namtemp(1:ls)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    END IF
  END DO    !  End of dependent species for a parallel reactions
  ndepend(np,k) = ndepend(np,k) - ndependex(np,k) - ndependsurf(np,k)
END DO      !  Loop over parallel reactions
END DO        !  Loop over minerals

!  ***************END OF ION EXCHANGE**************

!   ***********RADIOACTIVE DECAY SECTION**************

!section='decay'
!CALL readblock(nin,nout,section,found,ncount)

!IF (found) THEN

!  ALLOCATE(idecay(ncomp))
!  idecay = 0
!  ALLOCATE(nisotope(ncomp))
!  nisotope = 0
!  ALLOCATE(decay_label(30,ncomp))
!  decay_label = ' '
!  ALLOCATE(half_life(30,ncomp))
!  half_life = 0.0

!  WRITE(*,*) ' Radioactive decay block found'
!  CALL read_decay(nout,ncomp,ndecay)

!  nisotope_max = 0
!  DO id = 1,ndecay
!    nisotope_max = MAX(nisotope_max,nisotope(id))
!  END DO

!  Find the non-zero stoichiometric coefficient for radioactive species in minerals

!  ALLOCATE(kdecay(nrct,ndecay))
!  kdecay = 0
!  ALLOCATE(nmindecay(ndecay))
!  nmindecay = 0

!  DO id = 1,ndecay
!    i = idecay(id)     !  Point to primary species so as to sweep stoichiometric coeffs
!    kd = 0
!    DO k = 1,nrct
!      IF (mumin(1,k,i) /= 0.0) THEN
!        kd = kd + 1
!        kdecay(kd,id) = k    ! Point to mineral number
!      END IF
!    END DO
!    nmindecay(id) = kd    !  Set up so one can sweep "nmindecay" instead of "nrct"
!  END DO

!  nmindecay_max = 0
!  DO id = 1,ndecay
!    nmindecay_max = MAX(nmindecay_max,nmindecay(id))
!  END DO

!  i = size(idecay,1)
!  ALLOCATE(workint1(i))
!  workint1 = idecay
!  DEALLOCATE(idecay)
!  ALLOCATE(idecay(ndecay))
!  IF (ndecay /= 0) idecay(1:ndecay) = workint1(1:ndecay)
!  DEALLOCATE(workint1)

!  i = size(nisotope,1)
!  ALLOCATE(workint1(i))
!  workint1 = nisotope
!  DEALLOCATE(nisotope)
!  ALLOCATE(nisotope(ndecay))
!  IF (ndecay /= 0) nisotope(1:ndecay) = workint1(1:ndecay)
!  DEALLOCATE(workint1)

!  ndim1 = nmindecay_max
!  ndim2 = ndecay
!  i = size(kdecay,1)
!  j = size(kdecay,2)
!  ALLOCATE(workint2(i,j))
!  workint2 = kdecay
!  DEALLOCATE(kdecay)
!  ALLOCATE(kdecay(ndim1,ndim2))
!  IF(ndim1 /= 0 .AND. ndim2 /= 0)  &
!      kdecay(1:ndim1,1:ndim2) = workint2(1:ndim1,1:ndim2)
!  DEALLOCATE(workint2)

!  ndim1 = nisotope_max
!  ndim2 = ndecay
!  i = size(half_life,1)
!  j = size(half_life,2)
!  ALLOCATE(work2(i,j))
!  work2 = half_life
!  DEALLOCATE(half_life)
!  ALLOCATE(half_life(ndim1,ndim2))
!  IF(ndim1 /= 0 .AND. ndim2 /= 0)half_life(1:ndim1,1:ndim2) = work2(1:ndim1,1:ndim2)
!  DEALLOCATE(work2)

!  ndim1 = nisotope_max
!  ndim2 = ndecay
!  i = size(decay_label,1)
!  j = size(decay_label,2)
!  ALLOCATE(workchar2(i,j))
!  workchar2 = decay_label
!  DEALLOCATE(decay_label)
!  ALLOCATE(decay_label(ndim1,ndim2))
!  IF(ndim1 /= 0.AND.ndim2 /= 0)decay_label(1:ndim1,1:ndim2) = workchar2(1:ndim1,1:ndim2)
!  DEALLOCATE(workchar2)

!ELSE
!  WRITE(*,*) ' No radioactive decay block found'
!END IF

!     **************END OF RADIOACTIVE DECAY****************

!  Find the number of potentials (surface complexes using electrostatic correction)

kPotential = .FALSE.

DO is = 1,nsurf
k = ksurf(is)
IF (iedl(is) == 0) THEN
  kPotential(k) = .TRUE.
END IF
END DO

npot = 0
DO k = 1,nrct
IF (kPotential(k) .eqv.  .TRUE.) THEN
   npot = npot + 1
   kpot(npot) = k
END IF
END DO

IF (ALLOCATED(surfcharge_init)) THEN
DEALLOCATE(surfcharge_init)
ALLOCATE(surfcharge_init(nrct))
ELSE
ALLOCATE(surfcharge_init(nrct))
END IF
IF (ALLOCATED(LogPotential_tmp)) THEN
DEALLOCATE(LogPotential_tmp)
ALLOCATE(LogPotential_tmp(nsurf))
ELSE
ALLOCATE(LogPotential_tmp(nsurf))
END IF
IF (ALLOCATED(islink)) THEN
DEALLOCATE(islink)
ALLOCATE(islink(nsurf_sec))
ELSE
ALLOCATE(islink(nsurf_sec))
END IF
IF (ALLOCATED(nptlink)) THEN
DEALLOCATE(nptlink)
ALLOCATE(nptlink(nsurf_sec))
ELSE
ALLOCATE(nptlink(nsurf_sec))
END IF


surfcharge_init = 0.0
LogPotential_tmp = 0.0

!  Link the various secondary surface complexes to a primary surface hydroxyl site

DO ns = 1,nsurf_sec
DO is = 1,nsurf
  IF (musurf(ns,is+ncomp) /= 0.0) THEN
    islink(ns) = is
  END IF
END DO
END DO

nptlink = 0

DO ns = 1,nsurf_sec
DO npt = 1,npot
!!    IF (islink(ns) == ispot(npt)) THEN
!!      nptlink(ns) = npt
!!    END IF
 IF (ksurf(islink(ns)) == kpot(npt)) THEN
   nptlink(ns) = npt
 END IF
END DO
END DO

neqn = ncomp + nsurf + nexchange + npot

!  Temporary arrays deallocated later in START98

if (.not.allocated(SkipAdjust)) ALLOCATE(SkipAdjust(mchem))
SkipAdjust = .FALSE.
if (.not.allocated(tempcond)) ALLOCATE(tempcond(mchem))
if (.not.allocated(rocond)) ALLOCATE(rocond(mchem))
if (.not.allocated(porcond)) ALLOCATE(porcond(mchem))
if (.not.allocated(SaturationCond)) ALLOCATE(SaturationCond(mchem))
if (.not.allocated(PressureCond)) ALLOCATE(PressureCond(mchem))
if (.not.allocated(ctot)) ALLOCATE(ctot(ncomp,mchem))
if (.not.allocated(guess)) ALLOCATE(guess(ncomp,mchem))
if (.not.allocated(itype)) ALLOCATE(itype(ncomp,mchem))
if (.not.allocated(c_surf)) ALLOCATE(c_surf(nsurf,mchem))
if (.not.allocated(guess_surf)) ALLOCATE(guess_surf(nsurf,mchem))
if (.not.allocated(gaspp)) ALLOCATE(gaspp(ncomp,mchem))
if (.not.allocated(totexch)) ALLOCATE(totexch(nexchange,mchem))
IF (ALLOCATED(cec)) THEN
DEALLOCATE(cec)
ALLOCATE(cec(nexchange,mchem))
ELSE
ALLOCATE(cec(nexchange,mchem))
END IF
if (.not.allocated(ncon)) ALLOCATE(ncon(ncomp,mchem))
if (.not.allocated(condlabel)) ALLOCATE(condlabel(mchem))
if (.not.allocated(condtitle)) ALLOCATE(condtitle(mchem))
if (.not.allocated(fsurftmp)) ALLOCATE(fsurftmp(neqn,neqn))
if (.not.allocated(equilibrate)) ALLOCATE(equilibrate(nc,mchem))

!  **** End of temporary arrays (deallocated below)  ******

!  **** Permanent arrays  ***********************

IF (ALLOCATED(site_density)) THEN
DEALLOCATE(site_density)
ALLOCATE(site_density(nsurf,mchem))
ELSE
ALLOCATE(site_density(nsurf,mchem))
END IF
IF (ALLOCATED(specific)) THEN
DEALLOCATE(specific)
ALLOCATE(specific(nrct,mchem))
ELSE
ALLOCATE(specific(nrct,mchem))
END IF
IF (ALLOCATED(voltemp)) THEN
DEALLOCATE(voltemp)
ALLOCATE(voltemp(nrct,mchem))
ELSE
ALLOCATE(voltemp(nrct,mchem))
END IF
IF (ALLOCATED(scond)) THEN
DEALLOCATE(scond)
ALLOCATE(scond(ncomp,mchem))
ELSE
ALLOCATE(scond(ncomp,mchem))
END IF
IF (ALLOCATED(spcond)) THEN
DEALLOCATE(spcond)
ALLOCATE(spcond(ncomp+nspec,mchem))
ELSE
ALLOCATE(spcond(ncomp+nspec,mchem))
END IF
IF (ALLOCATED(spcond10)) THEN
DEALLOCATE(spcond10)
ALLOCATE(spcond10(ncomp+nspec,mchem))
ELSE
ALLOCATE(spcond10(ncomp+nspec,mchem))
END IF
IF (ALLOCATED(spcondgas)) THEN
DEALLOCATE(spcondgas)
ALLOCATE(spcondgas(ngas,mchem))
ELSE
ALLOCATE(spcondgas(ngas,mchem))
END IF
IF (ALLOCATED(spcondgas10)) THEN
DEALLOCATE(spcondgas10)
ALLOCATE(spcondgas10(ngas,mchem))
ELSE
ALLOCATE(spcondgas10(ngas,mchem))
END IF
IF (ALLOCATED(spcondex)) THEN
DEALLOCATE(spcondex)
ALLOCATE(spcondex(nexchange+nexch_sec,mchem))
ELSE
ALLOCATE(spcondex(nexchange+nexch_sec,mchem))
END IF
IF (ALLOCATED(spcondex10)) THEN
DEALLOCATE(spcondex10)
ALLOCATE(spcondex10(nexchange+nexch_sec,mchem))
ELSE
ALLOCATE(spcondex10(nexchange+nexch_sec,mchem))
END IF
IF (ALLOCATED(spcondsurf)) THEN
DEALLOCATE(spcondsurf)
ALLOCATE(spcondsurf(nsurf+nsurf_sec,mchem))
ELSE
ALLOCATE(spcondsurf(nsurf+nsurf_sec,mchem))
END IF
IF (ALLOCATED(spcondsurf10)) THEN
DEALLOCATE(spcondsurf10)
ALLOCATE(spcondsurf10(nsurf+nsurf_sec,mchem))
ELSE
ALLOCATE(spcondsurf10(nsurf+nsurf_sec,mchem))
END IF

IF (ALLOCATED(LogPotentialInit)) THEN
DEALLOCATE(LogPotentialInit)
ALLOCATE(LogPotentialInit(nsurf,mchem))
ELSE
ALLOCATE(LogPotentialInit(nsurf,mchem))
END IF

IF (ALLOCATED(volin)) THEN
DEALLOCATE(volin)
ALLOCATE(volin(nrct,mchem))
ELSE
ALLOCATE(volin(nrct,mchem))
END IF
volin = 0.0d0
IF (ALLOCATED(MineralMoles)) THEN
DEALLOCATE(MineralMoles)
ALLOCATE(MineralMoles(nrct,mchem))
ELSE
ALLOCATE(MineralMoles(nrct,mchem))
END IF
MineralMoles = 0.0d0
IF (ALLOCATED(areain)) THEN
DEALLOCATE(areain)
ALLOCATE(areain(nrct,mchem))
ELSE
ALLOCATE(areain(nrct,mchem))
END IF
IF (ALLOCATED(iarea)) THEN
DEALLOCATE(iarea)
ALLOCATE(iarea(nrct,mchem))
ELSE
ALLOCATE(iarea(nrct,mchem))
END IF
IF (ALLOCATED(fexch)) THEN
DEALLOCATE(fexch)
ALLOCATE(fexch(neqn,neqn))
ELSE
ALLOCATE(fexch(neqn,neqn))
END IF

IF (ALLOCATED(totex)) THEN
DEALLOCATE(totex)
ALLOCATE(totex(nexchange))
ELSE
ALLOCATE(totex(nexchange))
END IF
IF (ALLOCATED(sumactivity)) THEN
DEALLOCATE(sumactivity)
ALLOCATE(sumactivity(nexchange))
ELSE
ALLOCATE(sumactivity(nexchange))
END IF
IF (ALLOCATED(tec)) THEN
DEALLOCATE(tec)
ALLOCATE(tec(nexchange))
ELSE
ALLOCATE(tec(nexchange))
END IF
IF (ALLOCATED(wt_aexch)) THEN
DEALLOCATE(wt_aexch)
ALLOCATE(wt_aexch(nexchange))
ELSE
ALLOCATE(wt_aexch(nexchange))
END IF
IF (ALLOCATED(aexch)) THEN
DEALLOCATE(aexch)
ALLOCATE(aexch(nexch_sec))
ELSE
ALLOCATE(aexch(nexch_sec))
END IF
IF (ALLOCATED(TotChargeSave)) THEN
DEALLOCATE(TotChargeSave)
ALLOCATE(TotChargeSave(mchem))
ELSE
ALLOCATE(TotChargeSave(mchem))
END IF


!ALLOCATE(ratio_isotope_init(nisotope_max,nmindecay_max,ndecay,mchem))


equilibrate = .false.

!      ctot = 0.0001
!      guess = 0.0
!      itype = 1
!      gaspp= 0.0
!      spcond = 0.0
!      spcond10 = 0.0
!      spcondgas = 0.0
!      spcondgas10 = 0.0
!      spcondsurf = 0.0
!      spcondsurf10 = 0.0

!      c_surf= 0.0
!      guess_surf= 0.0
!      spcondex = 0.0
!      spcondex10 = 0.0
    totexch = 0.0
    cec = 0.0

    site_density = 0.0
    specific= 0.0
!!      iedl = 0

    tempcond = tinit

!      *****************POROSITY SECTION***********************

section = 'porosity'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
!!  WRITE(*,*)
!!  WRITE(*,*) ' Porosity parameters block found'
!!  WRITE(*,*)

parchar = 'fix_porosity'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter fix_porosity not found
  CONTINUE                ! Use mineral volume fractions
ELSE
  IF (realjunk > 0.0) THEN
    jpor = -1
    constantpor = realjunk
    WRITE(*,5010) constantpor
    WRITE(*,*) ' No update of porosity'
    GO TO 5011  ! If constantpor found, ignore porosity update
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Porosity should be greater than zero'
    WRITE(*,*) ' Porosity value: ', realjunk
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

parchar = 'minimum_porosity'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter minimum_porosity not found
  CONTINUE                ! Use mineral volume fractions
ELSE
  IF (realjunk > 0.0 .AND. realjunk <= 1.0) THEN
    MinimumPorosity = realjunk
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Minimum porosity should be greater than zero and less than or equal to 1'
    WRITE(*,*) ' Minimum porosity: ', realjunk
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

PoreFill = 0.0d0
parchar = 'porosity_exponent'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter minimum_porosity not found
  CONTINUE
ELSE
  IF (realjunk >= 0.0d0) THEN
    PoreFill = realjunk
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Porosity exponent should be greater than or equal to zero'
    WRITE(*,*) ' Porosity exponent: ', realjunk
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

PoreThreshold = 1.0d0
parchar = 'porosity_threshold'
parfind = ' '
realjunk = 0.0d0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter minimum_porosity not found
  CONTINUE                ! Use mineral volume fractions
ELSE
  IF (realjunk > 0.0 .AND. realjunk <= 1.0d0) THEN
    PoreThreshold = realjunk
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Porosity threshold should be greater than zero and less than 1.0'
    WRITE(*,*) ' Porosity threshold: ', realjunk
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END IF

parchar = 'read_porosityfile'
parfind = ' '
PorosityFile = ' '
CALL readFileName(nout,lchar,parchar,parfind,dumstring,section,PorosityFileFormat)
IF (parfind == ' ') THEN  !
  PorosityFile = ' '             ! Use default
!!  Check to make sure the user is not using the old designator "read_porosity"
  parchar = 'read_porosity'
  parfind = ' '
  CALL readFileName(nout,lchar,parchar,parfind,dumstring,section,PorosityFileFormat)
  IF (parfind == 'read_porosity') THEN
    WRITE(*,*)
    WRITE(*,*) 'Keyword "read_porosity" now obsolete--use "read_porosityfile"'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
ELSE
  PorosityFile = dumstring
  CALL stringlen(PorosityFile,ls)
  WRITE(*,*) ' Reading porosity from file: ',PorosityFile(1:ls)
!!!    jpor = 2
!!!    porosity_update = .FALSE.
!!!  New feature:  check to see if porosity update is set to true with read of porosity.
!!!	 If so, renormalize volume fractions to porosity that is read in.
  parchar = 'porosity_update'
  parfind = ' '
  porosity_update = .false.
  CALL read_logical(nout,lchar,parchar,parfind,porosity_update)
  IF (porosity_update) THEN
    jpor = 3
  ELSE
    jpor = 2
  END IF
END IF

IF (PorosityFile == ' ') THEN
  parchar = 'porosity_update'
  parfind = ' '
  porosity_update = .false.
  CALL read_logical(nout,lchar,parchar,parfind,porosity_update)
  IF (porosity_update) THEN
    jpor = 1
  ELSE
    jpor = 0
  END IF
END IF


5011   CONTINUE

ELSE

WRITE(*,*) ' Porosity parameters not found'
WRITE(*,*) ' Using defaults'
jpor = 1
MinimumPorosity = 1.0E-14

END IF


!  ****************GEOCHEMICAL INPUT****************

!  Now check here that appropriate species are present if a solute-based density module is to be used

MeanSalt = 0

IF (DensityModule == 'sodium_chloride') THEN
DO i = 1,ncomp+nspec
  IF (ulab(i) == 'Na+') then
    MeanSalt(1) = i
  END IF
END DO
DO i = 1,ncomp+nspec
  IF (ulab(i) == 'Cl-') then
    MeanSalt(2) = i
  END IF
END DO
IF (MeanSalt(1) == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Cation for mean solute concentration not found'
WRITE(*,*) ' ---> Could not find sodium in species list'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
IF (MeanSalt(2) == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Anion for mean solute concentration not found'
WRITE(*,*) ' ---> Could not find chloride in species list'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
ELSE IF (DensityModule == 'sodium_nitrate') THEN
DO i = 1,ncomp+nspec
  IF (ulab(i) == 'Na+') then
    MeanSalt(1) = i
  END IF
END DO
DO i = 1,ncomp+nspec
  IF (ulab(i) == 'NO3-') then
    MeanSalt(2) = i
  END IF
END DO
IF (MeanSalt(1) == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Cation for mean solute concentration not found'
WRITE(*,*) ' ---> Could not find sodium in species list'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
IF (MeanSalt(2) == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Anion for mean solute concentration not found'
WRITE(*,*) ' ---> Could not find nitrate in species list'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
ELSE IF (DensityModule == 'calcium_nitrate') THEN
DO i = 1,ncomp+nspec
  IF (ulab(i) == 'Ca++') then
    MeanSalt(1) = i
  END IF
END DO
DO i = 1,ncomp+nspec
  IF (ulab(i) == 'NO3-') then
    MeanSalt(2) = i
  END IF
END DO
IF (MeanSalt(1) == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Cation for mean solute concentration not found'
WRITE(*,*) ' ---> Could not find calcium in species list'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
IF (MeanSalt(2) == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Anion for mean solute concentration not found'
WRITE(*,*) ' ---> Could not find nitrate in species list'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
ELSE IF (DensityModule == 'potassium_nitrate') THEN
DO i = 1,ncomp+nspec
  IF (ulab(i) == 'K+') then
    MeanSalt(1) = i
  END IF
END DO
DO i = 1,ncomp+nspec
  IF (ulab(i) == 'NO3-') then
    MeanSalt(2) = i
  END IF
END DO
IF (MeanSalt(1) == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Cation for mean solute concentration not found'
WRITE(*,*) ' ---> Could not find potassium in species list'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
IF (MeanSalt(2) == 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' Anion for mean solute concentration not found'
WRITE(*,*) ' ---> Could not find nitrate in species list'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
ELSE IF (DensityModule == 'martin') THEN
DO i=1,ncomp+nspec
  if (ulab(i) == 'Tracer') then
    MeanSalt(1) = i
    MeanSalt(2) = i
  end if
END DO
ELSE
CONTINUE
CALL stringlen(DensityModule,ls)
WRITE(*,*)
WRITE(*,*) ' Current density module: ', DensityModule(1:ls)
WRITE(*,*)

END IF


!  Read in geochemical conditions

WRITE(iunit2,*)
WRITE(iunit2,*) ' ******  GEOCHEMICAL CONDITIONS INPUT  ******'
WRITE(iunit2,*)

WRITE(*,*)
WRITE(*,*) ' ---> READING GEOCHEMICAL CONDITIONS  '
WRITE(*,*)

REWIND nin
!DO i = 1,ncomp+nspec
!  WRITE(*,*) ulab(i)
!END DO

IF (ALLOCATED(SolidDensityFrom)) THEN
DEALLOCATE(SolidDensityFrom)
ALLOCATE(SolidDensityFrom(mchem))
ELSE
ALLOCATE(SolidDensityFrom(mchem))
END IF
IF (ALLOCATED(SolidDensity)) THEN
DEALLOCATE(SolidDensity)
ALLOCATE(SolidDensity(mchem))
ELSE
ALLOCATE(SolidDensity(mchem))
END IF
IF (ALLOCATED(SolidSolutionRatio)) THEN
DEALLOCATE(SolidSolutionRatio)
ALLOCATE(SolidSolutionRatio(mchem))
ELSE
ALLOCATE(SolidSolutionRatio(mchem))
END IF
IF (ALLOCATED(OneOverMassFraction)) THEN
DEALLOCATE(OneOverMassFraction)
ALLOCATE(OneOverMassFraction(mchem))
ELSE
ALLOCATE(OneOverMassFraction(mchem))
END IF
IF (ALLOCATED(conversion)) THEN
DEALLOCATE(conversion)
ALLOCATE(conversion(mchem))
ELSE
ALLOCATE(conversion(mchem))
END IF
IF (ALLOCATED(unitsflag)) THEN
DEALLOCATE(unitsflag)
ALLOCATE(unitsflag(mchem))
ELSE
ALLOCATE(unitsflag(mchem))
END IF

IF (Duan .OR. Duan2006) THEN
IF (ALLOCATED(GasPressureTotalInit)) THEN
  DEALLOCATE(GasPressureTotalInit)
  ALLOCATE(GasPressureTotalInit(mchem))
ELSE
  ALLOCATE(GasPressureTotalInit(mchem))
END IF
END IF

unitsflag = 1

SaturationCond = FixSaturation

CALL find_condition(nin,nout,found,phfound,ncomp,  &
  nspec,nrct,nkin,ngas,nexchange,nsurf,ndecay,   &
  ph,guessph,constraint,nchem,unitsflag,jpor,    &
  DensityModule,RunningPest)

IF (nchem > mchem) THEN
WRITE(*,*)
WRITE(*,*) '  Nchem dimensioned too small'
WRITE(*,*)
READ(*,*)
STOP
END IF
IF (nexchange > mexch) THEN
WRITE(*,*)
WRITE(*,*) ' Mexch dimensioned too small'
WRITE(*,*) ' Nexchange = ',nexchange
WRITE(*,*)
READ(*,*)
STOP
END IF

IF (nchem < 1) THEN
WRITE(*,*) 'You have not specified any geochemical conditions'
READ(*,*)
STOP
END IF
IF (nchem > 1) THEN
!  WRITE(*,*) 'Number of geochemical conditions specified = ', nchem
WRITE(iunit2,*) 'Number of geochemical conditions specified = ', nchem
END IF

DO k = 1,nchem

!  DO i = 1,ncomp
!    WRITE(*,*)
!    WRITE(*,*) ' Species = ',ulab(i)
!    WRITE(*,*) ' itype = ',itype(i,k)
!    WRITE(*,*)
!  END DO

dumstring = condlabel(k)
CALL stringlen(dumstring,ls)
!  WRITE(*,*)
!  WRITE(*,*) ' Condition Label: ',dumstring(1:ls)
!  WRITE(*,*) ' Condition title ',condtitle(k)
!  WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,*) ' Condition Label: ',dumstring(1:ls)
!  WRITE(iunit2,*) ' Condition title ',condtitle(k)
WRITE(iunit2,*)

!  WRITE(*,*)
!  WRITE(*,596)
!  DO i = 1,ncomp
!    WRITE(*,595) ulab(i),itype(i,k),guess(i,k), ctot(i,k),constraint(i,k)
!  END DO

WRITE(iunit2,*)
WRITE(iunit2,596)
DO i = 1,ncomp
  WRITE(iunit2,595) ulab(i),itype(i,k),guess(i,k), ctot(i,k),constraint(i,k)
END DO

portemp = porcond(k)

WRITE(iunit2,*)
WRITE(iunit2,599) portemp
WRITE(iunit2,*)
WRITE(iunit2,*)
!  WRITE(*,*)
!  WRITE(*,599) portemp
!  WRITE(*,*)
!  WRITE(*,*)

DO i = 1,ncomp
  ncon(i,k) = constraint(i,k)
END DO

END DO

!  Find master species and pH (if present)

ikmast = 0
DO ik = 1,ncomp+nspec
IF(ulab(ik) == master) THEN
  ikmast = ik
END IF
END DO
IF (ikmast == 0) THEN
ikmast = 1
END IF

dumstring = ulab(ikmast)
CALL stringlen(dumstring,ls)

!WRITE(*,*)
!WRITE(*,*) ' Master variable = ', dumstring(1:ls)
!WRITE(*,*)


ikph = 0
ikFe2 = 0
ikFe3 = 0
ikO2 = 0
DO ik = 1,ncomp+nspec

IF(ulab(ik) == 'H+' .OR. ulab(ik) == 'h+') THEN
  ikph = ik
END IF
IF (ulab(ik) == 'Fe++') THEN
  ikFe2 = ik
END IF
IF (ulab(ik) == 'Fe+++') THEN
  ikFe3 = ik
END IF
IF(ulab(ik) == 'O2(aq)' .OR. ulab(ik) == 'o2(aq)') THEN
  ikO2 = ik
END IF
IF (ulab(ik) == 'Na+') THEN
  ikNa = ik
END IF
IF (ulab(ik) == 'Cl-') THEN
  ikCl= ik
END IF
END DO

!!IF (ikFe2 /= 0 .and. ikFe3 /= 0 .and. iko2 /= 0) THEN


ikTracer = 0
DO ik = 1,ncomp+nspec
IF(ulab(ik) == 'Tracer' .OR. ulab(ik) == 'Tracer(aq)' .or. ulab(ik)=='tracer') THEN
  ikTracer = ik
END IF
END DO

H2Opresent = .FALSE.
DO i = 1,ncomp
IF(ulab(i) == 'H2O' .OR. ulab(ik) == 'h2o') THEN
  H2Opresent = .TRUE.
END IF
END DO


!  *****************PEST BLOCK***********************

section = 'pest'
CALL readblock(nin,nout,section,found,ncount)

NPestExchange = 0
NPestSurface = 0

IF (found) THEN

pest = .TRUE.

CreatePestInstructionFile = .FALSE.
parchar = 'CreatePestInstructionFile'
parfind = ' '
CALL read_logical(nout,lchar,parchar,parfind,CreatePestInstructionFile)

parchar = 'createpestexchangefile'
parfind = ' '
CALL ReadFileNameOnly(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN
 PestExchangeOutputFile = 'PestExchange.out'             ! Use default
ELSE
 PestExchangeOutputFile = dumstring
END IF

CALL PestExchange(nout,ncomp,nchem)

parchar = 'createpestsurfacefile'
parfind = ' '
CALL ReadFileNameOnly(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN
 PestSurfaceOutputFile = 'PestSurface.out'             ! Use default
ELSE
 PestSurfaceOutputFile = dumstring
END IF

CALL PestSurface(nout,ncomp,nchem)

ELSE
pest = .FALSE.
END IF

!  *****************ISOTOPES BLOCK***********************

IF (ALLOCATED(IsotopeMineralRare)) THEN
DEALLOCATE(IsotopeMineralRare)
END IF
ALLOCATE(IsotopeMineralRare(nrct))
IF (ALLOCATED(IsotopeMineralCommon)) THEN
DEALLOCATE(IsotopeMineralCommon)
END IF
ALLOCATE(IsotopeMineralCommon(nrct))

IsotopeMineralRare = .FALSE.
IsotopeMineralCommon = .FALSE.

section = 'isotopes'
CALL readblock(nin,nout,section,found,ncount)

IF (ALLOCATED(IsotopePrimaryRare)) THEN
DEALLOCATE(IsotopePrimaryRare)
END IF
ALLOCATE(IsotopePrimaryRare(ncomp))
IF (ALLOCATED(IsotopePrimaryCommon)) THEN
DEALLOCATE(IsotopePrimaryCommon)
END IF
ALLOCATE(IsotopePrimaryCommon(ncomp))
IsotopePrimaryRare = 0
IsotopePrimaryCommon = 0

IF (found) THEN

CALL read_Isotopes(nout,ncomp,nspec,nrct,nsurf,nexchange,ngas)

NoFractionationDissolution = .FALSE.
parchar = 'NoFractionationDissolution'
parfind = ' '
CALL read_logical(nout,lchar,parchar,parfind,NoFractionationDissolution)

CALL isotope_read_series(nout,nx,ny,nz)

i = size(jxisotopeseries,1)
ALLOCATE(workint1(i))
workint1 = jxisotopeseries
DEALLOCATE(jxisotopeseries)
ALLOCATE(jxisotopeseries(nisotopeseries))
IF(nisotopeseries /= 0) jxisotopeseries(1:nisotopeseries) = workint1(1:nisotopeseries)
DEALLOCATE(workint1)

i = size(jyisotopeseries,1)
ALLOCATE(workint1(i))
workint1 = jyisotopeseries
DEALLOCATE(jyisotopeseries)
ALLOCATE(jyisotopeseries(nisotopeseries))
IF(nisotopeseries /= 0)jyisotopeseries(1:nisotopeseries) = workint1(1:nisotopeseries)
DEALLOCATE(workint1)

i = size(jzisotopeseries,1)
ALLOCATE(workint1(i))
workint1 = jzisotopeseries
DEALLOCATE(jzisotopeseries)
ALLOCATE(jzisotopeseries(nisotopeseries))
IF(nisotopeseries /= 0)jzisotopeseries(1:nisotopeseries) = workint1(1:nisotopeseries)
DEALLOCATE(workint1)

ELSE


END IF
!!   ********* END OF ISOTOPES BLOCK ********************

! skip what follows if alquimia is defined
#ifndef ALQUIMIA
!     ********SPECIATION OF GEOCHEMICAL CONDITIONS******

!  First, call the initialization routine so that the
!  the geochemical conditions can be written to the output file

WRITE(*,*)
WRITE(*,*) ' ---> STARTING SPECIATION OF GEOCHEMICAL CONDITIONS  '
WRITE(*,*)

WRITE(iunit2,*)
WRITE(iunit2,*) '  ********  SPECIATION OF GEOCHEMICAL CONDITIONS  ********'
WRITE(iunit2,*)

ALLOCATE(AqueousToBulkCond(nchem))

iinit = 1
DO nco = 1,nchem

dumstring = condlabel(nco)
CALL stringlen(dumstring,ls)

tempc = tempcond(nco)

!!  sumpor = 0.0
!!  DO k = 1,nkin
!!    sumpor = sumpor + volin(k,nco)
!!  END DO
!!  IF (jpor == -1) THEN
!!    portemp = constantpor
!!  ELSE
!!    portemp = 1.0 - sumpor
!!  END IF

portemp = porcond(nco)
PressureTemp = PressureCond(nco)

CALL keqcalc2_init(ncomp,nrct,nspec,ngas,nsurf_sec,tempc,PressureTemp)

DO i = 1,ncomp
  namtemp = ulab(i)
  sptmp10(i) = guess(i,nco)
  sptmp(i) = DLOG(sptmp10(i))
END DO
DO ix = 1,nexchange
  spextmp10(ix) = 1.0
  spextmp(ix) = DLOG(spextmp10(ix))
END DO
DO is = 1,nsurf
  spsurftmp10(is) = guess_surf(is,nco)
  spsurftmp(is) = DLOG(spsurftmp10(is))
END DO

gamtmp = 0.0

LogPotential_tmp = 0.0
spgastmp = -100.0d0
spgastmp10 = 1.0D-35

CALL species_init(ncomp,nspec)
CALL gases_init(ncomp,ngas,tempc)
CALL surf_init(ncomp,nspec,nsurf,nsurf_sec,nchem)
CALL exchange_init(ncomp,nspec,nexchange,nexch_sec,nchem)
CALL totconc_init(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,nco)
CALL totgas_init(ncomp,nspec,ngas)
CALL totsurf_init(ncomp,nsurf,nsurf_sec)
CALL totexchange_init(ncomp,nexchange,nexch_sec,nco)

IF(.NOT. RunningPest) THEN
  WRITE(*,*)
  WRITE(*,*) ' --> Initializing Condition:  ',dumstring(1:ls)
  WRITE(*,*)
END IF

CALL equilib_co2(ncomp,nspec,nrct,ngas,nsurf,igamma,ikph,  &
    nco,nexchange,nexch_sec,nsurf_sec,npot,neqn,tempc,portemp,    &
    DensityModule,ipest,PestUnit,pest)

!! Convert from bars pressure to n/V (mol/m^3) by converting to Pascals, then dividing by RT
tk = tempc + 273.15d0
denmol = 1.e05/(8.314*tk)   ! P/RT = n/V, with pressure converted from bars to Pascals

spgastmp10 = spgastmp10*denmol
spgastmp = DLOG(spgastmp10)

DO ik = 1,ncomp+nspec
  spcond(ik,nco) = sptmp(ik)
  spcond10(ik,nco) = sptmp10(ik)
END DO
DO kk = 1,ngas
  spcondgas(kk,nco) = spgastmp(kk)
  spcondgas10(kk,nco) = spgastmp10(kk)
END DO

DO i = 1,ncomp
  scond(i,nco) = stmp(i)
END DO

IF (DensityModule /= 'temperature') THEN
!   Calculate the correction for the mass fraction of water:  kg_solution/kg_water
  MeanSaltConcentration = 0.001*(wtaq(MeanSalt(1))*scond(MeanSalt(1),nco) +   &
          wtaq(MeanSalt(2))*scond(MeanSalt(2),nco))
  MassFraction = 1.0/(1.0 + MeanSaltConcentration)
ELSE
  MassFraction = 1.0
END IF

AqueousToBulkCond(nco) = rocond(nco)*SaturationCond(nco)*porcond(nco)*MassFraction

DO nex = 1,nexchange+nexch_sec
  spcondex(nex,nco) = spextmp(nex)
  spcondex10(nex,nco) = spextmp10(nex)
END DO
DO is = 1,nsurf
  spcondsurf(is,nco) = spsurftmp(is)
  spcondsurf10(is,nco) = spsurftmp10(is)
  IF (iedl(is) == 0) THEN     !!  Electrostatic option
    LogPotentialInit(is,nco) = LogPotential_tmp(is)
  END IF
END DO
DO ns = 1,nsurf_sec
  is = ns + nsurf
  spcondsurf(is,nco) = spsurftmp(is)
  spcondsurf10(is,nco) = spsurftmp10(is)
END DO

!  Map the species to an array dimensioned to number of geochemical conditions

END DO


!  If ICOMPLETE = 1 (DATABASE SWEEP), STOP HERE

IF (icomplete == 1) THEN
WRITE(*,*)
WRITE(*,*) '    *****DATABASE SWEEP SPECIFIED*****'
WRITE(*,*) '  **PROGRAM STOPS AFTER INITIALIZATION**'
WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,*) ' *****DATABASE SWEEP SPECIFIED*****'
WRITE(iunit2,*) ' **PROGRAM STOPS AFTER INITIALIZATION**'
WRITE(iunit2,*)
WRITE(*,*)
WRITE(*,*) '               SPECIATION OF '
WRITE(*,*) '     INITIAL AND BOUNDARY CONDITIONS '
WRITE(*,*) '          SUCCESSFULLY COMPLETED'
WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,*) '               SPECIATION OF  '
WRITE(iunit2,*) '     INITIAL AND BOUNDARY CONDITIONS '
WRITE(iunit2,*) '          SUCCESSFULLY COMPLETED'
WRITE(iunit2,*)
STOP
END IF

WRITE(*,*)
WRITE(*,*) '               SPECIATION OF '
WRITE(*,*) '     INITIAL AND BOUNDARY CONDITIONS '
WRITE(*,*) '          SUCCESSFULLY COMPLETED'
WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,*) '               SPECIATION OF  '
WRITE(iunit2,*) '     INITIAL AND BOUNDARY CONDITIONS '
WRITE(iunit2,*) '          SUCCESSFULLY COMPLETED'
WRITE(iunit2,*)

IF (ispeciate == 1) STOP

#endif
! end of block to skip for ALQUIMIA

!  ***************STOP HERE WHEN DATABASE SWEEP IS DONE***************

!  Check dimensions for transport (beyond database sweep)

IF (ncomp > mcomp) THEN
WRITE(*,*)
WRITE(*,*) ' Primary species dimension too small! ',ncomp
WRITE(*,*) ' Redimension MCOMP in params.inc'
WRITE(*,*)
READ(*,*)
STOP
END IF
IF (nspec > mspec) THEN
WRITE(*,*)
WRITE(*,*) ' Secondary species dimension too small! ',nspec
WRITE(*,*) ' Redimension MSPEC in params.inc'
WRITE(*,*)
READ(*,*)
STOP
END IF
IF (nrct > mrct) THEN
WRITE(*,*)
WRITE(*,*) ' Mineral dimension too small! ',nrct
WRITE(*,*) ' Redimension MRCT in params.inc'
WRITE(*,*)
READ(*,*)
STOP
END IF
IF (ngas > mgas) THEN
WRITE(*,*)
WRITE(*,*) ' Gas dimension too small! ',ngas
WRITE(*,*) ' Redimension MGAS in params.inc'
WRITE(*,*)
READ(*,*)
STOP
END IF
IF (nexchange > mexch) THEN
WRITE(*,*)
WRITE(*,*) ' Exchange site dimension too small! ',nexchange
WRITE(*,*) ' Redimension MEXCH in params.inc'
WRITE(*,*)
READ(*,*)
STOP
END IF
IF (nexch_sec > mexch_sec) THEN
WRITE(*,*)
WRITE(*,*) ' Exchange secondary species dimension too small! ',nexch_sec
WRITE(*,*) ' Redimension MEXCH_SEC in params.inc'
WRITE(*,*)
READ(*,*)
STOP
END IF
IF (nsurf > msurf) THEN
WRITE(*,*)
WRITE(*,*) ' Surface hydroxyl dimension too small! ',nsurf
WRITE(*,*) ' Redimension MSURF in params.inc'
WRITE(*,*)
READ(*,*)
STOP
END IF
IF (nsurf_sec > msurf_sec) THEN
WRITE(*,*)
WRITE(*,*) ' Surface secondary species dimension too small! ',nsurf_sec
WRITE(*,*) ' Redimension MSURF_SEC in params.inc'
WRITE(*,*)
READ(*,*)
STOP
END IF

IF (ALLOCATED(halfsataq)) THEN
DEALLOCATE(halfsataq)
ALLOCATE(halfsataq(ncomp+nspec,maqkin))
ELSE
ALLOCATE(halfsataq(ncomp+nspec,maqkin))
END IF
IF (ALLOCATED(termMonod)) THEN
DEALLOCATE(termMonod)
ALLOCATE(termMonod(ncomp+nspec,maqkin))
ELSE
ALLOCATE(termMonod(ncomp+nspec,maqkin))
END IF
IF (ALLOCATED(rinhibitaq)) THEN
DEALLOCATE(rinhibitaq)
ALLOCATE(rinhibitaq(ncomp+nspec,maqkin))
ELSE
ALLOCATE(rinhibitaq(ncomp+nspec,maqkin))
END IF
IF (ALLOCATED(imonodaq)) THEN
DEALLOCATE(imonodaq)
ALLOCATE(imonodaq(ncomp+nspec,maqkin))
ELSE
ALLOCATE(imonodaq(ncomp+nspec,maqkin))
END IF
IF (ALLOCATED(inhibitaq)) THEN
DEALLOCATE(inhibitaq)
ALLOCATE(inhibitaq(ncomp+nspec,maqkin))
ELSE
ALLOCATE(inhibitaq(ncomp+nspec,maqkin))
END IF
IF (ALLOCATED(itot_monodaq)) THEN
DEALLOCATE(itot_monodaq)
ALLOCATE(itot_monodaq(ncomp,maqkin))
ELSE
ALLOCATE(itot_monodaq(ncomp,maqkin))
END IF
IF (ALLOCATED(itot_inhibitaq)) THEN
DEALLOCATE(itot_inhibitaq)
ALLOCATE(itot_inhibitaq(ncomp,maqkin))
ELSE
ALLOCATE(itot_inhibitaq(ncomp,maqkin))
END IF
IF (ALLOCATED(iaqtype)) THEN
DEALLOCATE(iaqtype)
ALLOCATE(iaqtype(maqkin))
ELSE
ALLOCATE(iaqtype(maqkin))
END IF
IF (ALLOCATED(nmonodaq)) THEN
DEALLOCATE(nmonodaq)
ALLOCATE(nmonodaq(maqkin))
ELSE
ALLOCATE(nmonodaq(maqkin))
END IF
IF (ALLOCATED(ninhibitaq)) THEN
DEALLOCATE(ninhibitaq)
ALLOCATE(ninhibitaq(maqkin))
ELSE
ALLOCATE(ninhibitaq(maqkin))
END IF
IF (ALLOCATED(keqkin)) THEN
DEALLOCATE(keqkin)
ALLOCATE(keqkin(maqkin))
ELSE
ALLOCATE(keqkin(maqkin))
END IF
IF (ALLOCATED(mukin)) THEN
DEALLOCATE(mukin)
ALLOCATE(mukin(maqkin,ncomp))
ELSE
ALLOCATE(mukin(maqkin,ncomp))
END IF

halfsataq = 0.0
rinhibitaq = 0.0
iaqtype = 0
nmonodaq = 0
ninhibitaq = 0
imonodaq = 0
inhibitaq = 0
itot_monodaq = 0
itot_inhibitaq = 0
keqkin = 0.0
mukin = 0.0


!         ***********AQUEOUS KINETICS SECTION**************

ikin = 0
nreactkinmax = 0
section='aqueous_kinetics'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN

WRITE(*,*) ' Aqueous kinetics block found'
!!CALL read_kinetics(nout,ncomp,nspec,nrct,ikin,data1)
CALL read_kinetics_Bio(nout,ncomp,nspec,nrct,ikin,nkin,data2)

!  If radioactive decay equations are present, check minerals for corresponding
!    radiogenic isotopes

nmonodaqmax = 0
ninhibitaqmax = 0
DO ir = 1,ikin
  nreactkinmax = MAX0(nreactkinmax,nreactkin(ir))
  DO kk = 1,nmonodaq(ir)
    nmonodaqmax = MAX0(nmonodaqmax,nmonodaq(ir))
    ninhibitaqmax = MAX0(ninhibitaq(ir),ninhibitaqmax)
  END DO
END DO

IF (ALLOCATED(nrad_decay)) THEN
  DEALLOCATE(nrad_decay)
  ALLOCATE(nrad_decay(ikin))
ELSE
  ALLOCATE(nrad_decay(ikin))
END IF
IF (ALLOCATED(iraddecay)) THEN
  DEALLOCATE(iraddecay)
  ALLOCATE(iraddecay(ikin))
ELSE
  ALLOCATE(iraddecay(ikin))
END IF

iraddecay = 0
nrad_decay = 0

ALLOCATE(list_tmp(100))
ALLOCATE(krad_tmp(50,ikin))
ALLOCATE(nprad_tmp(50,ikin))


DO ir = 1,ikin
  IF (iaqtype(ir) == 4) THEN
    DO i = 1,ncomp
      IF (dependk(i,1,ir) == 1.0) then
        iraddecay(ir) = i            !  This flags the isotopic parent in radioactive decay rxn
      END IF
    END DO
  END IF
END DO

DO ir = 1,ikin
  IF (iaqtype(ir) == 4) THEN
    i = iraddecay(ir)
    DO k = 1,nrct
      DO np = 1,nreactmin(k)
        IF (mumin(np,k,i) > 0.0) THEN
          dumstring = umin(k)
          CALL stringlen(dumstring,ls)
          WRITE(*,*)
          WRITE(*,*) ' Mineral found with radioactive decay reaction: ', dumstring(1:ls)
          WRITE(*,*) ' Parallel reaction:  ', np
          WRITE(*,*)
          nrad_decay(ir) = nrad_decay(ir) + 1
          kd = nrad_decay(ir)
          krad_tmp(kd,ir) = k
          nprad_tmp(kd,ir) = np
          irsave = ir
!           Find the daughter element (if there is one and if it hasn't been loaded yet)
5103      continue
          idaughter = 0
          DO i2 = 1,ncomp
            IF (mukin(irsave,i2) == 1.0) THEN
              idaughter= i2
            END IF
          END DO
          IF (idaughter /= 0 .AND. mumin(np,k,idaughter) == 0.0) THEN
            DaughterFound = .FALSE.
            DO ir2 = 1,ikin
              IF (idaughter == iraddecay(ir2)) THEN
                dumstring = ulab(idaughter)
                CALL stringlen(dumstring,ls)
                WRITE(*,*)
                WRITE(*,*) ' Decay reaction involving daughter found: ', dumstring(1:ls)
                WRITE(*,*)
                nrad_decay(ir2) = nrad_decay(ir2) + 1
                kd = nrad_decay(ir2)
                krad_tmp(kd,ir2) = k
                nprad_tmp(kd,ir2) = np
                irsave = ir2
                DaughterFound = .TRUE.
              END IF
            END DO
            IF (DaughterFound) THEN
              GO TO 5103
            ELSE
              GO TO 5104
            END IF
          ELSE
            GO TO 5104
          END IF
!            GO TO 5103
5104      CONTINUE
        END IF
      END DO
    END DO
  END IF
END DO

nradmax = 0
DO ir = 1,ikin
  nradmax = MAX0(nradmax,nrad_decay(ir))
END DO

IF (ALLOCATED(kradpoint)) THEN
  DEALLOCATE(kradpoint)
  ALLOCATE(kradpoint(nradmax,ikin))
ELSE
  ALLOCATE(kradpoint(nradmax,ikin))
END IF
IF (ALLOCATED(npradpoint)) THEN
  DEALLOCATE(npradpoint)
  ALLOCATE(npradpoint(nradmax,ikin))
ELSE
  ALLOCATE(npradpoint(nradmax,ikin))
END IF

kradpoint = 0
npradpoint = 1

DO ir = 1,ikin
  DO kd = 1,nrad_decay(ir)
    kradpoint(kd,ir) = krad_tmp(kd,ir)
    npradpoint(kd,ir) = nprad_tmp(kd,ir)
  END DO
END DO

DEALLOCATE(list_tmp)
DEALLOCATE(krad_tmp)
DEALLOCATE(nprad_tmp)

ELSE
WRITE(*,*) ' No aqueous kinetics block found'
END IF



!     **************END OF AQUEOUS KINETICS****************
!   *************  Kd SECTION  ************************

IF (ALLOCATED(distrib)) THEN
DEALLOCATE(distrib)
END IF
ALLOCATE(distrib(ncomp))
distrib = 0.0

section = 'retardation'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
!!  WRITE(*,*)
!!  WRITE(*,*) ' Retardation parameters block found'
!!  WRITE(*,*)


CALL read_kd(nout,ncomp,nretard)

ELSE

WRITE(*,*) ' Retardation parameters not found'
WRITE(*,*) ' Assuming NO retardation (Kd = 0) '

END IF


!   ***************************************************

! skip what follows if alquimia is defined
#ifndef ALQUIMIA
!************DISCRETIZATION****************************
!  Check for discretization block
!  If absent, assume a reaction path calculation
!******************************************************

section = 'discretization'
CALL readblock(nin,nout,section,found,ncount)

OutputDistanceScale = 1.0

IF (found) THEN
!!  WRITE(*,*) ' Discretization block found'
!!  WRITE(*,*)

dist_scale = 1.0d0
CALL units_distance(nout,section,dist_scale)

OutputDistanceScale = dist_scale

!*****************
!  X discretization
!******************

parchar = 'xzones'
parfind = ' '
realmult = 0.0
CALL read_multpar(nout,lchar,parchar,parfind,realmult,lenarray,section)
IF (parfind == ' ') THEN
!          write(*,*) ' No discretization in X direction'
!          write(*,*) ' Setting nx = 1'
  nzonex = 0
  nx = 1
ELSE
  IF (realmult(1) <= 0.0) THEN
    WRITE(*,*) ' X discretization specified, but no info given'
    WRITE(*,*) ' Setting nx = 1'
    nzonex = 0
    nx = 1
  ELSE
!            write(*,*) ' X discretization'
!  Check to see if there are an even number values given
!  One for the number of cells, the second for the spacing
    IF (MOD(lenarray,2) == 0) THEN
      nzonex = lenarray/2
      DO i = 1,lenarray,2
        ii = (i+1)/2
        nvx(ii) = realmult(i)
        dxxt(ii) = realmult(i+1)
!                write(*,*) nvx(ii),dxxt(ii)
      END DO
!              write(*,*)  ' Number of X zones = ',nzonex
    ELSE
      lenarray = lenarray - 1
      IF (lenarray > 1) THEN
        WRITE(*,*) ' Both the number of cells and the grid'
        WRITE(*,*) '   spacing should be specified--'
        WRITE(*,*) '   Using only complete pairs'
        nzonex = lenarray/2
        WRITE(*,*)  ' Number of X zones = ',nzonex
        WRITE(*,*)
        DO i = 1,lenarray,2
          ii = (i+1)/2
          nvx(ii) = realmult(i)
          dxxt(ii) = realmult(i+1)
!                  write(*,*) nvx(ii),dxxt(ii)
        END DO
      ELSE
        WRITE(*,*) ' Both the number of cells and the grid'
        WRITE(*,*) '   spacing should be specified'
        WRITE(*,*) ' Only one value provided'
        READ(*,*)
        STOP
      END IF
    END IF
  END IF
END IF

!!!  End of X discretiztion
!!!  ********************************************



!!!IF (.NOT. ReadGridVolumes) THEN
!*****************
!  Y discretization
!******************

parchar = 'yzones'
parfind = ' '
realmult = 0.0
CALL read_multpar(nout,lchar,parchar,parfind,realmult,lenarray,section)
IF (parfind == ' ') THEN
!          write(*,*) ' No discretization in Y direction'
!          write(*,*) ' Setting ny = 1'
!          write(*,*)
  nzoney = 0
  ny = 1
ELSE
  IF (realmult(1) <= 0.0) THEN
    WRITE(*,*) ' Y discretization specified, but no info given'
    WRITE(*,*) ' Setting ny = 1'
    WRITE(*,*)
    nzoney = 0
    ny = 1
  ELSE
!            write(*,*) ' Y discretization'
!  Check to see if there are an even number values given
!  One for the number of cells, the second for the spacing
    IF (MOD(lenarray,2) == 0) THEN
      nzoney = lenarray/2
      DO i = 1,lenarray,2
        ii = (i+1)/2
        nvy(ii) = realmult(i)
        dyyt(ii) = realmult(i+1)
!                write(*,*) nvy(ii),dyyt(ii)
      END DO
!              write(*,*)  ' Number of Y zones = ',nzoney
    ELSE
      lenarray = lenarray - 1
      IF (lenarray > 1) THEN
        WRITE(*,*) ' Both the number of cells and the grid'
        WRITE(*,*) '   spacing should be specified--'
        WRITE(*,*) '   Using only complete pairs'
        nzoney = lenarray/2
        WRITE(*,*)  ' Number of Y zones = ',nzoney
        WRITE(*,*)
        DO i = 1,lenarray,2
          ii = (i+1)/2
          nvy(ii) = realmult(i)
          dyyt(ii) = realmult(i+1)
!                  write(*,*) nvy(ii),dyyt(ii)
        END DO
      ELSE
        WRITE(*,*) ' Both the number of cells and the grid'
        WRITE(*,*) '   spacing should be specified'
        WRITE(*,*) ' Only one value provided'
        READ(*,*)
        STOP
      END IF
    END IF
  END IF
END IF

!*****************
!  Z discretization
!******************

parchar = 'zzones'
parfind = ' '
realmult = 0.0
CALL read_multpar(nout,lchar,parchar,parfind,realmult,lenarray,section)

      if (parfind.eq.' ') then
        write(*,*) ' No discretization in Z direction'
        write(*,*) ' Setting nz = 1'
        write(*,*)
        nzonez = 0
        nz = 1
      else
        if (realmult(1).le.0.0) then
          write(*,*) ' Z discretization specified, but no info given'
          write(*,*) ' Setting nz = 1'
          write(*,*)
          nzonez = 0
          nz = 1
        else
          write(*,*) ' Z discretization'

!  Check to see if there are an even number values given
!  One for the number of cells, the second for the spacing

         if (mod(lenarray,2).eq.0) then

            nzonez = lenarray/2
            do i = 1,lenarray,2
              ii = (i+1)/2
              nvz(ii) = realmult(i)
              dzzt(ii) = realmult(i+1)
            end do

         else

            lenarray = lenarray - 1
            if (lenarray.gt.1) then
              write(*,*) ' Both the number of cells and the grid'
              write(*,*) '   spacing should be specified--'
              write(*,*) '   Using only complete pairs'
              nzonez = lenarray/2
              write(*,*)  ' Number of Z zones = ',nzonez
              write(*,*)
              do i = 1,lenarray,2
                ii = (i+1)/2
                nvz(ii) = realmult(i)
                dzzt(ii) = realmult(i+1)
                write(*,*) nvz(ii),dzzt(ii)
              end do

            else

              write(*,*) ' Both the number of cells and the grid'
              write(*,*) '   spacing should be specified'
              write(*,*) ' Only one value provided'
              READ(*,*)
              STOP

            endif
          endif
        endif
      endif

!!! END IF

!  WRITE(*,*)
!  WRITE(*,*)
!  WRITE(*,*) ' No Z discretization allowed in GIMRT option at present time'
!  WRITE(*,*)
!  nz = 1
!  nzonez = 0

 dxxt = dxxt/dist_scale
 dyyt = dyyt/dist_scale
 dzzt = dzzt/dist_scale

 dist_scale = 1.0d0

!!! ExportGridLocations  true
ExportGridLocations = .FALSE.
parchar = 'ExportGridLocations'
parfind = ' '
CALL read_logical(nout,lchar,parchar,parfind,ExportGridLocations)

IF (ExportGridLocations) THEN
WRITE(*,*)
WRITE(*,*) 'Exporting grid locations'

OPEN(unit=98,file='GridLocations.txt',status='unknown')

END IF

!*************

ELSE
WRITE(*,*)
WRITE(*,*) ' Failed to find discretization block'
WRITE(*,*)
nx = 1
ny = 1
nz = 1
nzonex = 0
nzoney = 0
nzonez = 0
nxyz = 1
END IF

! end of block to skip for ALQUIMIA
#else
! ALQUIMIA is defined

! alquimia single-cell chemistry
nx = 1
ny = 1
nz = 1
nzonex = 0
nzoney = 0
nzonez = 0
nxyz = 1

#endif
! end of ALQUIMIA block



#ifndef ALQUIMIA

!*****************************************************
!  Write out information on discretization

!      write(iunit2,*)
!      write(iunit2,1016) nzonex
!      write(iunit2,1018)
!      do i = 1,nzonex
!        write(iunit2,1017) nvx(i),dxxt(i)
!      end do
!      write(iunit2,*)

!      write(iunit2,*)
!      write(iunit2,1026) nzoney
!      write(iunit2,1018)
!      do i = 1,nzoney
!        write(iunit2,1017) nvy(i),dyyt(i)
!      end do
!      write(iunit2,*)

!      write(iunit2,*)
!      write(iunit2,1036) nzonez
!      write(iunit2,1018)
!      do i = 1,nzonez
!        write(iunit2,1017) nvz(i),dzzt(i)
!      end do
!      write(iunit2,*)

!      write(*,*)
!      write(*,1016) nzonex
!      write(*,1018)
!      do i = 1,nzonex
!        write(*,1017) nvx(i),dxxt(i)
!      end do
!      write(*,*)

!      write(*,*)
!      write(*,1026) nzoney
!      write(*,1018)
!      do i = 1,nzoney
!        write(*,1017) nvy(i),dyyt(i)
!      end do
!      write(*,*)

!      write(*,*)
!      write(*,1036) nzonez
!      write(*,1018)
!      do i = 1,nzonez
!        write(*,1017) nvz(i),dzzt(i)
!      end do
!      write(*,*)

IF (nzonex == 0) THEN
!         write(*,*)
!         write(*,*) '   No zones in X direction specified'
!         write(*,*) '   Setting NX = 1 '
ELSE
nsum = 0
DO i = 1,nzonex
  nsum = nsum + nvx(i)
END DO
END IF

IF (nzoney == 0) THEN
!        write(*,*)
!        write(*,*) '  No zones in Y direction specified'
!        write(*,*) '  Setting NY = 1'
ELSE
nsum = 0
DO i = 1,nzoney
  nsum = nsum + nvy(i)
END DO
END IF

IF (nzonez == 0) THEN
!        write(*,*)
!        write(*,*) '  No zones in Z direction specified'
!        write(*,*) '  Setting NZ = 1'
ELSE
nsum = 0
DO i = 1,nzonez
  nsum = nsum + nvz(i)
END DO
END IF

!  Now calculate NX, NY, and NZ

IF (nzonex == 0) THEN
nx = 1

IF (ALLOCATED(x)) THEN
  DEALLOCATE(x)
  ALLOCATE(x(0:nx))
ELSE
  ALLOCATE(x(0:nx))
END IF
IF (ALLOCATED(dxx)) THEN
  DEALLOCATE(dxx)
  ALLOCATE(dxx(-1:nx+2))
ELSE
  ALLOCATE(dxx(-1:nx+2))
END IF


dxx(-1) = 1.0d0
dxx(0) = 1.0d0
dxx(1) = 1.0d0
dxx(2) = 1.0d0
dxx(3) = 1.0d0
x(0) = 0.0d0
x(1) = 0.5d0*dxx(1)

ELSE

nx = 0
DO i = 1,nzonex
  IF (nvx(i) == 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Zone with 0 grid cells specified in X dir.'
    READ(*,*)
    STOP
  END IF
  DO jx = 1,nvx(i)
    nx = nx + 1
  END DO
END DO

!  Allocate "X" arrays that depend on NX here

IF (ALLOCATED(x)) THEN
  DEALLOCATE(x)
  ALLOCATE(x(0:nx))
ELSE
  ALLOCATE(x(0:nx))
END IF
IF (ALLOCATED(dxx)) THEN
  DEALLOCATE(dxx)
  ALLOCATE(dxx(-1:nx+2))
ELSE
  ALLOCATE(dxx(-1:nx+2))
END IF

jxx = 0
DO i = 1,nzonex
  DO jx = 1,nvx(i)
    jxx = jxx + 1
    dxx(jxx) = dxxt(i)
  END DO
END DO

dxx(-1) = dxx(1)
dxx(0) = dxx(1)
dxx(nx+1) = dxx(nx)
dxx(nx+2) = dxx(nx)

x(0) = 0.0d0
x(1) = 0.5d0*dxx(1)
DO jx = 2,nx
  x(jx) = x(jx-1) + 0.5d0*(dxx(jx)+dxx(jx-1))
END DO

END IF

1021 FORMAT(2X,'Cell',2X,'Distance (m)')

!      write(*,*) '  ****DISCRETIZATION****'
!      write(*,*)
!      write(*,*) '  NX = ',nx
!      write(*,*)
!      write(*,1021)
!      do jx = 1,nx
!        write(*,1019) jx,x(jx)
!      end do

!      write(*,*)

!      write(iunit2,*) '  ****DISCRETIZATION****'
!      write(iunit2,*)
!      write(iunit2,*) '  NX = ',nx
!      write(iunit2,*)
!      write(iunit2,1021)
!      do jx = 1,nx
!        write(iunit2,1019) jx,x(jx)
!      end do

!      write(iunit2,*)


1019 FORMAT(1X,i3,1X,1PE12.3)

IF (nzoney == 0) THEN
ny = 1

IF (ALLOCATED(y)) THEN
  DEALLOCATE(y)
  ALLOCATE(y(0:ny))
ELSE
  ALLOCATE(y(0:ny))
END IF
IF (ALLOCATED(dyy)) THEN
  DEALLOCATE(dyy)
  ALLOCATE(dyy(-1:ny+2))
ELSE
  ALLOCATE(dyy(-1:ny+2))
END IF


dyy(-1) = 1.0d0
dyy(0) = 1.0d0
dyy(1) = 1.0d0
dyy(2) = 1.0d0
dyy(3) = 1.0d0
y(0) = 0.0d0
y(1) = 0.5d0*dyy(1)
ELSE

ny = 0
DO i = 1,nzoney
  IF (nvy(i) == 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Zone with 0 grid cells specified in Y dir.'
    READ(*,*)
    STOP
  END IF
  DO jy = 1,nvy(i)
    ny = ny + 1
  END DO
END DO

!  Allocate "Y" arrays that depend on NY here

IF (ALLOCATED(y)) THEN
  DEALLOCATE(y)
  ALLOCATE(y(0:ny))
ELSE
  ALLOCATE(y(0:ny))
END IF
IF (ALLOCATED(dyy)) THEN
  DEALLOCATE(dyy)
  ALLOCATE(dyy(-1:ny+2))
ELSE
  ALLOCATE(dyy(-1:ny+2))
END IF

jyy = 0
DO i = 1,nzoney
  DO jy = 1,nvy(i)
    jyy = jyy + 1
    dyy(jyy) = dyyt(i)
  END DO
END DO

dyy(-1) = dyy(1)
dyy(0) = dyy(1)
dyy(ny+1) = dyy(ny)
dyy(ny+2) = dyy(ny)

y(0) = 0.0d0
y(1) = 0.5d0*dyy(1)
DO jy = 2,ny
  y(jy) = y(jy-1) + 0.5d0*dyy(jy-1) + 0.5d0*dyy(jy)
END DO

END IF

!      write(*,*)
!      write(*,*) '  NY = ',ny
!      write(*,*)
!      write(*,1021)
!      do jy = 1,ny
!        write(*,1019) jy,y(jy)
!      end do

!      write(*,*)

!      write(iunit2,*)
!      write(iunit2,*) '  NY = ',ny
!      write(iunit2,*)
!      write(iunit2,1021)
!      do jy = 1,ny
!        write(iunit2,1019) jy,y(jy)
!      end do

!      write(iunit2,*)

IF (nzonez == 0) THEN
nz = 1

IF (ALLOCATED(z)) THEN
  DEALLOCATE(z)
  ALLOCATE(z(0:nz))
ELSE
  ALLOCATE(z(0:nz))
END IF
IF (ALLOCATED(dzz)) THEN
  DEALLOCATE(dzz)
  ALLOCATE(dzz(0:nx+1,0:ny+1,-1:nz+2))
ELSE
  ALLOCATE(dzz(0:nx+1,0:ny+1,-1:nz+2))
END IF

dzz = 1.0d0
z(0) = 0.0d0
z(1) = 0.5d0*dzz(1,1,1)
ELSE

nz = 0
DO i = 1,nzonez
  IF (nvz(i) == 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Zone with 0 grid cells specified in Z dir.'
    READ(*,*)
    STOP
  END IF
  DO jz = 1,nvz(i)
    nz = nz + 1
  END DO
END DO

!  Allocate "Z" arrays that depend on NZ here

IF (ALLOCATED(z)) THEN
  DEALLOCATE(z)
  ALLOCATE(z(0:nz))
ELSE
  ALLOCATE(z(0:nz))
END IF
IF (ALLOCATED(dzz)) THEN
  DEALLOCATE(dzz)
  ALLOCATE(dzz(0:nx+1,0:ny+1,-1:nz+2))
ELSE
  ALLOCATE(dzz(0:nx+1,0:ny+1,-1:nz+2))
END IF

jzz = 0
DO i = 1,nzonez
  DO jz = 1,nvz(i)
    jzz = jzz + 1
    dzz(:,:,jzz) = dzzt(i)
  END DO
END DO

dzz(:,:,-1) = dzz(:,:,1)
dzz(:,:,0) = dzz(:,:,1)
dzz(:,:,nz+1) = dzz(:,:,nz)
dzz(:,:,nz+2) = dzz(:,:,nz)

z(0) = 0.0d0
z(1) = 0.5*dzz(1,1,1)
DO jz = 2,nz
  z(jz) = z(jz-1) + 0.5*(dzz(1,1,jz)+dzz(1,1,jz-1))
END DO

END IF

WRITE(*,*)
WRITE(*,*) '  NZ = ',nz
WRITE(*,*)
!!WRITE(*,1021)
!!DO jz = 1,nz
!!  WRITE(*,1019) jz,z(jz)
!!END DO

!!WRITE(*,*)

WRITE(iunit2,*)
WRITE(iunit2,*) '  NZ = ',nz
WRITE(iunit2,*)
WRITE(iunit2,1021)
DO jz = 1,nz
WRITE(iunit2,1019) jz,z(jz)
END DO
WRITE(iunit2,*)
#endif
nxyz = nx*ny*nz

IF (streamtube) THEN

!!!  Now check to see if Volumes will be read for 1D streamtube case (if so, skip Y and Z discretization)
ReadGridVolumes = .FALSE.
GridVolumeFileFormat = 'FullForm'
parchar = 'read_volumes'
parfind = ' '
GridVolumeFile = ' '
CALL readFileName(nout,lchar,parchar,parfind,dumstring,section,GridVolumeFileFormat)
IF (parfind == ' ') THEN  !
  GridVolumeFile = 'streamtube.vol'            ! Use default
  GridVolumeFileFormat = 'FullForm'
ELSE
  GridVolumeFile = dumstring
  CALL stringlen(GridVolumeFile,ls)
  WRITE(*,*) ' Reading grid volumes from file: ',GridVolumeFile(1:ls)
  ReadGridVolumes = .TRUE.
END IF

!!!  NOTE:  Assumes only 1D file, since this is for a 1D streamtube

IF (GridVolumeFile /= ' ') THEN
ALLOCATE(work1(nx))
INQUIRE(FILE=GridVolumeFile,EXIST=ext)

IF (.NOT. ext) THEN
  CALL stringlen(GridVolumeFile,ls)
  WRITE(*,*)
  WRITE(*,*) ' GridVolumeFile file not found: ', GridVolumeFile(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

OPEN(UNIT=52,FILE=GridVolumeFile,STATUS='OLD',ERR=6001)
FileTemp = GridVolumeFile
CALL stringlen(FileTemp,FileNameLength)

IF (GridVolumeFileFormat == 'ContinuousRead') THEN

  READ(52,*,END=1020) (work1(jx),jx=1,nx)

ELSE IF (GridVolumeFileFormat == 'SingleColumn') THEN

      DO jx= 1,nx
        READ(52,*,END=1020) work1(jx)
      END DO

ELSE IF (GridVolumeFileFormat == 'FullForm') THEN

    DO jx= 1,nx
      READ(52,*,END=1020) xdum,work1(jx)
    END DO

ELSE IF (GridVolumeFileFormat == 'Unformatted') THEN

  READ(52,END=1020) work1

ELSE

!!! If no file format is given, assume FullForm (with dummy X)

    DO jx= 1,nx
      READ(52,*,END=1020) xdum,work1(jx)
    END DO

END IF
CLOSE(UNIT=52)

DO jx = 1,nx
  dzz(jx,1,1) = work1(jx)
END DO

!!! Make ghost cells the same as first grid cell within domain
dzz(0,1,1) = dzz(1,1,1)
dzz(nx+1,1,1) = dzz(nx,1,1)

END IF
IF (ExportGridLocations) THEN

!!!  do jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      write(98,*) x(jx)*1.0D+6 ,y(jy)*1.0D+6
    END DO
  END DO
!!!  END DO

close(unit=98,status='keep')

END IF
END IF

!!!  End of DISCRETIZATION

IF (nxyz > nx .OR. nxyz == 1) THEN
IF (spherical) THEN
  WRITE(*,*)
  WRITE(*,*) ' Spherical coordinates allowed presently only for 1D (in X) case'
  WRITE(*,*) ' Aborting run'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
END IF

IF (nxyz == 1) THEN
ihindmarsh = 0
petscon = .FALSE.
END IF

IF (nx > 1 .AND. ny > 1 .AND. gimrt) THEN      !  2D or greater problem
IF (SolveHindmarsh .OR. ihindmarsh == 1) THEN
  WRITE(*,*)
  WRITE(*,*) ' Hindmarsh block tridiagonal solver cannot be used for 2D/3D problems'
  WRITE(*,*) ' Switching to PETSc '
  WRITE(*,*) ' Return to continue'
  READ(*,*)
  ihindmarsh = 0
  petscon = .TRUE.
END IF
END IF

IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN               !  3D problem, use OS3D
IF (gimrt) THEN
  WRITE(*,*)
  WRITE(*,*) ' 3D not currently available with the global implicit (gimrt) option'
  WRITE(*,*) ' Hit "RETURN" to proceed with the OS3D option'
  READ(*,*)
  os3d = .TRUE.
END IF
IF (os3d) THEN
  gimrt = .FALSE.
  petscon = .FALSE.
END IF
END IF

CALL REALLOCATE(ncomp,nspec,nrct,nkin,ngas,nsurf,nexchange,ikin,nexch_sec,nsurf_sec)
! biomass
call read_CatabolicPath(ncomp,nkin,ikin,data3)
! biomass end


IF (nradmax > 0) THEN
ALLOCATE(mumin_decay(ndim1,ndim2,ndim3,nx,ny,nz))
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      mumin_decay(:,:,:,jx,jy,jz) = mumin(:,:,:)
    END DO
  END DO
END DO
END IF

!  **************  ReALLOCATE 4D arrays  *****************************

!  *******  Real arrays  ****************

!ndim1 = nisotope_max
!ndim2 = nmindecay_max
!ndim3 = ndecay
!ndim4 = nchem
!i = size(ratio_isotope_init,1)
!j = size(ratio_isotope_init,2)
!k = size(ratio_isotope_init,3)
!l = size(ratio_isotope_init,4)
!ALLOCATE(work4(i,j,k,l))
!work4 = ratio_isotope_init
!DEALLOCATE(ratio_isotope_init)
!ALLOCATE(ratio_isotope_init(ndim1,ndim2,ndim3,ndim4))
!IF(ndim1 /= 0.AND.ndim2 /= 0.AND.ndim3 /= 0.AND.ndim4 /= 0) THEN
!  ratio_isotope_init(1:ndim1,1:ndim2,1:ndim3,1:ndim4) = work4(1:ndim1,1:ndim2,1:ndim3,1:ndim4)
!END IF
!DEALLOCATE(work4)

!   *****************  GLOBAL ARRAYS  ***********************
!    Allocate global arrays, mostly over the spatial domain
CALL GlobalArrayAllocation(ncomp,nspec,nkin,nrct,ngas,npot,nexchange,nexch_sec,nsurf,nsurf_sec,ikin,nx,ny,nz)
qx = 0.0d0
qy = 0.0d0
qz = 0.0d0
! Zhi Li

IF (ALLOCATED(us)) THEN
DEALLOCATE(us)
ALLOCATE(us(nx+1,ny,nz))
ELSE
ALLOCATE(us(nx+1,ny,nz))
END IF

IF (ALLOCATED(vs)) THEN
DEALLOCATE(vs)
ALLOCATE(vs(nx+1,ny,nz))
ELSE
ALLOCATE(vs(nx+1,ny,nz))
END IF
dspx = 0.0d0
dspy = 0.0d0
dspz = 0.0d0
qxgas = 0.0d0
qygas = 0.0d0
qzgas = 0.0d0
netflowx = 0.0d0
netDiffuseX = 0.0d0
satliq = FixSaturation
satliqold = FixSaturation
dstar = 0.0d0
por = 1.0d0
porin = 1.0d0
dxy = 1.0d0
exchangesites = 0
volfx = 0.0d0
VolumeLastTimeStep = 1.0E-10
dppt = 0.0d0
area = 0.0d0
xgram = 1.0d0

s = 0.0d0
sn = 0.0d0

si = 1.0d0
silog = 0.0d0
actenergy = 0.0d0
pre_rmin = 0.0d0
rmin = 0.0d0
jac_rmin = 0.0d0
d_sp = 0.0d0
pre_raq = 0.0d0
raq_tot = 0.0d0
rdkin = 0.0d0
raq = 0.0d0

!  These are deallocated at the end of START98
ALLOCATE(ndist(nxyz))
ALLOCATE(jxxlo(nxyz))
ALLOCATE(jxxhi(nxyz))
ALLOCATE(jyylo(nxyz))
ALLOCATE(jyyhi(nxyz))
ALLOCATE(jzzlo(nxyz))
ALLOCATE(jzzhi(nxyz))
ALLOCATE(jjfix(nxyz))
jxxlo = 1
jxxhi = 1
jyylo = 1
jyyhi = 1
jzzhi = 1
jzzlo = 1
jjfix = 0
t = tinit

! *****************************************************************

! skip what follows if alquimia is defined
#ifndef ALQUIMIA
!    ***************INTERNAL HETEROGENEITIES********************
section = 'initial_conditions'
CALL readblock(nin,nout,section,found,ncount)
IF (found) THEN
!!  WRITE(*,*)
!!  WRITE(*,*) ' Initial conditions block found'
!!  WRITE(*,*)
ELSE
WRITE(*,*)
WRITE(*,*) ' No initial conditions found'
WRITE(*,*) ' --->  No defaults available: Aborting run'
WRITE(*,*)
READ(*,*)
STOP
END IF

CALL read_het(nout,nchem,nhet,nx,ny,nz)

IF (ReadInitialConditions .and. InitialConditionsFile /= ' ') THEN

ALLOCATE(work3(nx,ny,nz))
INQUIRE(FILE=InitialConditionsFile,EXIST=ext)
IF (.NOT. ext) THEN
  CALL stringlen(InitialConditionsFile,ls)
  WRITE(*,*)
  WRITE(*,*) ' InitialConditionsFile not found: ', InitialConditionsFile(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

OPEN(UNIT=52,FILE=InitialConditionsFile,STATUS='OLD',ERR=6001)
FileTemp = InitialConditionsFile
CALL stringlen(FileTemp,FileNameLength)

IF (MontTerri) THEN

  nhet = 0
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx= 1,nx
        
        nhet = nhet + 1
        READ(52,*,END=1020) xdum,ydum,zdum, work3(jx,jy,jz)
!!!                            x    y    z    condition #             
        jinit(jx,jy,jz) = work3(jx,jy,jz) 
        activecell(jx,jy,jz) = 1

      END DO
    END DO
  END DO

  CLOSE(UNIT=52)

END IF

if (nmmLogical) then

  jz = 1
  ALLOCATE(stress(nx,ny,1))
  
  OPEN(UNIT=53,FILE='RoughFracture-Dec11.dat',STATUS='OLD',ERR=6001)
  allocate(work3b(nx,ny,nz))
  
  jinit = 2

  nhet = 0
  DO jy = 1,ny
    DO jx= 1,nx
      
      nhet = nhet + 1
      IF (SaltCreep) THEN
        READ(52,*,END=1020) xdum,ydum,zdum, work3(jx,jy,jz), xdum, ydum, zdum, xdum, ydum, xdum, stress(jx,jy,jz), zdum,   xdum
!!!                            x    y    bn    mt             sx    sy    txy   dx    dy    sig1  sig3            re-sig1  re-sig

        jinit(jx,jy,jz) = DNINT(work3(jx,jy,jz)) + 1
        
!!! --- > CalciteCreep *********************************
        
      ELSE IF (CalciteCreep) THEN
        
        READ(52,*,END=1020) xdum,ydum,zdum, work3(jx,jy,jz), xdum, ydum, zdum, xdum, ydum, xdum, stress(jx,jy,jz), zdum,   xdum
!!!                            x    y    bn    mt             sx    sy    txy   dx    dy    sig1  sig3            re-sig1  re-sig

        
!!!        READ(53,*,END=1020) xdum,ydum, work3(jx,jy,jz), zdum
!!!                            x    y    bn    mt   
        
        work3(jx,jy,jz) = work3(jx,jy,jz) + 1.0
        
        if (work3(jx,jy,jz) > 0.999d0 .and. work3(jx,jy,jz) < 1.001d0) then
          jinit(jx,jy,jz) = 1   
        end if
        
        READ(53,*,END=1020) xdum,ydum, work3b(jx,jy,jz), zdum
!!!                                  x    y       mt            junk
        
        if (work3b(jx,jy,jz) > 1.999d0 .and. work3b(jx,jy,jz) <2.001) then
          if (jinit(jx,jy,jz) == 1) then
            continue
          else
            jinit(jx,jy,jz) = 3 
          end if
        end if 
        
        if (jy > 120) then
          jinit(jx,jy,jz) = 2
!!!          stress(jx,jy,jz) = 2.50E06
        end if
        if (jy < 30) then
          jinit(jx,jy,jz) = 2
!!!          stress(jx,jy,jz) = 2.50E06
        end if
      
        
      ELSE IF (FractureNetwork) THEN

        READ(52,*,END=1020) xdum,ydum,zdum, work3(jx,jy,jz)
!!!                            x    y    bn    mt

        jinit(jx,jy,jz) = DNINT(work3(jx,jy,jz))

      ELSE
        CONTINUE
      ENDIF
      activecell(jx,jy,jz) = 1
    END DO
  END DO

  CLOSE(UNIT=52)
  CLOSE(UNIT=53)
  deallocate(work3b)
  
  stress = stress/2.0
  !!!stress = 10**5
  
  StressMaxVal= MaxVal(ABS(stress*1.0E-5))

  
  write(*,*)
  write(*,*) ' StressMaxVal (bars) =', StressMaxVal
  write(*,*)


END IF

END IF

IF (nhet == 0) THEN
WRITE(*,*)
WRITE(*,*) ' No initial conditions found'
WRITE(*,*) ' --->  No defaults available: Aborting run'
WRITE(*,*)
READ(*,*)
STOP
END IF


IF (nhet > 0 .and. .not. ReadInitialConditions) THEN
DO l = 1,nhet
  IF (jxxhi(l) > nx) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a corner at JX > NX'
    READ(*,*)
    STOP
  END IF
  IF (jyyhi(l) > ny) THEN
    write(*,*) ' Ny = ', ny
    write(*,*) ' jyyhi = ',jyyhi(l)
    WRITE(*,*)
    WRITE(*,*) 'You have specified a corner at JY > NY'
    READ(*,*)
    STOP
  END IF
  IF (jzzhi(l) > nz) THEN
    WRITE(*,*)
    WRITE(*,*) 'You have specified a corner at JZ > NZ'
    READ(*,*)
    STOP
  END IF
END DO

!  Save geochemical condition numbers in JINIT pointer array

jinit = 0
activecell = 1

DO l = 1,nhet
  ll = l + 1
  DO jz = jzzlo(l),jzzhi(l)
    DO jy = jyylo(l),jyyhi(l)
      DO jx = jxxlo(l),jxxhi(l)
        jinit(jx,jy,jz) = ndist(l)
        IF (jjfix(l) == 1) THEN
          activecell(jx,jy,jz) = 0
        ELSE
          activecell(jx,jy,jz) = 1
        END IF
      END DO
    END DO
  END DO
END DO

!! Or overwrite with read of geochemical conditions
IF (ReadGeochemicalConditions) THEN

  INQUIRE(FILE='InvadingCluster.dat',EXIST=ext)

  IF (ext) THEN

    OPEN(unit=98,file='InvadingCluster.dat',status='old')

    read(98,*) GridCoordinateX,GridCoordinateY
    do jy = 1,ny
      do jx = 1,nx
        if (jx == GridCoordinateX .AND. jy == GridCoordinateY) then
          jinit(jx,jy,1) = 2
          read(98,*,END=980) GridCoordinateX,GridCoordinateY
        else
          jinit(jx,jy,1) = 1
        endif
      end do
    end do

    close(unit=98,status='keep')
980   write(*,*) ' End of InvadingCluster.dat file'

  ELSE
    continue
  END IF

    END IF    !! End of ReadGeochemicalConditions = .TRUE.

  IF (ReadGautier) THEN

    IF (ALLOCATED(tortuosity)) THEN
      DEALLOCATE(tortuosity)
      ALLOCATE(tortuosity(nx,ny,nz))
    ELSE
      ALLOCATE(tortuosity(nx,ny,nz))
    END IF

    tortuosity = 0.1d0

  INQUIRE(FILE='FILE15.dat',EXIST=ext)

!!%% 1st column = averaged pixel number (65536 averaged pixels = 65536 rows)
!!%% 2nd column = most abundant mineral in the averaged pixel.
!!Number ID: 1= PORE/EMPTINESS ; 2= QUARTZ ; 3= CHLORITE ; 4=ILLITE/MICA ; 5=TI OXIDE; 6=KAOLINITE ; 7=ILLITE/SMECTITE ; 8=FE OXIDE
!!%% 3rd column = PORE percentage in the averaged pixel
!!%% 4th column = QUARTZ percentage in the averaged pixel
!!%% 5th column = CHLORITE percentage in the averaged pixel
!!%% 6th column = ILLITE/MICA percentage in the averaged pixel
!!%% 7th column = TI OXIDE percentage in the averaged pixel
!!%% 8th column = KAOLINITE percentage in the averaged pixel
!!%% 9th column = ILLITE/SMECTITE percentage in the averaged pixel
!!%% 10th column = FE OXIDE percentage in the averaged pixel

!! Quartz = 1
!! Chlorite = 2
!! Illite = 3
!! Kaolinite = 4
!! Smectite = 5
!! FeOxide = 6

  IF (ext) THEN

    OPEN(unit=98,file='FILE15.dat',status='old')

    jz = 1
    do jy = 1,ny
      do jx = 1,nx

         IF (PorosityRead == 0.00) THEN
           por(jx,jy,1) = 0.01d0
         ELSE
           por(jx,jy,1) = PorosityRead/100.0d0
         END IF
         IF (QuartzRead == 100.0) THEN
           volfx(1,jx,jy,1) = 0.99
         ELSE
           volfx(1,jx,jy,1) = QuartzRead/100.0d0
         END IF
         IF (ChloriteRead == 100.0) THEN
           volfx(2,jx,jy,1) = 0.99
         ELSE
           volfx(2,jx,jy,1) = ChloriteRead/100.0d0
           IF (volfx(2,jx,jy,1) > 0.50) THEN
!! Use estimate from FIB-SEM modeling
               tortuosity(jx,jy,jz) = 0.007
           END IF
         END IF

         IF (IlliteRead == 100.0) THEN
           volfx(3,jx,jy,1) = 0.99
         ELSE
           volfx(3,jx,jy,1) = IlliteRead/100.0d0
         END IF
         IF (KaoliniteRead == 100.0) THEN
           volfx(4,jx,jy,1) = 0.99
         ELSE
           volfx(4,jx,jy,1) = KaoliniteRead/100.0d0
         END IF
         IF (SmectiteRead == 100.0) THEN
           volfx(5,jx,jy,1) = 0.99
         ELSE
           volfx(5,jx,jy,1) = SmectiteRead/100.0d0
         END IF
         IF (FeOxideRead == 100.0) THEN
           volfx(6,jx,jy,1) = 0.99
         ELSE
           volfx(6,jx,jy,1) = FeOxideRead/100.0d0
         END IF
         IF (por(jx,jy,jz) < 0.05) THEN
           tortuosity(jx,jy,jz) = 0.001
         ELSE
           tortuosity(jx,jy,jz) = 0.1
         END IF

      end do
    end do


    close(unit=98,status='keep')

  ELSE
    continue
  END IF

    END IF    !! End of ReadGautier = .TRUE.


ELSE
IF (ReadInitialConditions) THEN
  WRITE(*,*)
  WRITE(*,*) ' Initial conditions read from file'
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Initial conditions must be specified'
  WRITE(*,*) ' No DEFAULT condition assumed'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
END IF


5005 FORMAT(' No initial condition specified at grid cell: ',i4,1x,i4,1x,i4)

DO jz = 1,nz
DO jy = 1,ny
  DO jx = 1,nx
    IF (jinit(jx,jy,jz) == 0) THEN
      WRITE(*,*)
      WRITE(*,5005) jx,jy,jz
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
  END DO
END DO
END DO

IF (jtemp == 2) THEN

INQUIRE(FILE=TFile,EXIST=ext)

IF (.NOT. ext) THEN
  CALL stringlen(TFile,ls)
  WRITE(*,*)
  WRITE(*,*) ' Temperature file for read not found: ', TFile(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

OPEN(UNIT=52,FILE=TFile,STATUS='OLD',ERR=6002)
FileTemp = TFile
CALL stringlen(FileTemp,FileNameLength)
IF (TemperatureFileFormat == 'ContinuousRead') THEN
  READ(52,*,END=1020) (((t(jx,jy,jz),jx=1,nx),jy=1,ny),jz=1,nz)
ELSE IF (TemperatureFileFormat == 'SingleColumn') THEN
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx= 1,nx
        READ(52,*,END=1020) t(jx,jy,jz)
      END DO
    END DO
  END DO
ELSE IF (TemperatureFileFormat == 'FullForm') THEN
  IF (ny > 1 .AND. nz > 1) THEN
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx= 1,nx
          READ(52,*,END=1020) xdum,ydum,zdum,t(jx,jy,jz)
        END DO
      END DO
    END DO
  ELSE IF (ny > 1 .AND. nz == 1) THEN
    jz = 1
    jy = 1
    DO jy = 1,ny
      DO jx= 1,nx
        READ(52,*,END=1020) xdum,ydum,t(jx,jy,jz)
      END DO
    END DO
  ELSE
    jz = 1
    jy = 1
    DO jx= 1,nx
      READ(52,*,END=1020) xdum,t(jx,jy,jz)
    END DO
  END IF
ELSE IF (TemperatureFileFormat == 'Unformatted') THEN
  READ(52,END=1020) t
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Temperature file format not recognized'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

CLOSE(UNIT=52)

ELSE IF (jtemp == 1) THEN                         !! Temperature gradient in X direction specified

WRITE(*,*)
WRITE(*,*) '  A temperature gradient has been specified'
WRITE(*,*) '  Gradient operating ONLY in X direction'
WRITE(*,*)
IF (tgrad == 0.0) THEN
  WRITE(*,*)
  WRITE(*,*) 'You have specified a non-isothermal run '
  WRITE(*,*) 'without a non-zero temperature gradient'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
        t(jx,jy,jz) = tinit + tgrad*x(jx)
    END DO
  END DO
END DO

! ************************************
! Edit by Lucien Stolze, June 2023
! Temperature time series
ELSEIF (jtemp == 3) THEN !! Temperature time series allocated to specific regions

  CALL read_tempregion(nout,nx,ny,nz,len(TFile),TFile,TemperatureFileFormat)

  IF (RunTempts) THEN
    t = t_default
    DO jz = 1,nz
      DO jy = 1,ny
      DO jx = 1,nx
      DO i = 1,nb_temp_ts
          IF (temp_region(jx,jy,jz) == reg_temp_ts(i)) THEN
            t(jx,jy,jz) = temp_ts(i,1)
            ! IF (jx == 1) THEN
            ! t(0,jy,jz) = t(jx,jy,jz)
            ! ENDIF
            ! IF (jx == nx) THEN
            ! t(nx+1,jy,jz) = t(jx,jy,jz)
            ! ENDIF
            ! IF (jy == 1) THEN
            ! t(jx,0,jz) = t(jx,jy,jz)
            ! ENDIF
            ! IF (jy == ny) THEN
            ! t(jx,ny+1,jz) = t(jx,jy,jz)
            ! ENDIF
            ! IF (jz == 1) THEN
            ! t(jx,jy,0) = t(jx,jy,jz)
            ! ENDIF
            ! IF (jz == nz) THEN
            ! t(jx,jy,nz+1) = t(jx,jy,jz)
            ! ENDIF
          ENDIF
        END DO
        !Allocate fixed temperature to the regions:
          DO j = 1,nb_temp_fix
            IF (temp_region(jx,jy,jz) == reg_temp_fix(j)) THEN
              t(jx,jy,jz) = temp_fix(j)
              ! IF (jx == 1) THEN
              ! t(0,jy,jz) = t(jx,jy,jz)
              ! ENDIF
              ! IF (jx == nx) THEN
              ! t(nx+1,jy,jz) = t(jx,jy,jz)
              ! ENDIF
              ! IF (jy == 1) THEN
              ! t(jx,0,jz) = t(jx,jy,jz)
              ! ENDIF
              ! IF (jy == ny) THEN
              ! t(jx,ny+1,jz) = t(jx,jy,jz)
              ! ENDIF
              ! IF (jz == 1) THEN
              ! t(jx,jy,0) = t(jx,jy,jz)
              ! ENDIF
              ! IF (jz == nz) THEN
              ! t(jx,jy,nz+1) = t(jx,jy,jz)
              ! ENDIF
            ENDIF
        END DO
      
      END DO
      END DO
      END DO
  
  ENDIF
! ************************************
! Finish edit by Lucien Stolze, June 2023

ELSEIF (jtemp == 0) THEN
  t = tinit
ELSE

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      t(jx,jy,jz) = tempcond(jinit(jx,jy,jz))
    END DO
  END DO
END DO

END IF

DO jz = 1,nz
DO jy = 1,ny
  DO jx = 1,nx
    CALL density(jx,jy,jz)
  END DO
END DO
END DO

!*************************************************************
! Edit by Toshiyuki Bandai, June 2023
! Edit by Lucien Stolze, June 2023
! compute the dynamics viscosity of water [Pa year] based on local temperature
! mu_water is defined in Flow module
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      mu_water(jx,jy,jz) = 10.0d0**(-4.5318d0 - 220.57d0/(149.39 - t(jx,jy,jz) - 273.15d0)) * 86400.0d0 * 365.0d0 ! 
      rho_water_2 = 0.99823d0 * 1.0E3
      !rho_water2 = 1000.0d0*(1.0d0 - (t(jx,jy,jz) + 288.9414d0) / (508929.2d0*(t(jx,jy,jz) + 68.12963d0))*(t(jx,jy,jz)-3.9863d0)**2.0d0)
    END DO
  END DO
END DO
!     END DO
!   END DO
! END DO

mu_water = 0.001*secyr

! End of Edit by Toshiyuki Bandai, June 2023
!*************************************************************

roOld = ro

!! WRITE(iunit2,*)
!! WRITE(iunit2,*) '  TEMPERATURE FIELD'
!! WRITE(iunit2,*)
!! WRITE(iunit2,*) ' Distance (m)   T (C)'
!! jy = 1
!! jz = 1
!! DO jx = 1,nx
!fp! set_index({#ident# jy #});
!fp! set_index({#ident# jz #});
!MATT HACK
!!   WRITE(iunit2,889) x(jx)
!!   WRITE(iunit2,889) t(jx,jy,jz)
!! END DO
!! WRITE(iunit2,*)

889   FORMAT(1X,f10.4,1X,f10.2)
WRITE(*,*)
WRITE(*,*) ' Number of heterogeneities = ', nhet
WRITE(*,*)

WRITE(iunit2,*)
WRITE(iunit2,*) 'Number of heterogeneities = ', nhet
WRITE(iunit2,*)
!  ******************* INITIALIZATION OF SPATIAL DOMAIN ***************************
sexch = 0.0
fexch = 0.0

IF (SaturationFile /= ' ') THEN
INQUIRE(FILE=SaturationFile,EXIST=ext)
IF (.NOT. ext) THEN
  CALL stringlen(SaturationFile,ls)
  WRITE(*,*)
  WRITE(*,*) ' Liquid saturation file not found: ', SaturationFile(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
OPEN(UNIT=52,FILE=SaturationFile,STATUS='OLD',ERR=6002)
FileTemp = SaturationFile
CALL stringlen(FileTemp,FileNameLength)
IF (SaturationFileFormat == 'ContinuousRead') THEN
  READ(52,*,END=1020) (((satliq(jx,jy,jz),jx=1,nx),jy=1,ny),jz=1,nz)
ELSE IF (SaturationFileFormat == 'SingleColumn') THEN
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx= 1,nx
        READ(52,*,END=1020) satliq(jx,jy,jz)
      END DO
    END DO
  END DO
ELSE IF (SaturationFileFormat == 'FullForm') THEN
  IF (ny > 1 .AND. nz > 1) THEN
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx= 1,nx
          READ(52,*,END=1020) xdum,ydum,zdum,satliq(jx,jy,jz)
        END DO
      END DO
    END DO
  ELSE IF (ny > 1 .AND. nz == 1) THEN
    jz = 1
    jy = 1
    DO jy = 1,ny
      DO jx= 1,nx
        READ(52,*,END=1020) xdum,ydum,satliq(jx,jy,jz)
      END DO
    END DO
  ELSE
    jz = 1
    jy = 1
    DO jx= 1,nx
      READ(52,*,END=1020) xdum,satliq(jx,jy,jz)
    END DO
  END IF
ELSE IF (SaturationFileFormat == 'Unformatted') THEN
  READ(52,END=1020) satliq
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Saturation file format not recognized'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
satliqold = satliq
MinSaturation = MINVAL(satliq(1:nx,1:ny,1:nz))
IF (MinSaturation < 1.0) THEN
  isaturate = 1
  WRITE(*,*)
  WRITE(*,*) ' Running as an unsaturated problem'
  WRITE(*,*)
END IF
CLOSE(UNIT=52)
END IF
satliq(0,1,1) = satliq(1,1,1)

IF (jpor /= 0 .AND. PorosityFile /= ' ') THEN
ALLOCATE(work3(nx,ny,nz))
INQUIRE(FILE=PorosityFile,EXIST=ext)
IF (.NOT. ext) THEN
  CALL stringlen(PorosityFile,ls)
  WRITE(*,*)
  WRITE(*,*) ' Porosity file not found: ', PorosityFile(1:ls)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
OPEN(UNIT=52,FILE=PorosityFile,STATUS='OLD',ERR=6001)
FileTemp = PorosityFile
CALL stringlen(FileTemp,FileNameLength)
IF (PorosityFileFormat == 'ContinuousRead') THEN
  READ(52,*,END=1020) (((work3(jx,jy,jz),jx=1,nx),jy=1,ny),jz=1,nz)
ELSE IF (PorosityFileFormat == 'SingleColumn') THEN
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx= 1,nx
        READ(52,*,END=1020) work3(jx,jy,jz)
      END DO
    END DO
  END DO
ELSE IF (PorosityFileFormat == 'FullForm') THEN
  IF (ny > 1 .AND. nz > 1) THEN
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx= 1,nx
          READ(52,*,END=1020) xdum,ydum,zdum,work3(jx,jy,jz)
        END DO
      END DO
    END DO
  ELSE IF (ny > 1 .AND. nz == 1) THEN
    jz = 1
    DO jy = 1,ny
      DO jx= 1,nx
        READ(52,*,END=1020) xdum,ydum,work3(jx,jy,jz)
      END DO
    END DO
  ELSE
    jz = 1
    jy = 1
    DO jx= 1,nx
      READ(52,*,END=1020) xdum,work3(jx,jy,jz)
    END DO
  END IF
ELSE IF (PorosityFileFormat == 'Unformatted') THEN
  READ(52,END=1020) work3
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Porosity file format not recognized'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
CLOSE(UNIT=52)
END IF

DO jz = 1,nz
DO jy = 1,ny
  DO jx = 1,nx
    j = (jz-1)*nx*ny + (jy-1)*nx + jx
    DO ik = 1,ncomp+nspec
      sp10(ik,jx,jy,jz) = spcond10(ik,jinit(jx,jy,jz))
      sp(ik,jx,jy,jz) = spcond(ik,jinit(jx,jy,jz))
    END DO
    DO i = 1,ncomp
      s(i,jx,jy,jz) = scond(i,jinit(jx,jy,jz))
      sn(i,jx,jy,jz) = scond(i,jinit(jx,jy,jz))
    END DO

    DO kk = 1,ngas
      spgas10(kk,jx,jy,jz) = spcondgas10(kk,jinit(jx,jy,jz))
      spgas(kk,jx,jy,jz) = spcondgas(kk,jinit(jx,jy,jz))
    END DO

    sum = 0.0
    DO k = 1,nrct
      IF (.NOT. ReadGautier) THEN
        volfx(k,jx,jy,jz) = volin(k,jinit(jx,jy,jz))
      END IF
      VolumeLastTimeStep(k,jx,jy,jz) = volfx(k,jx,jy,jz)
      area(k,jx,jy,jz) = areain(k,jinit(jx,jy,jz))
      sum = sum + volfx(k,jx,jy,jz)
    END DO

    IF (Duan .OR. Duan2006) THEN
!!      Save Residual Volume for CO2 fugacity calculation
      vrSave(jx,jy,jz) = vrInitial(jinit(jx,jy,jz))
    END IF

    IF (constantpor /= 0.0 .AND. .NOT. ReadGautier) THEN
       porin(jx,jy,jz) = constantpor
       por(jx,jy,jz) = constantpor
       porOld(jx,jy,jz) = por(jx,jy,jz)
    END IF


    IF (jpor == 2 .OR. jpor == 3) THEN                  !! Read porosity from file
      porin(jx,jy,jz) = work3(jx,jy,jz)
    ELSE
      porin(jx,jy,jz) = porcond(jinit(jx,jy,jz))
    END IF

!! Porosity calculated from aperture crashes in Hang version

    IF (SaturationFile == ' ') THEN
      satliq(jx,jy,jz) = SaturationCond(jinit(jx,jy,jz))
      satliqold(jx,jy,jz) = satliq(jx,jy,jz)
    END IF

    IF (jpor == 3) THEN    ! Renormalize volume fractions for case of porosity read from file and porosity update

      ScaleMineralVolumes = ( (1.0d0-porin(jx,jy,jz)) / (1.0d0-porcond(jinit(jx,jy,jz))) )
!!!        ScaleMineralVolumes = porcond(jinit(jx,jy,jz))/porin(jx,jy,jz)

      SumMineralVolume = 0.0d0
      DO k = 1,nrct
!!!            volin(k,jinit(jx,jy,jz)) = volin(k,jinit(jx,jy,jz)) * ScaleMineralVolumes
!!!            volfx(k,jx,jy,jz) = volin(k,jinit(jx,jy,jz))
          volfx(k,jx,jy,jz) = volin(k,jinit(jx,jy,jz)) * ScaleMineralVolumes
          SumMineralVolume = SumMineralVolume + volfx(k,jx,jy,jz)
          VolumeLastTimeStep(k,jx,jy,jz) = volfx(k,jx,jy,jz)
          area(k,jx,jy,jz) = areain(k,jinit(jx,jy,jz))
      END DO
      porin(jx,jy,jz) = 1.0d0-SumMineralVolume

    END IF

    IF (.NOT. ReadGautier) THEN
      por(jx,jy,jz) = porin(jx,jy,jz)
    END IF
    porOld(jx,jy,jz) = por(jx,jy,jz)
    IF (por(jx,jy,jz) <= 0.0) THEN
      WRITE(*,*)
      WRITE(*,*) '  You have specified a porosity < 0'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

    CALL density(jx,jy,jz)

!  Change the site concentrations to sites per bulk volume porous medium

    IF (DensityModule /= 'temperature') THEN
!       Calculate the correction for the mass fraction of water:  kg_solution/kg_water
      MeanSaltConcentration = 0.001d0*(wtaq(MeanSalt(1))*s(MeanSalt(1),jx,jy,jz) +   &
          wtaq(MeanSalt(2))*s(MeanSalt(2),jx,jy,jz))
      MassFraction = 1.0d0/(1.0d0 + MeanSaltConcentration)
    ELSE
      MassFraction = 1.0d0
    END IF

!!      convert = ro(jx,jy,jz)*porin(jx,jy,jz)*MassFraction
    convert = ro(jx,jy,jz)*porcond(jinit(jx,jy,jz))*SaturationCond(jinit(jx,jy,jz))*MassFraction

    DO ix = 1,nexchange
      exchangesites(ix,jx,jy,jz) = convert*totexch(ix,jinit(jx,jy,jz)) ! Now in equivalents/m3 por. med.
!!        exchangesites(ix,jx,jy,jz) = totexch(ix,jinit(jx,jy,jz)) ! Already in equivalents/m3 por. med.
    END DO

    do ix = 1,nexchange
      spex(ix,jx,jy,jz) = spcondex(ix,jinit(jx,jy,jz))
      spex10(ix,jx,jy,jz) = convert*spcondex10(ix,jinit(jx,jy,jz))  ! Now in eq/m3 por. med.
    end do
    DO ix = 1,nexch_sec
      spex10(ix+nexchange,jx,jy,jz) = convert*spcondex10(ix+nexchange,jinit(jx,jy,jz))  ! Now in eq/m3 por. med.
    END DO

    DO is = 1,nsurf+nsurf_sec
      spsurf10(is,jx,jy,jz) = convert*spcondsurf10(is,jinit(jx,jy,jz))
    END DO
    DO is = 1,nsurf
      spsurf(is,jx,jy,jz) = LOG(convert*spcondsurf10(is,jinit(jx,jy,jz)))
    END DO

  END DO
END DO
END DO

!*****************************
!Stolze Lucien, June 2023
!START: read mineral volume fraction and bulk surface area from file
IF (ALLOCATED(mineral_id)) THEN
DEALLOCATE(mineral_id)
ENDIF
ALLOCATE(mineral_id(50))

IF (ALLOCATED(mineral_name)) THEN
DEALLOCATE(mineral_name)
ENDIF
ALLOCATE(mineral_name(50))

IF (ALLOCATED(mineral_name_length)) THEN
DEALLOCATE(mineral_name_length)
ENDIF
ALLOCATE(mineral_name_length(50))

IF (ALLOCATED(volfracfile)) THEN
DEALLOCATE(volfracfile)
ENDIF
ALLOCATE(volfracfile(50))

IF (ALLOCATED(bsafile)) THEN
DEALLOCATE(bsafile)
ENDIF
ALLOCATE(bsafile(50))

IF (ALLOCATED(lfile_volfrac)) THEN
DEALLOCATE(lfile_volfrac)
ENDIF
ALLOCATE(lfile_volfrac(50))

IF (ALLOCATED(lfile_bsa)) THEN
DEALLOCATE(lfile_bsa)
ENDIF
ALLOCATE(lfile_bsa(50))

IF (ALLOCATED(FileFormatType_volfrac)) THEN
DEALLOCATE(FileFormatType_volfrac)
ENDIF
ALLOCATE(FileFormatType_volfrac(50))

IF (ALLOCATED(FileFormatType_bsa)) THEN
DEALLOCATE(FileFormatType_bsa)
ENDIF
ALLOCATE(FileFormatType_bsa(50))

IF (ALLOCATED(readmin_ssa)) THEN
  DEALLOCATE(readmin_ssa)
  ENDIF
  ALLOCATE(readmin_ssa(50))

mineral_index = 0
CALL read_mineralfile(nout,nx,ny,nz,readmineral,mineral_index,mineral_id,mineral_name,mineral_name_length,volfracfile,bsafile,lfile_volfrac,lfile_bsa,FileFormatType_volfrac,FileFormatType_bsa,readmin_ssa)

IF (readmineral) THEN
DO i = 1,mineral_index
        min_id = mineral_id(i)
        min_name = mineral_name(i)
        min_name_l = mineral_name_length(i)
        vv_file = volfracfile(i)
        vv_file_l = lfile_volfrac(i)
        vv_fileformat = FileFormatType_volfrac(i)
        bsa_file = bsafile(i)
        bsa_file_l = lfile_bsa(i)
        bsa_fileformat = FileFormatType_bsa(i)
        ssa_or_bsa = readmin_ssa(i)

        INQUIRE(FILE=vv_file,EXIST=ext)
        IF (.NOT. ext) THEN
          WRITE(*,*)
          WRITE(*,*) ' Volume fraction file ', vv_file(1:vv_file_l) ,' for ', min_name(1:min_name_l) ,' not found.'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        ! IF (ALLOCATED(VG_n)) THEN
        !   DEALLOCATE(VG_n)
        !   ALLOCATE(VG_n(nx, ny, nz))
        ! ELSE
        !   ALLOCATE(VG_n(nx, ny, nz))
        ! END IF
        OPEN(UNIT=23,FILE=vv_file,STATUS='old',ERR=8001)
        FileTemp = vv_file
        CALL stringlen(FileTemp,FileNameLength)
        IF (vv_fileformat == 'ContinuousRead') THEN
          READ(23,*,END=1020) (((volfx(min_id,jx,jy,jz),jx=1,nx),jy=1,ny),jz=1,nz)
        ELSE IF (vv_fileformat == 'SingleColumn') THEN
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) volfx(min_id,jx,jy,jz)
              END DO
            END DO
          END DO
        ELSE IF (vv_fileformat == 'FullForm') THEN
          IF (ny > 1 .AND. nz > 1) THEN
            DO jz = 1,nz
              DO jy = 1,ny
                DO jx= 1,nx
                  READ(23,*,END=1020) xdum,ydum,zdum,volfx(min_id,jx,jy,jz)
                END DO
              END DO
            END DO
          ELSE IF (ny > 1 .AND. nz == 1) THEN
            jz = 1
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) xdum,ydum,volfx(min_id,jx,jy,jz)
              END DO
            END DO
          ELSE
          jz = 1
          jy = 1
          DO jx= 1,nx
            READ(23,*,END=1020) xdum,volfx(min_id,jx,jy,jz)
          END DO
          END IF
        ELSE IF (vv_fileformat == 'Unformatted') THEN
        READ(23,END=1020) volfx
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Mineral volume fraction file format not recognized'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
  

        INQUIRE(FILE=bsa_file,EXIST=ext)
        IF (.NOT. ext) THEN
          WRITE(*,*)
          WRITE(*,*) ' Bulk surface area file ', bsa_file(1:bsa_file_l) ,' for ', min_name(1:min_name_l) ,' not found.'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        ! IF (ALLOCATED(VG_n)) THEN
        !   DEALLOCATE(VG_n)
        !   ALLOCATE(VG_n(nx, ny, nz))
        ! ELSE
        !   ALLOCATE(VG_n(nx, ny, nz))
        ! END IF
        OPEN(UNIT=23,FILE=bsa_file,STATUS='old',ERR=8001)
        FileTemp = bsa_file
        CALL stringlen(FileTemp,FileNameLength)
        IF (bsa_fileformat == 'ContinuousRead') THEN
          READ(23,*,END=1020) (((area(min_id,jx,jy,jz),jx=1,nx),jy=1,ny),jz=1,nz)
        ELSE IF (bsa_fileformat == 'SingleColumn') THEN
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) area(min_id,jx,jy,jz)
                IF (ssa_or_bsa == 1) THEN
                area(min_id,jx,jy,jz) = volfx(min_id,jx,jy,jz)*area(min_id,jx,jy,jz)*wtmin(min_id)/volmol(min_id)
                ENDIF
              END DO
            END DO
          END DO
        ELSE IF (bsa_fileformat == 'FullForm') THEN
          IF (ny > 1 .AND. nz > 1) THEN
            DO jz = 1,nz
              DO jy = 1,ny
                DO jx= 1,nx
                  READ(23,*,END=1020) xdum,ydum,zdum,area(min_id,jx,jy,jz)
                  IF (ssa_or_bsa == 1) THEN
                  area(min_id,jx,jy,jz) = volfx(min_id,jx,jy,jz)*area(min_id,jx,jy,jz)*wtmin(min_id)/volmol(min_id)
                  ENDIF
                END DO
              END DO
            END DO
          ELSE IF (ny > 1 .AND. nz == 1) THEN
            jz = 1
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) xdum,ydum,area(min_id,jx,jy,jz)
                IF (ssa_or_bsa == 1) THEN
                area(min_id,jx,jy,jz) = volfx(min_id,jx,jy,jz)*area(min_id,jx,jy,jz)*wtmin(min_id)/volmol(min_id)
                ENDIF
              END DO
            END DO
          ELSE
          jz = 1
          jy = 1
          DO jx= 1,nx
            READ(23,*,END=1020) xdum,area(min_id,jx,jy,jz)
            IF (ssa_or_bsa == 1) THEN
            area(min_id,jx,jy,jz) = volfx(min_id,jx,jy,jz)*area(min_id,jx,jy,jz)*wtmin(min_id)/volmol(min_id)
            ENDIF
          END DO
          END IF
        ELSE IF (bsa_fileformat == 'Unformatted') THEN
        READ(23,END=1020) area
        IF (ssa_or_bsa == 1) THEN
        area = volfx*area*wtmin(min_id)/volmol(min_id)
        ENDIF
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Bulk surface area file format not recognized'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
  !stop
ENDDO

ENDIF
!*****************************
!Stolze Lucien, June 2023
!END: read mineral volume fraction and bulk surface area from file

!! Or overwrite with read of geochemical conditions
IF (ReadGeochemicalConditions) THEN

!!    INQUIRE(FILE='calcite.dat',EXIST=ext)

!!    IF (ext) THEN

!! Hard wired here for calcite as Mineral 1

!!      OPEN(unit=99,file='calcite.dat',status='old')

!!      read(99,*) GridCoordinateX,GridCoordinateY
!!      do jy = 1,ny
!!        do jx = 1,nx
!!          if (jx == GridCoordinateX .AND. jy == GridCoordinateY) then
!!            volfx(1,jx,jy,1) = 0.05
!!            read(99,*,END=981) GridCoordinateX,GridCoordinateY
!!          else
!!            volfx(1,jx,jy,1) = 0.000
!!          endif
!!        end do
!!      end do
!!      close(unit=99,status='keep')
!!981   write(*,*) ' End of calcite.dat file'

!!    ELSE
!!      continue
!!    END IF

!! Plagioclase read (mineral 6)

  INQUIRE(FILE='plagioclase.dat',EXIST=ext)

  IF (ext) THEN

!! Hard wired here for plagioclase as Mineral 6

    OPEN(unit=97,file='plagioclase.dat',status='old')

    read(97,*) GridCoordinateX,GridCoordinateY
    do jy = 1,ny
      do jx = 1,nx
        if (jx == GridCoordinateX .AND. jy == GridCoordinateY) then
          if (jinit(jx,jy,1) == 2) then
            volfx(6,jx,jy,1) = 0.00
          else
            volfx(6,jx,jy,1) = 0.20
          end if
          read(97,*,END=979) GridCoordinateX,GridCoordinateY
        else
          volfx(6,jx,jy,1) = 0.000
        endif
      end do
    end do
    close(unit=97,status='keep')
979   write(*,*) ' End of plagiclase.dat file'

  ELSE
    continue
  END IF


END IF    !! End of ReadGeochemicalConditions = .TRUE.

IF (ALLOCATED(work3)) then
DEALLOCATE(work3)
END IF

! Fill in ghost cells

text = 'Porosity'
lowX = LBOUND(por,1)
lowY = LBOUND(por,2)
lowZ = LBOUND(por,3)
highX = UBOUND(por,1)
highY = UBOUND(por,2)
highZ = UBOUND(por,3)
call GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,por,TEXT)

text = 'Old Porosity'
lowX = LBOUND(porin,1)
lowY = LBOUND(porin,2)
lowZ = LBOUND(porin,3)
highX = UBOUND(porin,1)
highY = UBOUND(porin,2)
highZ = UBOUND(porin,3)
call GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,porin,TEXT)

text = 'Density'
lowX = LBOUND(ro,1)
lowY = LBOUND(ro,2)
lowZ = LBOUND(ro,3)
highX = UBOUND(ro,1)
highY = UBOUND(ro,2)
highZ = UBOUND(ro,3)
call GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,ro,TEXT)

text = 'Old Density'
lowX = LBOUND(roold,1)
lowY = LBOUND(roold,2)
lowZ = LBOUND(roold,3)
highX = UBOUND(roold,1)
highY = UBOUND(roold,2)
highZ = UBOUND(roold,3)
call GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,roold,TEXT)

text = 'Saturation'
lowX = LBOUND(satliq,1)
lowY = LBOUND(satliq,2)
lowZ = LBOUND(satliq,3)
highX = UBOUND(satliq,1)
highY = UBOUND(satliq,2)
highZ = UBOUND(satliq,3)
call GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,satliq,TEXT)

text = 'Old Saturation'
lowX = LBOUND(satliqold,1)
lowY = LBOUND(satliqold,2)
lowZ = LBOUND(satliqold,3)
highX = UBOUND(satliqold,1)
highY = UBOUND(satliqold,2)
highZ = UBOUND(satliqold,3)
call GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,satliqold,TEXT)

text = 'Gas Density'
lowX = LBOUND(rogas,1)
lowY = LBOUND(rogas,2)
lowZ = LBOUND(rogas,3)
highX = UBOUND(rogas,1)
highY = UBOUND(rogas,2)
highZ = UBOUND(rogas,3)
call GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,rogas,TEXT)

spold = sp
spexold = spex
spsurfold = spsurf
spnO2(1:nx,1:ny,1:nz) = sp(ikmast,1:nx,1:ny,1:nz)
spnnO2(1:nx,1:ny,1:nz) = sp(ikmast,1:nx,1:ny,1:nz)

!!!    ***************BOUNDARY CONDITIONS********************

section = 'boundary_conditions'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
WRITE(*,*) ' Boundary condition block found'
ELSE
WRITE(*,*)
WRITE(*,*) ' No boundary conditions found'
WRITE(*,*)
END IF

IF (ALLOCATED(JcByGrid)) THEN
DEALLOCATE(JcByGrid)
END IF
ALLOCATE(JcByGrid(0:nx+1,0:ny+1,0:nz+1))
JcByGrid = 0

CALL read_BoundaryConditionByZone(nout,nx,ny,nz,nBoundaryConditionZone)

IF (nBoundaryConditionZone > 0) THEN

!!  Initialize concentration at boundaries from various zones

DO l = 1,nBoundaryConditionZone

  DO jz = jzzBC_lo(l),jzzBC_hi(l)
    DO jy = jyyBC_lo(l),jyyBC_hi(l)
      DO jx = jxxBC_lo(l),jxxBC_hi(l)

        ConditionName = BoundaryConditionName(l)
        ConditionNameFound = .FALSE.
        ConditionNumber = 0

        DO nco = 1,nchem
          IF (ConditionName == CondLabel(nco)) THEN
            ConditionNameFound = .TRUE.
            ConditionNumber = nco
          END IF
          IF (ConditionNameFound) THEN
            EXIT
          END IF
        END DO

        IF (.NOT. ConditionNameFound) THEN
          WRITE(*,*)
          WRITE(*,*) ' Condition block not found for boundary input read'
          WRITE(*,*) ' Looking for:  ', ConditionName
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF

        jinit(jx,jy,jz) = ConditionNumber

        JcByGrid(jx,jy,jz) = JcByBoundaryZone(l)

      END DO
    END DO
  END DO

END DO

DEALLOCATE (JcByBoundaryZone)

jx = 0
DO jz = 1,nz
  DO jy = 1,ny

    ConditionNumber = jinit(jx,jy,jz)

    IF (ConditionNumber == 0) THEN
      WRITE(*,*)
      WRITE(*,*) ' Portion of JX = 0 boundary is not initialized'
      WRITE(*,*) ' Missing initialization at jy = ',jy
      WRITE(*,*) ' Missing initialization at jz = ',jz
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF


    DO ik = 1,ncomp+nspec
      sp10(ik,jx,jy,jz)    = spcond10(ik,ConditionNumber)
      sp(ik,jx,jy,jz)      = spcond(ik,ConditionNumber)
    END DO

    DO i = 1,ncomp
      s(i,jx,jy,jz)        = scond(i,ConditionNumber)
    END DO

    DO kk = 1,ngas
      spgas10(kk,jx,jy,jz) = spcondgas10(kk,ConditionNumber)
      spgas(kk,jx,jy,jz)   = spcondgas(kk,ConditionNumber)
    END DO

    do ix = 1,nexchange
      spex(ix,jx,jy,jz)    = spcondex(ix,ConditionNumber)
    end do
    DO ix = 1,nexchange+nexch_sec
      spex10(ix+nexchange,jx,jy,jz) = convert*spcondex10(ix+nexchange,ConditionNumber)  ! Now in eq/m3 por. med.
    END DO

    DO is = 1,nsurf
      spsurf(is,jx,jy,jz)   = LOG(convert*spcondsurf10(is,ConditionNumber))
      IF (iedl(is) == 0) THEN
        LogPotential(is,jx,jy,jz) = LogPotentialInit(is,ConditionNumber)
      END IF
    END DO
    DO is = 1,nsurf+nsurf_sec
      spsurf10(is,jx,jy,jz) = convert*spcondsurf10(is,ConditionNumber)
    END DO

    t(jx,jy,jz) = tempcond(ConditionNumber)
    por(jx,jy,jz) = porcond(ConditionNumber)

  END DO
END DO

jx = nx+1
DO jz = 1,nz
  DO jy = 1,ny

    ConditionNumber = jinit(jx,jy,jz)

    IF (ConditionNumber == 0) THEN
      WRITE(*,*)
      WRITE(*,*) ' Portion of JX = nx+1 boundary is not initialized'
      WRITE(*,*) ' Missing initialization at jy = ',jy
      WRITE(*,*) ' Missing initialization at jz = ',jz
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

    DO ik = 1,ncomp+nspec
      sp10(ik,jx,jy,jz)    = spcond10(ik,ConditionNumber)
      sp(ik,jx,jy,jz)      = spcond(ik,ConditionNumber)
    END DO

    DO i = 1,ncomp
      s(i,jx,jy,jz)        = scond(i,ConditionNumber)
    END DO

    DO kk = 1,ngas
      spgas10(kk,jx,jy,jz) = spcondgas10(kk,ConditionNumber)
      spgas(kk,jx,jy,jz)   = spcondgas(kk,ConditionNumber)
    END DO

    do ix = 1,nexchange
      spex(ix,jx,jy,jz)    = spcondex(ix,ConditionNumber)
    end do
    DO ix = 1,nexchange+nexch_sec
      spex10(ix+nexchange,jx,jy,jz) = convert*spcondex10(ix+nexchange,ConditionNumber)  ! Now in eq/m3 por. med.
    END DO

    DO is = 1,nsurf
      spsurf(is,jx,jy,jz)   = LOG(convert*spcondsurf10(is,ConditionNumber))
      IF (iedl(is) == 0) THEN
        LogPotential(is,jx,jy,jz) = LogPotentialInit(is,ConditionNumber)
      END IF
    END DO
    DO is = 1,nsurf+nsurf_sec
      spsurf10(is,jx,jy,jz) = convert*spcondsurf10(is,ConditionNumber)
    END DO

    t(jx,jy,jz) = tempcond(ConditionNumber)
    por(jx,jy,jz) = porcond(ConditionNumber)

  END DO
END DO

IF (ny > 1) THEN

  jy = 0
  DO jz = 1,nz
    DO jx = 1,nx

      ConditionNumber = jinit(jx,jy,jz)

      IF (ConditionNumber == 0) THEN
        WRITE(*,*)
        WRITE(*,*) ' Portion of JY = 0 boundary is not initialized'
        WRITE(*,*) ' Missing initialization at jx = ',jx
        WRITE(*,*) ' Missing initialization at jz = ',jz
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

      DO ik = 1,ncomp+nspec
        sp10(ik,jx,jy,jz)    = spcond10(ik,ConditionNumber)
        sp(ik,jx,jy,jz)      = spcond(ik,ConditionNumber)
      END DO

      DO i = 1,ncomp
        s(i,jx,jy,jz)        = scond(i,ConditionNumber)
      END DO

      DO kk = 1,ngas
        spgas10(kk,jx,jy,jz) = spcondgas10(kk,ConditionNumber)
        spgas(kk,jx,jy,jz)   = spcondgas(kk,ConditionNumber)
      END DO

      do ix = 1,nexchange
        spex(ix,jx,jy,jz)    = spcondex(ix,ConditionNumber)
      end do
      DO ix = 1,nexchange+nexch_sec
        spex10(ix+nexchange,jx,jy,jz) = convert*spcondex10(ix+nexchange,ConditionNumber)  ! Now in eq/m3 por. med.
      END DO

      DO is = 1,nsurf
        spsurf(is,jx,jy,jz)   = LOG(convert*spcondsurf10(is,ConditionNumber))
        IF (iedl(is) == 0) THEN
          LogPotential(is,jx,jy,jz) = LogPotentialInit(is,ConditionNumber)
        END IF
      END DO
      DO is = 1,nsurf+nsurf_sec
        spsurf10(is,jx,jy,jz) = convert*spcondsurf10(is,ConditionNumber)
      END DO

      t(jx,jy,jz) = tempcond(ConditionNumber)
      por(jx,jy,jz) = porcond(ConditionNumber)

    END DO
  END DO

  jy = ny+1
  DO jz = 1,nz
    DO jx = 1,nx

      ConditionNumber = jinit(jx,jy,jz)

      IF (ConditionNumber == 0) THEN
        WRITE(*,*)
        WRITE(*,*) ' Portion of JY = ny+1 boundary is not initialized'
        WRITE(*,*) ' Missing initialization at jx = ',jx
        WRITE(*,*) ' Missing initialization at jz = ',jz
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

      DO ik = 1,ncomp+nspec
        sp10(ik,jx,jy,jz)    = spcond10(ik,ConditionNumber)
        sp(ik,jx,jy,jz)      = spcond(ik,ConditionNumber)
      END DO

      DO i = 1,ncomp
        s(i,jx,jy,jz)        = scond(i,ConditionNumber)
      END DO

      DO kk = 1,ngas
        spgas10(kk,jx,jy,jz) = spcondgas10(kk,ConditionNumber)
        spgas(kk,jx,jy,jz)   = spcondgas(kk,ConditionNumber)
      END DO

      do ix = 1,nexchange
        spex(ix,jx,jy,jz)    = spcondex(ix,ConditionNumber)
      end do
      DO ix = 1,nexchange+nexch_sec
        spex10(ix+nexchange,jx,jy,jz) = convert*spcondex10(ix+nexchange,ConditionNumber)  ! Now in eq/m3 por. med.
      END DO

      DO is = 1,nsurf
        spsurf(is,jx,jy,jz)   = LOG(convert*spcondsurf10(is,ConditionNumber))
        IF (iedl(is) == 0) THEN
          LogPotential(is,jx,jy,jz) = LogPotentialInit(is,ConditionNumber)
        END IF
      END DO
      DO is = 1,nsurf+nsurf_sec
        spsurf10(is,jx,jy,jz) = convert*spcondsurf10(is,ConditionNumber)
      END DO

      t(jx,jy,jz) = tempcond(ConditionNumber)
      por(jx,jy,jz) = porcond(ConditionNumber)


    END DO
  END DO

END IF   !!  End of NY block

IF (nz > 1) THEN

  jz = 0
  DO jy = 1,ny
    DO jx = 1,nx

      ConditionNumber = jinit(jx,jy,jz)

      IF (ConditionNumber == 0) THEN
        WRITE(*,*)
        WRITE(*,*) ' Portion of JZ = 0 boundary is not initialized'
        WRITE(*,*) ' Missing initialization at jx = ',jx
        WRITE(*,*) ' Missing initialization at jy = ',jy
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

      DO ik = 1,ncomp+nspec
        sp10(ik,jx,jy,jz)    = spcond10(ik,ConditionNumber)
        sp(ik,jx,jy,jz)      = spcond(ik,ConditionNumber)
      END DO

      DO i = 1,ncomp
        s(i,jx,jy,jz)        = scond(i,ConditionNumber)
      END DO

      DO kk = 1,ngas
        spgas10(kk,jx,jy,jz) = spcondgas10(kk,ConditionNumber)
        spgas(kk,jx,jy,jz)   = spcondgas(kk,ConditionNumber)
      END DO

      do ix = 1,nexchange
        spex(ix,jx,jy,jz)    = spcondex(ix,ConditionNumber)
      end do
      DO ix = 1,nexchange+nexch_sec
        spex10(ix+nexchange,jx,jy,jz) = convert*spcondex10(ix+nexchange,ConditionNumber)  ! Now in eq/m3 por. med.
      END DO

      DO is = 1,nsurf
        spsurf(is,jx,jy,jz)   = LOG(convert*spcondsurf10(is,ConditionNumber))
        IF (iedl(is) == 0) THEN
          LogPotential(is,jx,jy,jz) = LogPotentialInit(is,ConditionNumber)
        END IF
      END DO
      DO is = 1,nsurf+nsurf_sec
        spsurf10(is,jx,jy,jz) = convert*spcondsurf10(is,ConditionNumber)
      END DO

      t(jx,jy,jz) = tempcond(ConditionNumber)
      por(jx,jy,jz) = porcond(ConditionNumber)


    END DO
  END DO

  jz = nz+1
  DO jy = 1,ny
    DO jx = 1,nx

      ConditionNumber = jinit(jx,jy,jz)

      IF (ConditionNumber == 0) THEN
        WRITE(*,*)
        WRITE(*,*) ' Portion of JZ = nz+1 boundary is not initialized'
        WRITE(*,*) ' Missing initialization at jx = ',jx
        WRITE(*,*) ' Missing initialization at jy = ',jy
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

      DO ik = 1,ncomp+nspec
        sp10(ik,jx,jy,jz)    = spcond10(ik,ConditionNumber)
        sp(ik,jx,jy,jz)      = spcond(ik,ConditionNumber)
      END DO

      DO i = 1,ncomp
        s(i,jx,jy,jz)        = scond(i,ConditionNumber)
      END DO

      DO kk = 1,ngas
        spgas10(kk,jx,jy,jz) = spcondgas10(kk,ConditionNumber)
        spgas(kk,jx,jy,jz)   = spcondgas(kk,ConditionNumber)
      END DO

      do ix = 1,nexchange
        spex(ix,jx,jy,jz)    = spcondex(ix,ConditionNumber)
      end do
      DO ix = 1,nexchange+nexch_sec
        spex10(ix+nexchange,jx,jy,jz) = convert*spcondex10(ix+nexchange,ConditionNumber)  ! Now in eq/m3 por. med.
      END DO

      DO is = 1,nsurf
        spsurf(is,jx,jy,jz)   = LOG(convert*spcondsurf10(is,ConditionNumber))
        IF (iedl(is) == 0) THEN
          LogPotential(is,jx,jy,jz) = LogPotentialInit(is,ConditionNumber)
        END IF
      END DO
      DO is = 1,nsurf+nsurf_sec
        spsurf10(is,jx,jy,jz) = convert*spcondsurf10(is,ConditionNumber)
      END DO

      t(jx,jy,jz) = tempcond(ConditionNumber)
      por(jx,jy,jz) = porcond(ConditionNumber)


    END DO
  END DO

END IF    !! End of NZ block



ELSE   !! Conventional treatment of boundaries as corresponding to an entire boundary face

DO i = 1,6
  jc(i) = 2   ! Initialize (and make default) a flux boundary
END DO

IF (ALLOCATED(spb)) THEN
  DEALLOCATE(spb)
  ALLOCATE(spb(ncomp+nspec,6))
ELSE
  ALLOCATE(spb(ncomp+nspec,6))
END IF
IF (ALLOCATED(spbgas)) THEN
  DEALLOCATE(spbgas)
  ALLOCATE(spbgas(ngas,6))
ELSE
  ALLOCATE(spbgas(ngas,6))
END IF
IF (ALLOCATED(spexb)) THEN
  DEALLOCATE(spexb)
  ALLOCATE(spexb(nexchange+nexch_sec,6))
ELSE
  ALLOCATE(spexb(nexchange+nexch_sec,6))
END IF
IF (ALLOCATED(spsurfb)) THEN
  DEALLOCATE(spsurfb)
  ALLOCATE(spsurfb(nsurf+nsurf_sec,6))
ELSE
  ALLOCATE(spsurfb(nsurf+nsurf_sec,6))
END IF

spb = 0.0
spbgas = 0.0
spexb = 0.0
spsurfb = 0.0

! if (Richards) THEN
! CALL read_bound_richard(nout,nchem,nx,ny,nz,ncomp,nspec,ngas,nkin,nexchange,nexch_sec,  &
!    nsurf,nsurf_sec)
! ELSE
CALL read_bound(nout,nchem,nx,ny,nz,ncomp,nspec,ngas,nkin,nexchange,nexch_sec,  &
   nsurf,nsurf_sec)
!ENDIF

IF (ALLOCATED(AqueousToBulkCond)) THEN
  DEALLOCATE(AqueousToBulkCond)
END IF

IF (ALLOCATED(sbnd)) THEN
  DEALLOCATE(sbnd)
  ALLOCATE(sbnd(ncomp,6))
ELSE
  ALLOCATE(sbnd(ncomp,6))
END IF
IF (ALLOCATED(s_local)) THEN
  DEALLOCATE(s_local)
  ALLOCATE(s_local(ncomp))
ELSE
  ALLOCATE(s_local(ncomp))
END IF

sbnd = 0.0
s_local = 0.0

DO nbnd = 1,6
  CALL bdcalc(ncomp,nspec,nbnd)
  DO i = 1,ncomp
    sbnd(i,nbnd) = s_local(i)
  END DO
END DO

!!!  Set the boundaries

jx = 0
DO jz = 1,nz
  DO jy = 1,ny
    DO i = 1,ncomp
      s(i,jx,jy,jz) = sbnd(i,1)
    END DO
  END DO
END DO

jx = nx+1
DO jz = 1,nz
  DO jy = 1,ny
    DO i = 1,ncomp
      s(i,jx,jy,jz) = sbnd(i,2)
    END DO
  END DO
END DO

jy = 0
DO jz = 1,nz
  DO jx = 1,nx
    DO i = 1,ncomp
      s(i,jx,jy,jz) = sbnd(i,3)
    END DO
  END DO
END DO

jy = ny+1
DO jz = 1,nz
  DO jx = 1,nx
    DO i = 1,ncomp
      s(i,jx,jy,jz) = sbnd(i,4)
    END DO
  END DO
END DO

END IF

!  ****************END OF BOUNDARY CONDITIONS*************

!     ***********OUTPUT SECTION***************************
section = 'output'
CALL readblock(nin,nout,section,found,ncount)

OutputTimeUnits = 'years'

IF (found) THEN
WRITE(*,*) ' Output block found'
WRITE(*,*)

CALL units_timeOutput(nout,section,time_scale)

nstop = 0
parchar = 'spatial_profile'
parchar2 = 'spatial_profile_at_time'
parfind = ' '
realmult= 0.0
CALL read_snapshot(nout,lchar,parchar,parchar2,parfind,realmult,lenarray,section)

IF (parfind == 'spatial_profile' .OR. parfind == 'spatial_profile_at_time') THEN
  CONTINUE
ELSE
  parchar = 'read_snapshotfile'
  parfind = ' '
  CALL readFileName(nout,lchar,parchar,parfind,dumstring,section,SnapshotFileFormat)
  IF (parfind == 'read_snapshotfile') THEN
  lenarray=0
  OPEN(UNIT=23,FILE=dumstring,STATUS='old')
  FileTemp = dumstring
  CALL stringlen(FileTemp,FileNameLength)
  do while (ierr == 0)
    lenarray = lenarray + 1
    READ(23,*,iostat=IERR) realmult_dum(lenarray) 
  enddo
  realmult=realmult_dum(1:lenarray-1)
  ENDIF
ENDIF

IF (parfind == ' ') THEN
  WRITE(*,*) ' Timestepping off--initialization only'
  nstop = 0
ELSE
  IF (realmult(1) <= 0.0) THEN
    WRITE(*,*) ' Timestepping turned on, but no time provided'
    WRITE(*,*) ' Doing initialization only'
    nstop = 0
  ELSE
    IF (parfind == 'read_snapshotfile') THEN
    nstop = (lenarray-1)
    ELSE
    nstop = lenarray
    ENDIF

    IF (ALLOCATED(prtint)) THEN
      DEALLOCATE(prtint)
      ALLOCATE(prtint(nstop))
    ELSE
      ALLOCATE(prtint(nstop))
    END IF

    WRITE(*,*) ' Timestepping on--output files written'
    DO i = 1,nstop
      prtint(i) = realmult(i)
      WRITE(*,*) realmult(i)
    END DO
    WRITE(*,*)
  END IF
END IF

parchar2 = ''

OutputTimeScale = 1.0d0/time_scale

!  Convert units if necessary

DO i = 1,nstop
  prtint(i) = prtint(i)*time_scale
END DO

!!  Beginning of minseries definition ********************

!! This reads the file names and the location of the mineral time series curves
minseries = 0
CALL read_minseries(nout,minseries,nx,ny,nz)

DO ll = 1,minseries

  IF (jxminseries(ll) > nx) THEN
    WRITE(*,*)
    WRITE(*,*) '  You have specified a mineral series plot location > NX'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (jyminseries(ll) > ny) THEN
    WRITE(*,*)
    WRITE(*,*) '  You have specified a plot location > NY'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (jzminseries(ll) > nz) THEN
    WRITE(*,*)
    WRITE(*,*) '  You have specified a plot location > NZ'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END DO


!!  End of minseries definition ********************

!! This reads the file names and the location of the breakthrough curves
nseries = 0
CALL read_series(nout,nseries,nx,ny,nz)

iplotall = 0

!! This reads the species to be written to the breakthrough files

parchar = 'time_series_print'
parfind = ' '
lenarray = 0
nplot = 0

CALL read_multstring(nout,lchar,parchar,parfind,stringarray,lenarray,section)
IF (parfind == ' ') THEN    ! Assume default
  WRITE(*,*) ' Parfind is blank'
  WRITE(*,*)
  WRITE(*,*) ' Tracking all species at time series locations'
  WRITE(*,*)
  iplotall = 1
  nplot = ncomp
  DO ll = 1,nplot
    iplot(ll) = ll
  END DO
ELSE
  IF (stringarray(1) == ' ') THEN
    WRITE(*,*) ' First string is blank'
    WRITE(*,*)
    WRITE(*,*) ' Tracking all species at time series locations'
    WRITE(*,*)
    iplotall = 1
    nplot = ncomp
    DO ll = 1,nplot
      iplot(ll) = ll
    END DO
  ELSE IF (stringarray(1) == 'all') THEN
    iplotall = 1
    WRITE(*,*)
    WRITE(*,*) ' Tracking all species at time series locations'
    WRITE(*,*)
    nplot = ncomp
    DO ll = 1,nplot
      iplot(ll) = ll
    END DO
  ELSE

!  Check to see that strings match species names

    iplotph = 0
    nplot = lenarray
    DO ll = 1,nplot
      iplot(ll) = 0
      IF (stringarray(ll) == 'pH' .OR. stringarray(ll) == 'ph') THEN
        IF (ikph /= 0) THEN
          iplot(ll) = ikph
          iplotph = ll
        ELSE
          WRITE(*,*) '       ERROR '
          WRITE(*,*) ' pH not included in problem, but you are asking to write pH in'
          WRITE(*,*) '   time series file'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
!!!        ELSE IF (stringarray(ll) == 'minerals' .OR. stringarray(ll) == 'Minerals' .OR. stringarray(ll) == 'MINERALS') THEN
!!!          iplotMinerals = .TRUE.
      ELSE
        DO ik = 1,ncomp
          IF (ulab(ik) == stringarray(ll)) THEN
            iplot(ll) = ik
          END IF
        END DO
        IF (iplot(ll) == 0) THEN
          dumstring = stringarray(ll)
          CALL stringlen(dumstring,nlength)
          WRITE(*,*)
          WRITE(*,*) ' Plotting species not found in list: ',dumstring(1:nlength)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      END IF
    END DO
  END IF
END IF

IF (ALLOCATED(TimeSeriesSpecies)) THEN
  DEALLOCATE(TimeSeriesSpecies)
END IF
ALLOCATE(TimeSeriesSpecies(nplot))
TimeSeriesSpecies = stringarray(1:nplot)

parchar = 'time_series_units'
parfind = ' '
lenarray = 0

CALL read_multstring(nout,lchar,parchar,parfind,stringarray,lenarray,section)
IF (parfind == ' ') THEN    ! No units provided
  DO ll = 1,nplot
    stringarray(ll) = ' '
  END DO
  IF (ALLOCATED(TimeSeriesUnits)) THEN
    DEALLOCATE(TimeSeriesUnits)
  END IF
  ALLOCATE(TimeSeriesUnits(nplot))
  DO ll = 1,nplot
    TimeSeriesUnits(ll) = stringarray(ll)
  END DO
ELSE
  IF (ALLOCATED(TimeSeriesUnits)) THEN
    DEALLOCATE(TimeSeriesUnits)
  END IF
  ALLOCATE(TimeSeriesUnits(nplot))
  DO ll = 1,nplot
    TimeSeriesUnits(ll) = stringarray(ll)
  END DO
END IF

parchar = 'time_series_output'
parfind = ' '
CALL read_multpar(nout,lchar,parchar,parfind,realmult,lenarray,section)
IF (parfind == ' ') THEN    !  If specific output not found, then look for an interval

  parchar = 'time_series_interval'
  parfind = ' '
  intjunk = 0
  CALL read_integer(nout,lchar,parchar,parfind,intjunk,section)
  IF (parfind == ' ') THEN  ! Parameter to set time series interval not found
    interval = 1            ! Use default
  ELSE
    IF (intjunk == 0) THEN
      interval = 1         ! Use default
    ELSE
      interval = intjunk
    END IF
  END IF

  IF (ALLOCATED(OutputTime)) THEN
    DEALLOCATE(OutputTime)
  END IF
  ALLOCATE(OutputTime(1))

  OutputTime(1) = 0.0
ELSE
  NumOutputTimes = lenarray

  IF (ALLOCATED(OutputTime)) THEN
    DEALLOCATE(OutputTime)
  END IF
  ALLOCATE(OutputTime(NumOutputTimes+1))

  DO i = 1,NumOutputTimes
    OutputTime(i) = realmult(i)*time_scale
  END DO
  interval = 0
END IF

time_scale = 1.0d0

DO ll = 1,nseries
  IF (jxseries(ll) > nx) THEN
    WRITE(*,*)
    WRITE(*,*) '  You have specified a plot location > NX'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (jyseries(ll) > ny) THEN
    WRITE(*,*)
    WRITE(*,*) '  You have specified a plot location > NY'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (jzseries(ll) > nz) THEN
    WRITE(*,*)
    WRITE(*,*) '  You have specified a plot location > NZ'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END DO

DO ll = 1,nisotopeseries
  IF (jxisotopeseries(ll) > nx) THEN
    WRITE(*,*)
    WRITE(*,*) '  You have specified an isotope plot location > NX'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (jyisotopeseries(ll) > ny) THEN
    WRITE(*,*)
    WRITE(*,*) '  You have specified an isotope plot location > NY'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  IF (jzisotopeseries(ll) > nz) THEN
    WRITE(*,*)
    WRITE(*,*) '  You have specified an isotope plot location > NZ'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END DO

ELSE

WRITE(*,*)
WRITE(*,*) '        NO OUTPUT TIMES SPECIFIED'
WRITE(*,*) '   TIMESTEPPING OFF--INITIALIZATION ONLY'
WRITE(*,*)
nstop = 0
IF (ALLOCATED(OutputTime)) THEN
  DEALLOCATE(OutputTime)
END IF
ALLOCATE(OutputTime(1))
OutputTime(1) = 0.0
END IF

MakeMovie = .FALSE.
parchar = 'MakeMovie'
parfind = ' '
CALL read_logical(nout,lchar,parchar,parfind,MakeMovie)

call read_AqueousFluxSeries(nout,ncomp,nx,ny,nz)

call read_FluxWeightedConcentration(nout,ncomp,nx,ny,nz)

i = size(jxseries,1)
ALLOCATE(workint1(i))
workint1 = jxseries
DEALLOCATE(jxseries)
ALLOCATE(jxseries(nseries))
IF(nseries /= 0) jxseries(1:nseries) = workint1(1:nseries)
DEALLOCATE(workint1)

i = size(jyseries,1)
ALLOCATE(workint1(i))
workint1 = jyseries
DEALLOCATE(jyseries)
ALLOCATE(jyseries(nseries))
IF(nseries /= 0)jyseries(1:nseries) = workint1(1:nseries)
DEALLOCATE(workint1)

i = size(jzseries,1)
ALLOCATE(workint1(i))
workint1 = jzseries
DEALLOCATE(jzseries)
ALLOCATE(jzseries(nseries))
IF(nseries /= 0)jzseries(1:nseries) = workint1(1:nseries)
DEALLOCATE(workint1)

!!!  ***********************
!!!  Add jxminseries
i = size(jxminseries,1)
ALLOCATE(workint1(i))
workint1 = jxminseries
DEALLOCATE(jxminseries)
ALLOCATE(jxminseries(minseries))
IF(minseries /= 0) jxminseries(1:minseries) = workint1(1:minseries)
DEALLOCATE(workint1)

i = size(jyminseries,1)
ALLOCATE(workint1(i))
workint1 = jyminseries
DEALLOCATE(jyminseries)
ALLOCATE(jyminseries(minseries))
IF(minseries /= 0)jyminseries(1:minseries) = workint1(1:minseries)
DEALLOCATE(workint1)

i = size(jzminseries,1)
ALLOCATE(workint1(i))
workint1 = jzminseries
DEALLOCATE(jzminseries)
ALLOCATE(jzminseries(minseries))
IF(minseries /= 0)jzminseries(1:minseries) = workint1(1:minseries)
DEALLOCATE(workint1)

!!!  *********************

i = size(jxAqueousFluxSeries_lo,1)
ALLOCATE(workint1(i))
workint1 = jxAqueousFluxSeries_lo
DEALLOCATE(jxAqueousFluxSeries_lo)
ALLOCATE(jxAqueousFluxSeries_lo(nAqueousFluxSeriesFile))
IF(nAqueousFluxSeriesFile /= 0) jxAqueousFluxSeries_lo(1:nAqueousFluxSeriesFile) = workint1(1:nAqueousFluxSeriesFile)
DEALLOCATE(workint1)

i = size(jyAqueousFluxSeries_lo,1)
ALLOCATE(workint1(i))
workint1 = jyAqueousFluxSeries_lo
DEALLOCATE(jyAqueousFluxSeries_lo)
ALLOCATE(jyAqueousFluxSeries_lo(nAqueousFluxSeriesFile))
IF(nAqueousFluxSeriesFile /= 0) jyAqueousFluxSeries_lo(1:nAqueousFluxSeriesFile) = workint1(1:nAqueousFluxSeriesFile)
DEALLOCATE(workint1)

i = size(jzAqueousFluxSeries_lo,1)
ALLOCATE(workint1(i))
workint1 = jzAqueousFluxSeries_lo
DEALLOCATE(jzAqueousFluxSeries_lo)
ALLOCATE(jzAqueousFluxSeries_lo(nAqueousFluxSeriesFile))
IF(nseries /= 0) jzAqueousFluxSeries_lo(1:nAqueousFluxSeriesFile) = workint1(1:nAqueousFluxSeriesFile)
DEALLOCATE(workint1)

i = size(jxAqueousFluxSeries_hi,1)
ALLOCATE(workint1(i))
workint1 = jxAqueousFluxSeries_hi
DEALLOCATE(jxAqueousFluxSeries_hi)
ALLOCATE(jxAqueousFluxSeries_hi(nAqueousFluxSeriesFile))
IF(nAqueousFluxSeriesFile /= 0) jxAqueousFluxSeries_hi(1:nAqueousFluxSeriesFile) = workint1(1:nAqueousFluxSeriesFile)
DEALLOCATE(workint1)

i = size(jyAqueousFluxSeries_hi,1)
ALLOCATE(workint1(i))
workint1 = jyAqueousFluxSeries_hi
DEALLOCATE(jyAqueousFluxSeries_hi)
ALLOCATE(jyAqueousFluxSeries_hi(nAqueousFluxSeriesFile))
IF(nAqueousFluxSeriesFile /= 0) jyAqueousFluxSeries_hi(1:nAqueousFluxSeriesFile) = workint1(1:nAqueousFluxSeriesFile)
DEALLOCATE(workint1)

i = size(jzAqueousFluxSeries_hi,1)
ALLOCATE(workint1(i))
workint1 = jzAqueousFluxSeries_hi
DEALLOCATE(jzAqueousFluxSeries_hi)
ALLOCATE(jzAqueousFluxSeries_hi(nAqueousFluxSeriesFile))
IF(nAqueousFluxSeriesFile /= 0) jzAqueousFluxSeries_hi(1:nAqueousFluxSeriesFile) = &
              workint1(1:nAqueousFluxSeriesFile)
DEALLOCATE(workint1)


i = size(jxFluxWeightedConcentration_lo,1)
ALLOCATE(workint1(i))
workint1 = jxFluxWeightedConcentration_lo
DEALLOCATE(jxFluxWeightedConcentration_lo)
ALLOCATE(jxFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile))
IF(nFluxWeightedConcentrationFile /= 0) jxFluxWeightedConcentration_lo(1:nFluxWeightedConcentrationFile) = &
              workint1(1:nFluxWeightedConcentrationFile)
DEALLOCATE(workint1)

i = size(jyFluxWeightedConcentration_lo,1)
ALLOCATE(workint1(i))
workint1 = jyFluxWeightedConcentration_lo
DEALLOCATE(jyFluxWeightedConcentration_lo)
ALLOCATE(jyFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile))
IF(nFluxWeightedConcentrationFile /= 0) jyFluxWeightedConcentration_lo(1:nFluxWeightedConcentrationFile) = &
              workint1(1:nFluxWeightedConcentrationFile)
DEALLOCATE(workint1)

i = size(jzFluxWeightedConcentration_lo,1)
ALLOCATE(workint1(i))
workint1 = jzFluxWeightedConcentration_lo
DEALLOCATE(jzFluxWeightedConcentration_lo)
ALLOCATE(jzFluxWeightedConcentration_lo(nFluxWeightedConcentrationFile))
IF(nFluxWeightedConcentrationFile /= 0) jzFluxWeightedConcentration_lo(1:nFluxWeightedConcentrationFile) = &
              workint1(1:nFluxWeightedConcentrationFile)
DEALLOCATE(workint1)

i = size(jxFluxWeightedConcentration_hi,1)
ALLOCATE(workint1(i))
workint1 = jxFluxWeightedConcentration_hi
DEALLOCATE(jxFluxWeightedConcentration_hi)
ALLOCATE(jxFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile))
IF(nFluxWeightedConcentrationFile /= 0) jxFluxWeightedConcentration_hi(1:nFluxWeightedConcentrationFile) = &
              workint1(1:nFluxWeightedConcentrationFile)
DEALLOCATE(workint1)

i = size(jyFluxWeightedConcentration_hi,1)
ALLOCATE(workint1(i))
workint1 = jyFluxWeightedConcentration_hi
DEALLOCATE(jyFluxWeightedConcentration_hi)
ALLOCATE(jyFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile))
IF(nFluxWeightedConcentrationFile /= 0) jyFluxWeightedConcentration_hi(1:nFluxWeightedConcentrationFile) = &
              workint1(1:nFluxWeightedConcentrationFile)
DEALLOCATE(workint1)

i = size(jzFluxWeightedConcentration_hi,1)
ALLOCATE(workint1(i))
workint1 = jzFluxWeightedConcentration_hi
DEALLOCATE(jzFluxWeightedConcentration_hi)
ALLOCATE(jzFluxWeightedConcentration_hi(nFluxWeightedConcentrationFile))
IF(nFluxWeightedConcentrationFile /= 0) jzFluxWeightedConcentration_hi(1:nFluxWeightedConcentrationFile) = &
              workint1(1:nFluxWeightedConcentrationFile)
DEALLOCATE(workint1)

!  ***************END OF OUTPUT********************************

!  *****************PEST BLOCK***********************

IF (ALLOCATED(PestScale)) THEN
DEALLOCATE(PestScale)
ALLOCATE(PestScale(ncomp))
ELSE
ALLOCATE(PestScale(ncomp))
END IF
PestScale = 1.0d0

section = 'pest'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN

CALL PestScaleOutput(nout)

END IF

WRITE(iunit2,*)
WRITE(iunit2,*) 'Parameters for this run:'
WRITE(*,*)
WRITE(*,*) 'Parameters for this run:'
WRITE(iunit2,*)
WRITE(*,*)
WRITE(iunit2,1011) tstep,delt
WRITE(iunit2,*)
IF (jpor == 1) THEN
WRITE(iunit2,*) '--> Porosity calculated from mineral volume fractions'
WRITE(iunit2,*) '----> Porosity is updated due to mineral dissolution and precipitation reactions'
WRITE(*,*)      '--> Porosity calculated from mineral volume fractions'
WRITE(*,*)      '----> Porosity is updated due to mineral dissolution and precipitation reactions'
ELSE IF (jpor == 2) THEN
WRITE(iunit2,*) '--> Porosity is read from file'
WRITE(iunit2,*) '----> No update of porosity'
WRITE(*,*)      '--> Porosity is read from file'
WRITE(*,*)      '----> No update of porosity'
ELSE IF (jpor == 3) THEN
WRITE(iunit2,*) '--> Porosity is read from file'
WRITE(iunit2,*) '----> With update of porosity after renormalization of mineral volume fractions'
WRITE(*,*)      '--> Porosity is read from file'
WRITE(*,*)      '--> With update of porosity after renormalization of mineral volume fractions'
ELSE IF (jpor == 0) THEN
WRITE(iunit2,*) '--> Porosity calculated from mineral volume fractions'
WRITE(iunit2,*) '----> No update of porosity'
WRITE(*,*)      '--> Porosity calculated from mineral volume fractions'
WRITE(*,*)      '----> No update of porosity'
ELSE
WRITE(iunit2,*) '--> Porosity set by "fix_porosity" or "set_porosity" keywords'
WRITE(iunit2,*) '----> No update of porosity'
WRITE(*,*)      '--> Porosity set by "fix_porosity" or "set_porosity" keywords'
WRITE(*,*)      '----> No update of porosity'
END IF
IF (igamma == 0) THEN
WRITE(iunit2,*) '--> Unit activity coefficients'
WRITE(*,*) '--> Unit activity coefficients'
ELSE IF (igamma == 2) THEN
WRITE(iunit2,*) '--> Extended Debye-Huckel activity model used'
WRITE(iunit2,*) '--> Dependence on activity coefficients'
WRITE(iunit2,*) '    NEGLECTED in Jacobian calculation'
WRITE(*,*) '--> Extended Debye-Huckel activity model used'
WRITE(*,*) '--> Dependence on activity coefficients'
WRITE(*,*) '    NEGLECTED in Jacobian calculation'
ELSE IF (igamma == 3) THEN
WRITE(iunit2,*) '--> Extended Debye-Huckel activity model used'
WRITE(iunit2,*) '--> Activity coefficients only computed'
WRITE(iunit2,*) '    at beginning of time step '
WRITE(*,*) '--> Extended Debye-Huckel activity model used'
WRITE(*,*) '--> Activity coefficients only computed'
WRITE(*,*) '    at beginning of time step '
ELSE
WRITE(iunit2,*) '--> Extended Debye-Huckel activity model used'
WRITE(iunit2,*) '--> Dependence on activity coefficients'
WRITE(iunit2,*) '    INCLUDED in Jacobian calculation'
WRITE(*,*) '--> Extended Debye-Huckel activity model used'
WRITE(*,*) '--> Dependence on activity coefficients'
WRITE(*,*) '    INCLUDED in Jacobian calculation'
END IF
WRITE(iunit2,*)
IF (jtemp == 0) THEN
tgradprt = 0.0
WRITE(iunit2,2002) tinit
WRITE(iunit2,2012) tgradprt
ELSE
tgradprt = tgrad
WRITE(iunit2,2002) tinit
WRITE(iunit2,2012) tgradprt
END IF
WRITE(iunit2,*)
WRITE(iunit2,1013)
WRITE(iunit2,1014) (prtint(i),i=1,nstop)

! ****************************************************************
IF (nxyz == 1) THEN
IF (ipath == 1) THEN
  WRITE(*,*)
  WRITE(*,*) '   Running in reaction path mode '
  WRITE(*,*) '   No update of mineral volume fractions '
  WRITE(*,*)
  WRITE(iunit2,*)
  WRITE(iunit2,*) '   Running in reaction path mode '
  WRITE(iunit2,*) '   No update of mineral volume fractions '
  WRITE(iunit2,*)
ELSE IF (ipath == 0) THEN
  WRITE(*,*)
  WRITE(*,*) '  Running in batch mode'
  WRITE(*,*)
  WRITE(iunit2,*)
  WRITE(iunit2,*) '  Running in batch mode'
  WRITE(iunit2,*)
ELSE
  WRITE(*,*)
  WRITE(*,*) '  Must specify either batch or reaction path'
  WRITE(*,*) '  mode when NXYZ = 1 (set ipath = 1 or 0)'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
END IF

!  *****************MODFLOW BLOCK***********************

nwells = 0
ncnh = 0
nrivers = 0
ndrains = 0

MODFLOWfile = ' '

WRITE(*,*) ' Reading MODFLOW block'

section = 'modflow'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
WRITE(*,*) ' MODFLOW block found'

CALL read_MODFLOWfile(nout,lfile,mxwell,mxrivr,mxdrn)

IF (.NOT. modflow .OR. MODFLOWfile == ' ') THEN
  WRITE(*,*)
  WRITE(*,*) ' Modflow file not specified '
  WRITE(*,*) ' MODFLOW keyword block should not be included without a Modflow *.hff file'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF


IF (ALLOCATED(qcnh)) THEN
  DEALLOCATE(qcnh)
  ALLOCATE(qcnh(ndimdummy))  !  Read this one from the *.hff file
ELSE
  ALLOCATE(qcnh(ndimdummy))  !  Read this one from the *.hff file
END IF
IF (ALLOCATED(q)) THEN
  DEALLOCATE(q)
  ALLOCATE(q(mxwell))
ELSE
  ALLOCATE(q(mxwell))
END IF
IF (ALLOCATED(qriver)) THEN
  DEALLOCATE(qriver)
  ALLOCATE(qriver(mxrivr))
ELSE
  ALLOCATE(qriver(mxrivr))
END IF
IF (ALLOCATED(qdrain)) THEN
  DEALLOCATE(qdrain)
  ALLOCATE(qdrain(mxdrn))
ELSE
  ALLOCATE(qdrain(mxdrn))
END IF
IF (ALLOCATED(cnhIn)) THEN
  DEALLOCATE(cnhIn)
  ALLOCATE(cnhIn(ndimdummy))
ELSE
  ALLOCATE(cnhIn(ndimdummy))
END IF
IF (ALLOCATED(wellIn)) THEN
  DEALLOCATE(wellIn)
  ALLOCATE(wellIn(mxwell))
ELSE
  ALLOCATE(wellIn(mxwell))
END IF
IF (ALLOCATED(riverIn)) THEN
  DEALLOCATE(riverIn)
  ALLOCATE(riverIn(mxrivr))
ELSE
  ALLOCATE(riverIn(mxrivr))
END IF
IF (ALLOCATED(jxWellLoc)) THEN
  DEALLOCATE(jxWellLoc)
  ALLOCATE(jxWellLoc(mxwell))
ELSE
  ALLOCATE(jxWellLoc(mxwell))
END IF
IF (ALLOCATED(jyWellLoc)) THEN
  DEALLOCATE(jyWellLoc)
  ALLOCATE(jyWellLoc(mxwell))
ELSE
  ALLOCATE(jyWellLoc(mxwell))
END IF
IF (ALLOCATED(jzWellLoc)) THEN
  DEALLOCATE(jzWellLoc)
  ALLOCATE(jzWellLoc(mxwell))
ELSE
  ALLOCATE(jzWellLoc(mxwell))
END IF
IF (ALLOCATED(jxRiverLoc)) THEN
  DEALLOCATE(jxRiverLoc)
  ALLOCATE(jxRiverLoc(mxrivr))
ELSE
  ALLOCATE(jxRiverLoc(mxrivr))
END IF
IF (ALLOCATED(jyRiverLoc)) THEN
  DEALLOCATE(jyRiverLoc)
  ALLOCATE(jyRiverLoc(mxrivr))
ELSE
  ALLOCATE(jyRiverLoc(mxrivr))
END IF
IF (ALLOCATED(jzRiverLoc)) THEN
  DEALLOCATE(jzRiverLoc)
  ALLOCATE(jzRiverLoc(mxrivr))
ELSE
  ALLOCATE(jzRiverLoc(mxrivr))
END IF
IF (ALLOCATED(jxDrainLoc)) THEN
  DEALLOCATE(jxDrainLoc)
  ALLOCATE(jxDrainLoc(mxdrn))
ELSE
  ALLOCATE(jxDrainLoc(mxdrn))
END IF
IF (ALLOCATED(jyDrainLoc)) THEN
  DEALLOCATE(jyDrainLoc)
  ALLOCATE(jyDrainLoc(mxdrn))
ELSE
  ALLOCATE(jyDrainLoc(mxdrn))
END IF
IF (ALLOCATED(jzDrainLoc)) THEN
  DEALLOCATE(jzDrainLoc)
  ALLOCATE(jzDrainLoc(mxdrn))
ELSE
  ALLOCATE(jzDrainLoc(mxdrn))
END IF
IF (ALLOCATED(jxHeadLoc)) THEN
  DEALLOCATE(jxHeadLoc)
  ALLOCATE(jxHeadLoc(ndimdummy))
ELSE
  ALLOCATE(jxHeadLoc(ndimdummy))
END IF
IF (ALLOCATED(jyHeadLoc)) THEN
  DEALLOCATE(jyHeadLoc)
  ALLOCATE(jyHeadLoc(ndimdummy))
ELSE
  ALLOCATE(jyHeadLoc(ndimdummy))
END IF
IF (ALLOCATED(jzHeadLoc)) THEN
  DEALLOCATE(jzHeadLoc)
  ALLOCATE(jzHeadLoc(ndimdummy))
ELSE
  ALLOCATE(jzHeadLoc(ndimdummy))
END IF

call ModScan(nx,ny,nz,cnhIn,wellIn,riverIn,ndimdummy,mxwell,mxrivr,mxdrn)

REWIND 1

IF (nwells /= mxwell) THEN
  WRITE(*,*)
  WRITE(*,*) ' "Nwells" read from "hff" file should match MXWELL from MODFLOW "wel" file'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

IF (ndrains /= mxdrn) THEN
  WRITE(*,*)
  WRITE(*,*) ' "Ndrains" read from "hff" file should match MXDRN from MODFLOW "drn" file'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

IF (nrivers/= mxrivr) THEN
  WRITE(*,*)
  WRITE(*,*) ' "Nrivers" read from "hff" file should match MXRIVR from MODFLOW "riv" file'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

! Check to see how many geochemical conditions are needed

NeedWellCondition = 0
DO i = 1,nwells
  IF (wellIn(I)) THEN
    NeedWellCondition = NeedWellCondition + 1
  END IF
END DO

NeedRiverCondition = 0
DO i = 1,nrivers
  IF (riverIn(I)) THEN
    NeedRiverCondition = NeedRiverCondition + 1
  END IF
END DO

NeedHeadCondition = 0
DO i = 1,ncnh
  IF (cnhIn(I)) THEN
    NeedHeadCondition = NeedHeadCondition + 1
  END IF
END DO

IF (ALLOCATED(qcnh)) THEN
  DEALLOCATE(qcnh)
END IF
ALLOCATE(qcnh(ncnh))


! read in dimensions
READ(52,*)
READ(52,*)
READ(52,*) nztemp,nytemp,nxtemp,nstress,modflowTimeUnits
READ(52,*)
READ(52,*)
READ(52,*)
IF (nztemp /= nz) THEN
  WRITE(*,*)
  WRITE(*,*) ' NZ does not match in MODFLOW and CRUNCH'
  WRITE(*,*) ' NZ in MODFLOW:  ', nztemp
  WRITE(*,*) ' NZ in CRUNCH:   ', nz
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
IF (nytemp /= ny) THEN
  WRITE(*,*)
  WRITE(*,*) ' NY does not match in MODFLOW and CRUNCH'
  WRITE(*,*) ' NY in MODFLOW:  ', nytemp
  WRITE(*,*) ' NY in CRUNCH:   ', ny
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
IF (nxtemp /= nx) THEN
  WRITE(*,*)
  WRITE(*,*) ' NX does not match in MODFLOW and CRUNCH'
  WRITE(*,*) ' NX in MODFLOW:  ', nxtemp
  WRITE(*,*) ' NX in CRUNCH:   ', nx
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

! Calculate conversion from MODFLOW time units to years

IF (ModFlowTimeUnits == 1) THEN        !  MODFLOW time units in seconds
  ModFlowCnv = 365.0*24.0*60.0*60.0
ELSE IF (ModFlowTimeUnits == 2) THEN   !  MODFLOW time units in minutes
  ModFlowCnv = 365.0*24.0*60.0
ELSE IF (ModFlowTimeUnits == 3) THEN   !  MODFLOW time units in hours
  ModFlowCnv = 365.0*24.0
ELSE IF (ModFlowTimeUnits == 4) THEN   !  MODFLOW time units in days
  ModFlowCnv = 365.0
ELSE IF (ModFlowTimeUnits == 5) THEN   !  MODFLOW time units in years
  ModFlowCnv = 1.0
ELSE
  WRITE(*,*)
  WRITE(*,*) ' MODFLOW time units read from *.bas file not recognized'
  WRITE(*,*) ' MODFLOW time units flag: ',ModFlowTimeUnits
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

! Read in cell types from MODFLOW .bas file

!fp! auto_par_loops=0;
DO jz=1,nz
  DO jy=1,ny
    READ (52,*) ActiveCell(1:nx,jy,jz)
  END DO
END DO
!fp! auto_par_loops=1;

READ(52,*) realjunk
READ(52,*)

ALLOCATE(work3(nx,ny,nz))

!fp! auto_par_loops=0;
DO jz=1,nz
  DO jy=1,ny
    READ (52,*) work3(1:nx,jy,jz)
  END DO
END DO
!fp! auto_par_loops=1;

DEALLOCATE(work3)

ALLOCATE(perlen(nstress))
ALLOCATE(tsmult(nstress))
ALLOCATE(nstp(nstress))

DO i = 1,nstress
  READ(52,*,err=3000) perlen(i),nstp(i),tsmult(i)
END DO

CLOSE(52)

ALLOCATE(jxTemp(ndimdummy))
ALLOCATE(jyTemp(ndimdummy))
ALLOCATE(jzTemp(ndimdummy))
ALLOCATE(conditionNum(ndimdummy))

!  Read wells

nparams = 0
jxTemp = 0
jyTemp = 0
jzTemp = 0
modflowstring = 'well'
CALL stringlen(modflowstring,lenstring)

CALL readModFlowParameters(nout,nchem,nparams,jxTemp,jyTemp,jzTemp,  &
            conditionNum,modflowstring,lenstring)

IF (NeedWellCondition > nparams) THEN
  WRITE(*,*)
  WRITE(*,*) ' Number of wells needing geochemical conditions does not match CRUNCH input file'
  WRITE(*,*) ' Number of conditions needed:   ', NeedWellCondition
  WRITE(*,*) ' Number of conditions provided: ', nparams
  WRITE(*,*) ' Need well conditions at the following locations: '
  WRITE(*,*)
  DO i = 1,nwells
    IF (wellIn(i)) THEN
      WRITE(*,1501) jxWellLoc(i), jyWellLoc(i), jzWellLoc(i)
    END IF
  END DO
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

1501 FORMAT('well  ', i4,1x,i4,1x,i4)
1502 FORMAT('river  ', i4,1x,i4,1x,i4)
1503 FORMAT('constant_head  ', i4,1x,i4,1x,i4)

IF (ALLOCATED(WellCondition)) THEN
  DEALLOCATE(WellCondition)
END IF
ALLOCATE(WellCondition(nwells))

! Spin through to check again

DO i = 1,nwells                !  Read from MODFLOW *.wel file
  DO j = 1,nparams             !  Read from CRUNCH input file
     IF (jxWellLoc(i)==jxTemp(j) .AND.                  &
         jyWellLoc(i)==jyTemp(j) .AND.                  &
         jzWellLoc(i)==jzTemp(j)       )  THEN
        WellCondition(i) = conditionNum(j)
       GO TO 501
     END IF
  END DO
  IF (wellIn(i)) THEN          !  Based on info in *.hff file, source requires a geochemical condition
    WRITE(*,*)
    WRITE(*,*) ' Well needs a geochemical condition--not found in CRUNCH input file'
    WRITE(*,*) ' Well number: ',i
    WRITE(*,*) ' Well location: ',jxWellLoc(i),jyWellLoc(i),jzWellLoc(i)
    WRITE(*,*)
    READ(*,*)
    STOP
  ELSE
    CONTINUE
  END IF
  501 CONTINUE
END DO


!  Read drains

nparams = 0
jxTemp = 0
jyTemp = 0
jzTemp = 0
modflowstring = 'drain'
CALL stringlen(modflowstring,lenstring)

CALL readModFlowParameters(nout,nchem,nparams,jxTemp,jyTemp,jzTemp,  &
            conditionNum,modflowstring,lenstring)

!  Read rivers

nparams = 0
jxTemp = 0
jyTemp = 0
jzTemp = 0
modflowstring = 'river'
CALL stringlen(modflowstring,lenstring)

CALL readModFlowParameters(nout,nchem,nparams,jxTemp,jyTemp,jzTemp,  &
            conditionNum,modflowstring,lenstring)


IF (NeedRiverCondition > nparams) THEN
  WRITE(*,*)
  WRITE(*,*) ' Number of rivers needing geochemical conditions does not match CRUNCH input file'
  WRITE(*,*) ' Number of conditions needed:   ', NeedRiverCondition
  WRITE(*,*) ' Number of conditions provided: ', nparams
  WRITE(*,*) ' Need river conditions at the following locations: '
  WRITE(*,*)
  DO i = 1,nrivers
    IF (riverIn(i)) THEN
      WRITE(*,1502) jxRiverLoc(i), jyRiverLoc(i), jzRiverLoc(i)
    END IF
  END DO
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

IF (ALLOCATED(RiverCondition)) THEN
  DEALLOCATE(RiverCondition)
END IF
ALLOCATE(RiverCondition(nrivers))

DO i = 1,nrivers                !  Read from MODFLOW *.wel file
  DO j = 1,nparams             !  Read from CRUNCH input file
     IF (jxRiverLoc(i)==jxTemp(j) .AND.                  &
         jyRiverLoc(i)==jyTemp(j) .AND.                  &
         jzRiverLoc(i)==jzTemp(j)       )  THEN
        RiverCondition(i) = conditionNum(j)
       GO TO 502
     END IF
  END DO
  IF (riverIn(i)) THEN          !  Based on info in *.hff file, source requires a geochemical condition
    WRITE(*,*)
    WRITE(*,*) ' River needs a geochemical condition--not found in CRUNCH input file'
    WRITE(*,*) ' River number: ',i
    WRITE(*,*) ' River location: ',jxRiverLoc(i),jyRiverLoc(i),jzRiverLoc(i)
    WRITE(*,*)
    READ(*,*)
    STOP
  ELSE
    CONTINUE
  END IF
  502 CONTINUE
END DO

!  Read constant head

nparams = 0
jxTemp = 0
jyTemp = 0
jzTemp = 0
modflowstring = 'constant_head'
CALL stringlen(modflowstring,lenstring)

CALL readModFlowParameters(nout,nchem,nparams,jxTemp,jyTemp,jzTemp,  &
            conditionNum,modflowstring,lenstring)

IF (NeedHeadCondition > nparams) THEN
  WRITE(*,*)
  WRITE(*,*) ' Number of constant heads needing geochemical conditions does not match CRUNCH input file'
  WRITE(*,*) ' Number of conditions needed:   ', NeedHeadCondition
  WRITE(*,*) ' Number of conditions provided: ', nparams
  WRITE(*,*) ' Need constant head conditions at the following locations: '
  WRITE(*,*)
  DO i = 1,ncnh
    IF (cnhIn(i)) THEN
      WRITE(*,1503) jxHeadLoc(i), jyHeadLoc(i), jzHeadLoc(i)
    END IF
  END DO
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

IF (ALLOCATED(HeadCondition)) THEN
  DEALLOCATE(HeadCondition)
END IF
ALLOCATE(HeadCondition(ncnh))

DO i = 1,ncnh                  !  Read from MODFLOW *.wel file
  DO j = 1,nparams             !  Read from CRUNCH input file
     IF (jxHeadLoc(i)==jxTemp(j) .AND.                  &
         jyHeadLoc(i)==jyTemp(j) .AND.                  &
         jzHeadLoc(i)==jzTemp(j)       )  THEN
        HeadCondition(i) = conditionNum(j)
       GO TO 503
     END IF
  END DO
  IF (cnhIn(i)) THEN          !  Based on info in *.hff file, source requires a geochemical condition
    WRITE(*,*)
    WRITE(*,*) ' Head needs a geochemical condition--not found in CRUNCH input file'
    WRITE(*,*) ' Head number: ',i
    WRITE(*,*) ' Head location: ',jxHeadLoc(i),jyHeadLoc(i),jzHeadLoc(i)
    WRITE(*,*)
    READ(*,*)
    STOP
  ELSE
    CONTINUE
  END IF
  503 CONTINUE
END DO

DEALLOCATE(jxTemp)
DEALLOCATE(jyTemp)
DEALLOCATE(jzTemp)
DEALLOCATE(conditionNum)

! **** Now read the geochemical condition for MODFLOW recharge

CALL readRecharge(nout,RechargeCondition)

! *****************************************

!  IF (RechargeCondition == 0) THEN
!    WRITE(*,*)
!    WRITE(*,*) ' Geochemical condtion for recharge not found in MODFLOW block'
!    WRITE(*,*)
!    STOP
!  END IF

!  Allocate the various arrays needed for the MODFLOW option

IF (ALLOCATED(qevt)) THEN
  DEALLOCATE(qevt)
END IF
ALLOCATE(qevt(nx,ny))
IF (ALLOCATED(wtype)) THEN
  DEALLOCATE(wtype)
END IF
ALLOCATE(wtype(nwells))
IF (ALLOCATED(htype)) THEN
  DEALLOCATE(htype)
END IF
ALLOCATE(htype(ncnh))
IF (ALLOCATED(rtype)) THEN
  DEALLOCATE(rtype)
END IF
ALLOCATE(rtype(nrivers))
IF (ALLOCATED(lrecharge)) THEN
  DEALLOCATE(lrecharge)
END IF
ALLOCATE(lrecharge(nx,ny))
IF (ALLOCATED(levt)) THEN
  DEALLOCATE(levt)
END IF
ALLOCATE(levt(nx,ny))

!  Now, initialize the suckers...

!!  qrecharge = 0.0
q = 0.0
qcnh = 0.0
qevt = 0.0
qriver = 0.0
qdrain = 0.0

wtype = ' '
htype = ' '
rtype = ' '


ALLOCATE(dtTemp(1000))

CALL readModFlowStress(nout,perlen,nstp,tsmult,ndimdummy)

ntt = 0
DO i = 1,nstress
  IF (nstp(i) == 1) THEN
    ntt = ntt + 1
    dtTemp(ntt) = perlen(i)
  ELSE IF (nstp(i) > 1) THEN
    geometricSum = 0.0
    DO k = 1,nstp(i)
      geometricSum = geometricSum + (tsmult(i)**(k-1))
    END DO
    dtBase = perlen(i)/geometricSum
    DO k = 1,nstp(i)
      ntt = ntt + 1
      dtTemp(ntt) = dtBase*(tsmult(i)**(k-1))
    END DO
    checkSum = 0.0
    DO k = 1,nstp(i)
      checkSum = checkSum + dtTemp(i)
    END DO
    IF (checkSum < (perlen(i)-eps) .OR. checksum > (perlen(i)+eps) ) THEN
      WRITE(*,*)
      WRITE(*,*) ' Modflow time steps calculated incorrectly'
      WRITE(*,*) ' Length of stress period:                  ',perlen(i)
      WRITE(*,*) ' Calculated length from sum of time steps: ',checkSum
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF

  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Number of time steps within a stress period cannot be less than 1'
    WRITE(*,*) ' Stress period number: ',i
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
END DO

IF (ALLOCATED(dtModFlow)) THEN
  DEALLOCATE(dtModFlow)
END IF
ALLOCATE(dtModFlow(ntt))

dtModFlow(1:ntt) = dtTemp(1:ntt)
DEALLOCATE(dtTemp)

ntModFlow = ntt

DEALLOCATE(nstp)
DEALLOCATE(perlen)
DEALLOCATE(tsmult)

END IF


!  *****************END OF MODFLOW BLOCK***********************


!  *****************FLOW BLOCK***********************

!  Initialize flow variables

readvelocity = .FALSE.
readgasvelocity = .FALSE.
WRITE(*,*) ' Reading flow block'

section = 'flow'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
  WRITE(*,*) ' Flow block found'
END IF

velocityfile = ' '
gasvelocityfile = ' '



IF (ALLOCATED(qrecharge)) THEN
DEALLOCATE(qrecharge)
END IF
ALLOCATE(qrecharge(nx,ny))
qrecharge = 0

IF (found) THEN

!  Initialize pressure and pumping rate first


IF (ALLOCATED(intbndgas)) THEN
  DEALLOCATE(intbndgas)
END IF
ALLOCATE(intbndgas(5,nx,ny,nz))
intbndgas = 0

IF (ALLOCATED(gaspump)) THEN
  DEALLOCATE(gaspump)
  ALLOCATE(gaspump(5,nx,ny,nz))
ELSE
  ALLOCATE(gaspump(5,nx,ny,nz))

END IF

IF (isaturate == 1) THEN
  gaspump = 0.0d0
END IF

CALL units_time(nout,section,time_scale)
CALL units_distance(nout,section,dist_scale)

CALL read_constantflow(nout,nx,ny,nz,constant_flow,qxinit,qyinit,qzinit)
CALL read_constantgasflow(nout,nx,ny,nz,constant_gasflow,  &
  qxgasinit,qygasinit,qzgasinit)


  pumptimeseries = .FALSE.

  !CALL read_pumptimeseriesfile(nout,nx,ny,nz,pumptimeseriesfile,lfile,pumptimeseries,PumptimeseriesFileFormat)
  !CALL read_pumplocationsfile(nout,nx,ny,nz,pumplocationsfile,lfile2,PumplocationsFileFormat)
  !IF (pumptimeseries) THEN
  !
  !CALL  read_pump_timeseries2(nout,nx,ny,nz,nchem,lfile,pumptimeseriesfile,PumptimeseriesFileFormat,lfile2,pumplocationsfile,PumplocationsFileFormat)
  !
  !ELSE
  !  CALL read_pump(nout,nx,ny,nz,nchem)
  !ENDIF
  
  CALL read_pumptimeseries(nout,nx,ny,nz)

  IF (pumptimeseries) THEN

  CALL read_pumplocations(nout,nx,ny,nz,nchem)
  ELSE
    CALL read_pump(nout,nx,ny,nz,nchem)
  ENDIF

!!!  IF (isaturate == 1) THEN
  CALL read_gaspump(nout,nx,ny,nz,nchem,ngaspump)
!!!  END IF

irecharge = 0

IF (.NOT. modflow) THEN
  CALL read_infiltration(nout,nx,ny,nz)
END IF

!!  Convert pumping rate from liters/sec to m**3/yr
!!  qg = qg*secyr/1000.0d0                  !!  Converting from l/sec to m**3/yr

!!!  IF (isaturate == 1) THEN
  gaspump = gaspump*secyr/1000.0d0                  !!  Converting from l/sec to m**3/yr
!!!  END IF

!  NOTE: Specification of constant flow field overrides all other
!        instructions (e.g., file read or flow calculation)

IF (constant_flow) THEN
  WRITE(*,*)
  WRITE(*,*) ' Constant flow specified'
  WRITE(*,*)
ELSE

!  No constant flow field specified, so look for file read or
!     flow calculation

  readvelocity = .false.
  CALL read_flowfile(nout,nx,ny,nz,constant_flow,  &
      qxinit,qyinit,qzinit,velocityfile,lfile,VelocityFileFormat)

  IF (velocityfile /= ' ') THEN
    readvelocity = .TRUE.
    WRITE(*,*)
    WRITE(*,*) ' Velocities to be read from file: ',velocityfile(1:lfile)
    WRITE(*,*)

  END IF

!  *******************FLOW CALCULATION****************************

!  Check to see if flow should be calculated--only do so if either constant_flow = .FALSE., or if
!  readvelocity = .FALSE.

  CalculateFlow = .FALSE.
  IF (.NOT. ReadVelocity) THEN
    parchar = 'calculate_flow'
    parfind = ' '
    CALL read_logical(nout,lchar,parchar,parfind,CalculateFlow)
  END IF

  flow_if: IF (CalculateFlow .AND. nxyz > 1) THEN

      ! Select NS solver or Darcy solver, added by Zhi Li
      NavierStokes = .FALSE.
      parchar = 'ns_solve'
      parfind = ' '
      CALL read_logical(nout,lchar,parchar,parfind,NavierStokes)
      ! ***************************************************

      ! Select Richards solver by Toshiyuki Bandai, 2023 May
      Richards = .FALSE.
      parchar = 'Richards'
      parfind = ' '
      CALL read_logical(nout,lchar,parchar,parfind,Richards)
      !

      IF (Richards) THEN
        isaturate = 1
      ENDIF

      ! True if you want print statement on the screen from the Richards solver
      Richards_print = .FALSE.
      parchar = 'Richards_print'
      parfind = ' '
      CALL read_logical(nout,lchar,parchar,parfind,Richards_print)
      
      ! True if steady-state Richards solver is used to obtain the initial condition
      Richards_steady = .FALSE.
      parchar = 'Richards_steady'
      parfind = ' '
      CALL read_logical(nout,lchar,parchar,parfind,Richards_steady)
      
      ! True if the n parameter in the van Genuchten model is used as the input data
      vg_is_n = .TRUE.
      parchar = 'vg_is_n'
      parfind = ' '
      CALL read_logical(nout,lchar,parchar,parfind,vg_is_n)
      
      ! True if the primary variable psi in the Richards equation is pressure head [L] or not. If false, the input values for the initial and boundary conditions, and vg_alpha are interpreted as in terms of pressure [Pa].  
      psi_is_head = .TRUE.
      parchar = 'psi_is_head'
      parfind = ' '
      CALL read_logical(nout,lchar,parchar,parfind,psi_is_head)
      
      IF (.NOT. psi_is_head .AND. ABS(dist_scale - 1.0d0) > 1.0d-5) THEN
        WRITE(*,*) 'dist_scale must be set to 1.0 when psi_is_head = .FALSE.'
        STOP
      END IF
      
      ! True if the theta_s is the same as the porosity value
      theta_s_is_porosity = .TRUE.
      parchar = 'theta_s_is_porosity'
      parfind = ' '
      CALL read_logical(nout,lchar,parchar,parfind,theta_s_is_porosity)
      
      ! True if the input to the theta_r is the residual saturation value
      theta_r_is_S_r = .FALSE.
      parchar = 'theta_r_is_S_r'
      parfind = ' '
      CALL read_logical(nout,lchar,parchar,parfind,theta_r_is_S_r)
      
      ! End of Edit by Toshiyuki Bandai, 2023 May
      ! ***************************************************
    
    IF (ALLOCATED(harx)) THEN
      DEALLOCATE(harx)
      ALLOCATE(harx(0:nx,1:ny,1:nz))
    ELSE
      ALLOCATE(harx(0:nx,1:ny,1:nz))
    END IF
    IF (ALLOCATED(hary)) THEN
      DEALLOCATE(hary)
      ALLOCATE(hary(1:nx,0:ny,1:nz))
    ELSE
      ALLOCATE(hary(1:nx,0:ny,1:nz))
    END IF
    IF (ALLOCATED(harz)) THEN
      DEALLOCATE(harz)
      ALLOCATE(harz(1:nx,1:ny,0:nz))
    ELSE
      ALLOCATE(harz(1:nx,1:ny,0:nz))
    END IF

    IF (ALLOCATED(perminx)) THEN
      DEALLOCATE(perminx)
      ALLOCATE(perminx(0:nx+1,1:ny,1:nz))
    ELSE
      ALLOCATE(perminx(0:nx+1,1:ny,1:nz))
    END IF

    IF (ALLOCATED(perminy)) THEN
      DEALLOCATE(perminy)
      ALLOCATE(perminy(1:nx,0:ny+1,1:nz))
    ELSE
      ALLOCATE(perminy(1:nx,0:ny+1,1:nz))
    END IF

    IF (ALLOCATED(perminz)) THEN
      DEALLOCATE(perminz)
      ALLOCATE(perminz(1:nx,1:ny,0:nz+1))
    ELSE
      ALLOCATE(perminz(1:nx,1:ny,0:nz+1))
    END IF

    IF (ALLOCATED(permx)) THEN
      DEALLOCATE(permx)
      ALLOCATE(permx(0:nx+1,1:ny,1:nz))
    ELSE
      ALLOCATE(permx(0:nx+1,1:ny,1:nz))
    END IF

    IF (ALLOCATED(permxOld)) THEN
      DEALLOCATE(permxOld)
      ALLOCATE(permxOld(0:nx+1,1:ny,1:nz))
    ELSE
      ALLOCATE(permxOld(0:nx+1,1:ny,1:nz))
    END IF

    IF (ALLOCATED(permy)) THEN
      DEALLOCATE(permy)
      ALLOCATE(permy(1:nx,0:ny+1,1:nz))
    ELSE
      ALLOCATE(permy(1:nx,0:ny+1,1:nz))
    END IF

    IF (ALLOCATED(permyOld)) THEN
      DEALLOCATE(permyOld)
      ALLOCATE(permyOld(1:nx,0:ny+1,1:nz))
    ELSE
      ALLOCATE(permyOld(1:nx,0:ny+1,1:nz))
    END IF

    IF (ALLOCATED(permz)) THEN
      DEALLOCATE(permz)
      ALLOCATE(permz(1:nx,1:ny,0:nz+1))
    ELSE
      ALLOCATE(permz(1:nx,1:ny,0:nz+1))
    END IF

    IF (ALLOCATED(permzOld)) THEN
      DEALLOCATE(permzOld)
      ALLOCATE(permzOld(1:nx,1:ny,0:nz+1))
    ELSE
      ALLOCATE(permzOld(1:nx,1:ny,0:nz+1))
    END IF

    pres = 0.0
    perminx = 0.0
    perminy = 0.0
    perminz = 0.0

    WRITE(*,*)
    WRITE(*,*) ' Flow will be calculated'
    WRITE(*,*)
    
   
  ! ********************************************
  ! Edit by Toshiyuki Bandai 2023 May
  Toshi_allocate: IF (Richards) THEN
  ! allocate state variable for the Richards equation
  ! water potential psi
    IF (ALLOCATED(psi)) THEN
      DEALLOCATE(psi)
      ALLOCATE(psi(nx, ny, nz))
    ELSE
      ALLOCATE(psi(nx, ny, nz))
    END IF
     
    ! volumetric water content theta
    IF (ALLOCATED(theta)) THEN
      DEALLOCATE(theta)
      ALLOCATE(theta(nx, ny, nz))
    ELSE
      ALLOCATE(theta(nx, ny, nz))
    END IF
    
    IF (ALLOCATED(theta_prev)) THEN
      DEALLOCATE(theta_prev)
      ALLOCATE(theta_prev(nx, ny, nz))
    ELSE
      ALLOCATE(theta_prev(nx, ny, nz))
    END IF

    ! derivative of volumetric water content theta
    IF (ALLOCATED(dtheta)) THEN
      DEALLOCATE(dtheta)
      ALLOCATE(dtheta(nx, ny, nz))
    ELSE
      ALLOCATE(dtheta(nx, ny, nz))
    END IF

    ! head values
    IF (ALLOCATED(head)) THEN
      DEALLOCATE(head)
      ALLOCATE(head(nx, ny, nz))
    ELSE
      ALLOCATE(head(nx, ny, nz))
    END IF

    ! relative permeability
    IF (ALLOCATED(kr)) THEN
      DEALLOCATE(kr)
      ALLOCATE(kr(nx, ny, nz))
    ELSE
      ALLOCATE(kr(nx, ny, nz))
    END IF

    ! derivative of relative permeability
    IF (ALLOCATED(dkr)) THEN
      DEALLOCATE(dkr)
      ALLOCATE(dkr(nx, ny, nz))
    ELSE
      ALLOCATE(dkr(nx, ny, nz))
    END IF
     
    ! allocate and read van-Genuchten parameters
    VG_error = 0
    
  ! saturated water content theta_s (=porosity)
    
    IF (theta_s_is_porosity) THEN
      ! theta_s is the same as the porosity, so no need to read vg_theta_s
      IF (ALLOCATED(theta_s)) THEN
        DEALLOCATE(theta_s)
        ALLOCATE(theta_s(nx, ny, nz))
      ELSE
        ALLOCATE(theta_s(nx, ny, nz))
      END IF
      
      FORALL (jx=1:nx, jy=1:ny, jz=1:nz)
        theta_s(jx,jy,jz) = por(jx,jy,jz)
      END FORALL
    ELSE
      parchar = 'vg_theta_s'
      CALL read_vanGenuchten_parameters(nout, lchar, parchar, parfind, section, nx, ny, nz, VG_error)
      IF (VG_error == 1) THEN
        WRITE(*,*)
        WRITE(*,*) ' Error in reading van Genuchten parameters for ', parchar
        WRITE(*,*)
        STOP
      END IF
    
    ! if theta_s does not match porosity, make a warning
      DO jx = 1, nx
        DO jy = 1, ny
          DO jz = 1, nz
            IF (theta_s(jx,jy,jz) /= por(jx,jy,jz)) THEN
              WRITE(*,*)
              WRITE(*,*) ' Warning: theta_s /= porosity at (',jx,',',jy,',',jz,')'
              WRITE(*,*) ' theta_s = ',theta_s(jx,jy,jz)
              WRITE(*,*) ' porosity = ',por(jx,jy,jz)
              WRITE(*,*)
            END IF
          END DO
        END DO
      END DO
    END IF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! residual water content (theta_r)
  !!!!!!!!!!!!!!!!!!!!!!!!!
    parfind = ' '
    parchar = 'vg_theta_r'
    CALL read_vanGenuchten_parameters(nout, lchar, parchar, parfind, section, nx, ny, nz, VG_error)
    IF (VG_error == 1) THEN
      WRITE(*,*)
      WRITE(*,*) ' Error in reading van Genuchten parameters for ', parchar
      WRITE(*,*)
      STOP
    END IF

    IF (parfind == ' ') THEN
      ! ************************************
      ! Edit by Lucien Stolze, June 2023
      ! Read wcr parameter array from file
      CALL read_wcrfile(nout,nx,ny,nz,wcrfile,lfile,readwcr,wcrFileFormat)
      IF (readwcr) then
        wcrfile(1:lfile) = wcrfile(1:lfile)
        INQUIRE(FILE=wcrfile,EXIST=ext)
        IF (.NOT. ext) THEN
          WRITE(*,*)
          WRITE(*,*) ' wcr file not found: ',wcrfile(1:lfile)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        IF (ALLOCATED(theta_r)) THEN
          DEALLOCATE(theta_r)
          ALLOCATE(theta_r(nx, ny, nz))
        ELSE
          ALLOCATE(theta_r(nx, ny, nz))
        END IF
        OPEN(UNIT=23,FILE=wcrfile,STATUS='old')
        FileTemp = wcrfile
        CALL stringlen(FileTemp,FileNameLength)
        IF (wcrFileFormat == 'ContinuousRead') THEN
          READ(23,*,END=1020) (((theta_r(jx,jy,jz),jx=0,nx+1),jy=1,ny),jz=1,nz)
        ELSE IF (wcrFileFormat == 'SingleColumn') THEN
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) theta_r(jx,jy,jz)
              END DO
            END DO
          END DO
        ELSE IF (wcrFileFormat == 'FullForm') THEN
          IF (ny > 1 .AND. nz > 1) THEN
            DO jz = 1,nz
              DO jy = 1,ny
                DO jx= 1,nx
                  READ(23,*,END=1020) xdum,ydum,zdum,theta_r(jx,jy,jz)
                END DO
              END DO
            END DO
          ELSE IF (ny > 1 .AND. nz == 1) THEN
            jz = 1
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) xdum,ydum,theta_r(jx,jy,jz)
              END DO
            END DO
          ELSE
          jz = 1
          jy = 1
          DO jx= 1,nx
            READ(23,*,END=1020) xdum,theta_r(jx,jy,jz)
          END DO
          END IF
        ELSE IF (wcrFileFormat == 'Unformatted') THEN
        READ(23,END=1020) theta_r
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' wcr file format not recognized'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      ELSE
      WRITE(*,*)
      WRITE(*,*) 'Information on residual water content must be provided.'
      WRITE(*,*)
      READ(*,*)
      STOP
    ENDIF

    ENDIF

  ! the input value is actually residual saturation, convert it to residual water content
    IF (theta_r_is_S_r) THEN
      theta_r = theta_r*theta_s
    END IF


  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! alpha parameter in the van Genuchten model
  !!!!!!!!!!!!!!!!!!!!!!!!!
    parfind = ' '
    parchar = 'vg_alpha'
    CALL read_vanGenuchten_parameters(nout, lchar, parchar, parfind, section, nx, ny, nz, VG_error)
    IF (VG_error == 1) THEN
      WRITE(*,*)
      WRITE(*,*) ' Error in reading van Genuchten parameters for ', parchar
      WRITE(*,*)
      STOP
    END IF

    IF (parfind == ' ') THEN
    ! ************************************
    ! Edit by Lucien Stolze, June 2023
    ! Read alpha parameter array from file
    CALL read_vgafile(nout,nx,ny,nz,vgafile,lfile,readvga,vgaFileFormat)
    IF (readvga) then
      vgafile(1:lfile) = vgafile(1:lfile)
      INQUIRE(FILE=vgafile,EXIST=ext)
      IF (.NOT. ext) THEN
        WRITE(*,*)
        WRITE(*,*) ' vga file not found: ',vgafile(1:lfile)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
      IF (ALLOCATED(VG_alpha)) THEN
        DEALLOCATE(VG_alpha)
        ALLOCATE(VG_alpha(nx, ny, nz))
      ELSE
        ALLOCATE(VG_alpha(nx, ny, nz))
      END IF
      OPEN(UNIT=23,FILE=vgafile,STATUS='old')
      FileTemp = vgafile
      CALL stringlen(FileTemp,FileNameLength)
      IF (vgaFileFormat == 'ContinuousRead') THEN
        READ(23,*,END=1020) (((VG_alpha(jx,jy,jz),jx=0,nx+1),jy=1,ny),jz=1,nz)
      ELSE IF (vgaFileFormat == 'SingleColumn') THEN
        DO jz = 1,nz
          DO jy = 1,ny
            DO jx= 1,nx
              READ(23,*,END=1020) VG_alpha(jx,jy,jz)
            END DO
          END DO
        END DO
      ELSE IF (vgaFileFormat == 'FullForm') THEN
        IF (ny > 1 .AND. nz > 1) THEN
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) xdum,ydum,zdum,VG_alpha(jx,jy,jz)
              END DO
            END DO
          END DO
        ELSE IF (ny > 1 .AND. nz == 1) THEN
          jz = 1
          DO jy = 1,ny
            DO jx= 1,nx
              READ(23,*,END=1020) xdum,ydum,VG_alpha(jx,jy,jz)
            END DO
          END DO
        ELSE
        jz = 1
        jy = 1
        DO jx= 1,nx
          READ(23,*,END=1020) xdum,VG_alpha(jx,jy,jz)
        END DO
        END IF
      ELSE IF (vgaFileFormat == 'Unformatted') THEN
      READ(23,END=1020) VG_alpha
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' VGalpha file format not recognized'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
    ELSE
    WRITE(*,*)
    WRITE(*,*) 'Information on VG alpha must be provided.'
    WRITE(*,*)
    READ(*,*)
    STOP
  ENDIF

  ENDIF
    ! ************************************
    ! End edit by Lucien Stolze, June 2023


  ! convert unit
    VG_alpha = VG_alpha*dist_scale
    
    IF (.NOT. psi_is_head) THEN
      VG_alpha = VG_alpha*rho_water*9.80665d0
    END IF
    
  !!!!!!!!!!!!!!!!!!!!!!!!!
  ! n parameter in the van Genuchten model
  !!!!!!!!!!!!!!!!!!!!!!!!!
    parchar = 'vg_n'
    CALL read_vanGenuchten_parameters(nout, lchar, parchar, parfind, section, nx, ny, nz, VG_error)
    IF (VG_error == 1) THEN
      WRITE(*,*)
      WRITE(*,*) ' Error in reading van Genuchten parameters for ', parchar
      WRITE(*,*)
      STOP
    END IF

    IF (parfind == ' ') THEN
      ! ************************************
      ! Edit by Lucien Stolze, June 2023
      ! Read alpha parameter array from file
      CALL read_vgnfile(nout,nx,ny,nz,vgnfile,lfile,readvgn,vgnFileFormat)
      IF (readvgn) then
        vgnfile(1:lfile) = vgnfile(1:lfile)
        INQUIRE(FILE=vgnfile,EXIST=ext)
        IF (.NOT. ext) THEN
          WRITE(*,*)
          WRITE(*,*) ' vgn file not found: ',vgnfile(1:lfile)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        IF (ALLOCATED(VG_n)) THEN
          DEALLOCATE(VG_n)
          ALLOCATE(VG_n(nx, ny, nz))
        ELSE
          ALLOCATE(VG_n(nx, ny, nz))
        END IF
        OPEN(UNIT=23,FILE=vgnfile,STATUS='old')
        FileTemp = vgnfile
        CALL stringlen(FileTemp,FileNameLength)
        IF (vgnFileFormat == 'ContinuousRead') THEN
          READ(23,*,END=1020) (((VG_n(jx,jy,jz),jx=0,nx+1),jy=1,ny),jz=1,nz)
        ELSE IF (vgnFileFormat == 'SingleColumn') THEN
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) VG_n(jx,jy,jz)
              END DO
            END DO
          END DO
        ELSE IF (vgnFileFormat == 'FullForm') THEN
          IF (ny > 1 .AND. nz > 1) THEN
            DO jz = 1,nz
              DO jy = 1,ny
                DO jx= 1,nx
                  READ(23,*,END=1020) xdum,ydum,zdum,VG_n(jx,jy,jz)
                END DO
              END DO
            END DO
          ELSE IF (ny > 1 .AND. nz == 1) THEN
            jz = 1
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) xdum,ydum,VG_n(jx,jy,jz)
              END DO
            END DO
          ELSE
          jz = 1
          jy = 1
          DO jx= 1,nx
            READ(23,*,END=1020) xdum,VG_n(jx,jy,jz)
          END DO
          END IF
        ELSE IF (vgnFileFormat == 'Unformatted') THEN
        READ(23,END=1020) VG_n
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' VGn file format not recognized'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      ELSE
      WRITE(*,*)
      WRITE(*,*) 'Information on VG n must be provided.'
      WRITE(*,*)
      READ(*,*)
      STOP
    ENDIF

    ENDIF
      ! ************************************
      ! End edit by Lucien Stolze, June 2023
    
    IF (.NOT. vg_is_n) THEN
        ! the input value is actually the parameter m (m = 1 - 1/n), so convert it to n
        VG_n = 1.0d0/(1.0d0 - VG_n)
    END IF
    
  ! psi_s parameter in the modified van Genuchten model
    !parchar = 'vg_psi_s'
    !CALL read_vanGenuchten_parameters(nout, lchar, parchar, parfind, section, nx, ny, nz, VG_error)
    !IF (VG_error == 1) THEN
    !  WRITE(*,*)
    !  WRITE(*,*) ' Error in reading van Genuchten parameters for ', parchar
    !  WRITE(*,*)
    !  STOP
    !END IF
    !! convert unit
    !psi_s = psi_s/dist_scale
    
  END IF Toshi_allocate
  ! ***************************************************
  ! End of Edit by Toshiyuki Bandai 2023 May
     

!  Look for information on permeability, pressure, and pumping or injection wells
!  First, check to see whether permeability distribution is to be read from file

    readperm = .false.
    npermx = 0
    npermy = 0
    npermz = 0
    CALL read_permfile(nout,nx,ny,nz,permfile,lfile,readperm,PermFileFormat)

    IF (readperm) THEN

      permxfile = ' '
      permyfile = ' '
      permzfile = ' '

      WRITE(iunit2,*)
      WRITE(iunit2,*) '  Reading permeabilities from file: ',permfile(1:lfile)
      WRITE(iunit2,*)

!!!!!!!!!!
      IF (PermFileFormat == 'SingleFile3D') THEN

        INQUIRE(FILE=permfile,EXIST=ext)
        IF (.NOT. ext) THEN
          WRITE(*,*)
          WRITE(*,*) ' 3D permeability file not found: ',permfile(1:lfile)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF

        OPEN(UNIT=23,FILE=permfile,STATUS='old',ERR=8005)
        FileTemp = permfile
        CALL stringlen(FileTemp,FileNameLength)
        DO jz = 1,nz
          DO jy = 1,ny
            DO jx = 1,nx
              READ(23,*,END=1020) xdum,ydum,zdum,permx(jx,jy,jz),permdum,permy(jx,jy,jz)
              permx(jx,jy,jz) = 0.001*permx(jx,jy,jz)/(9.81*997.075)
              permy(jx,jy,jz) = 0.001*permy(jx,jy,jz)/(9.81*997.075)
!!               permz(jx,jy,jz) = permz(jx,jy,jz)/(9.81*997.075)
            END DO
          END DO
        END DO

        jz = 1
        do jy = 1,ny
          permx(0,jy,jz) = permx(1,jy,jz)
          permx(nx+1,jy,jz) = permx(nx,jy,jz)
        end do

        perminx = permx
        perminy = permy
        perminx = permx

        CLOSE(UNIT=23,STATUS='keep')
!!!!!!!
    ELSE

        permxfile(1:lfile) = permfile(1:lfile)
        permxfile(lfile+1:lfile+2) = '.x'
        INQUIRE(FILE=permxfile,EXIST=ext)
        IF (.NOT. ext) THEN
          WRITE(*,*)
          WRITE(*,*) ' X permeability file not found: ',permxfile(1:lfile+2)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        OPEN(UNIT=23,FILE=permxfile,STATUS='old',ERR=8001)
        FileTemp = permxfile
        CALL stringlen(FileTemp,FileNameLength)
        IF (PermFileFormat == 'ContinuousRead') THEN
          READ(23,*,END=1020) (((permx(jx,jy,jz),jx=0,nx+1),jy=1,ny),jz=1,nz)
        ELSE IF (PermFileFormat == 'SingleColumn') THEN
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx= 0,nx+1
                READ(23,*,END=1020) permx(jx,jy,jz)
              END DO
            END DO
          END DO
        ELSE IF (PermFileFormat == 'FullForm') THEN
          IF (ny > 1 .AND. nz > 1) THEN
            DO jz = 1,nz
              DO jy = 1,ny
                DO jx= 0,nx+1
                  READ(23,*,END=1020) xdum,ydum,zdum,permx(jx,jy,jz)
                END DO
              END DO
            END DO
          ELSE IF (ny > 1 .AND. nz == 1) THEN
            jz = 1
            DO jy = 1,ny
              DO jx= 0,nx+1
                READ(23,*,END=1020) xdum,ydum,permx(jx,jy,jz)
              END DO
            END DO
          ELSE
          jz = 1
          jy = 1
          DO jx= 0,nx+1
            READ(23,*,END=1020) xdum,permx(jx,jy,jz)
          END DO
        END IF
      ELSE IF (PermFileFormat == 'Unformatted') THEN
        READ(23,END=1020) permx
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' X Permeability file format not recognized'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

      CLOSE(UNIT=23,STATUS='keep')

      termA = 2.14074E-05

      do jy = 1,ny
          do jx = 1,nx
                     porTemp = 0.15
                     portemp1 = 0.15

              termB=permx(jx,jy,1)/termA
              term1 = portemp1**03 - termB*portemp1*portemp1 + 2.0*termB*portemp1 - termB
              term2 = 3.0*portemp1**02 - 2.0*termB*portemp1 + 2.0*termB
              portemp1 = portemp1 - term1/term2


                              term1 = portemp1**03 - termB*portemp1*portemp1 + 2.0*termB*portemp1 - termB
              term2 = 3.0*portemp1**02 - 2.0*termB*portemp1 + 2.0*termB
              portemp1 = portemp1 - term1/term2



                            term1 = portemp1**03 - termB*portemp1*portemp1 + 2.0*termB*portemp1 - termB
              term2 = 3.0*portemp1**02 - 2.0*termB*portemp1 + 2.0*termB
              portemp1 = portemp1 - term1/term2



                          term1 = portemp1**03 - termB*portemp1*portemp1 + 2.0*termB*portemp1 - termB
              term2 = 3.0*portemp1**02 - 2.0*termB*portemp1 + 2.0*termB
              portemp1 = portemp1 - term1/term2



              term1 = portemp1**03 - termB*portemp1*portemp1 + 2.0*termB*portemp1 - termB
              term2 = 3.0*portemp1**02 - 2.0*termB*portemp1 + 2.0*termB
              portemp1 = portemp1 - term1/term2



              term1 = portemp1**03 - termB*portemp1*portemp1 + 2.0*termB*portemp1 - termB
              term2 = 3.0*portemp1**02 - 2.0*termB*portemp1 + 2.0*termB
              portemp1 = portemp1 - term1/term2



              term1 = portemp1**03 - termB*portemp1*portemp1 + 2.0*termB*portemp1 - termB
              term2 = 3.0*portemp1**02 - 2.0*termB*portemp1 + 2.0*termB
              portemp1 = portemp1 - term1/term2

             por(jx,jy,1) = portemp1
             volfx(2,jx,jy,1) = 1.0 - portemp1 - 0.05

          end do
      end do

      IF (ny /= 1) THEN
        permyfile(1:lfile) = permfile(1:lfile)
        permyfile(lfile+1:lfile+2) = '.y'
        INQUIRE(FILE=permyfile,EXIST=ext)
        IF (.NOT. ext) THEN
          WRITE(*,*)
          WRITE(*,*) ' Y permeability file not found: ',permyfile(1:lfile+2)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        OPEN(UNIT=23,FILE=permyfile,STATUS='old',ERR=8002)
        FileTemp = permyfile
        CALL stringlen(FileTemp,FileNameLength)
        IF (PermFileFormat == 'ContinuousRead') THEN
          READ(23,*,END=1020) (((permy(jx,jy,jz),jx=1,nx),jy=0,ny+1),jz=1,nz)
        ELSE IF (PermFileFormat == 'SingleColumn') THEN
          DO jz = 1,nz
            DO jy = 0,ny+1
              DO jx= 1,nx
                READ(23,*,END=1020) permy(jx,jy,jz)
              END DO
            END DO
          END DO
        ELSE IF (PermFileFormat == 'FullForm') THEN
          IF (ny > 1 .AND. nz > 1) THEN
            DO jz = 1,nz
              DO jy = 0,ny+1
                DO jx= 1,nx
                  READ(23,*,END=1020) xdum,ydum,zdum,permy(jx,jy,jz)
                END DO
              END DO
            END DO
          ELSE IF (ny > 1 .AND. nz == 1) THEN
            jz = 1
            DO jy = 0,ny+1
              DO jx= 1,nx
                READ(23,*,END=1020) xdum,ydum,permy(jx,jy,jz)
              END DO
            END DO
          ELSE
            CONTINUE
          END IF
        ELSE IF (PermFileFormat == 'Unformatted') THEN
          READ(23,END=1020) permy
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Y Permeability file format not recognized'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        CLOSE(UNIT=23,STATUS='keep')
      END IF

      IF (nz /= 1) THEN
        permzfile(1:lfile) = permfile(1:lfile)
        permzfile(lfile+1:lfile+2) = '.z'
        INQUIRE(FILE=permzfile,EXIST=ext)
        IF (.NOT. ext) THEN
          WRITE(*,*)
          WRITE(*,*) ' Z permeability file not found: ',permzfile(1:lfile+2)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        OPEN(UNIT=23,FILE=permzfile,STATUS='old',ERR=8003)
        FileTemp = permzfile
        CALL stringlen(FileTemp,FileNameLength)
        IF (PermFileFormat == 'ContinuousRead') THEN
          READ(23,*,END=1020) (((permz(jx,jy,jz),jx=1,nx),jy=1,ny),jz=0,nz+1)
        ELSE IF (PermFileFormat == 'SingleColumn') THEN
          DO jz = 0,nz+1
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) permz(jx,jy,jz)
              END DO
            END DO
          END DO
        ELSE IF (PermFileFormat == 'FullForm') THEN
          IF (ny > 1 .AND. nz > 1) THEN
            DO jz = 0,nz+1
              DO jy = 1,ny
                DO jx= 1,nx
                  READ(23,*,END=1020) xdum,ydum,zdum,permz(jx,jy,jz)
                END DO
              END DO
            END DO
          ELSE
            CONTINUE
          END IF
        ELSE IF (PermFileFormat == 'Unformatted') THEN
          READ(23,END=1020) permz
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Z Permeability file format not recognized'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        CLOSE(UNIT=23,STATUS='keep')
      END IF

      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!

      permmaxX = 0.0d0
      permmaxy = 0.0d0
      permminX = 1.0d0
      permminY = 1.0d0

      do jz = 1,nz
        do jy = 1,ny
          do jx = 1,nx

            if (permx(jx,jy,jz) > permmaxX) then
              permmaxX = permx(jx,jy,jz)
            end if
            if (permY(jx,jy,jz) > permmaxY) then
              permmaxY = permY(jx,jy,jz)
            end if
            if (permx(jx,jy,jz) < permminX) then
              permminX = permx(jx,jy,jz)
            end if
            if (permy(jx,jy,jz) < permminy) then
              permminY = permy(jx,jy,jz)
            end if

          end do
        end do
      end do

      write(*,*) ' Internal MinMax in X: ', permminX,permmaxX
      write(*,*) ' Internal MinMax in Y: ', permminY,permmaxY
      write(*,*)


!!!        do jz = 1,nz
!!!          do jy = 1,ny
!!!            do jx = 1,nx
!!!
!!!                LogPermX = DLOG10(permx(jx,jy,jz))
!!!                NewLogPermX = (LogPermX - MidPointX)/2.0d0 + MidPointX
!!!                permx(jx,jy,jz) = 10**(NewLogPermX)
!!!
!!!                LogPermY = DLOG10(permY(jx,jy,jz))
!!!                NewLogPermY = (LogPermY - MidPointY)/2.0d0 + MidPointY
!!!                permy(jx,jy,jz) = 10**(NewLogPermY)
!!!
!!!           end do
!!!         end do
!!!        end do
!!!
!!!        permmaxX = 0.0d0
!!!        permmaxy = 0.0d0
!!!        permminX = 1.0d0
!!!        permminY = 1.0d0

!!!        do jz = 1,nz
!!!          do jy = 1,ny
!!!            do jx = 1,nx
!!!
!!!              if (permx(jx,jy,jz) > permmaxX) then
!!!                permmaxX = permx(jx,jy,jz)
!!!              end if
!!!              if (permY(jx,jy,jz) > permmaxY) then
!!!                permmaxY = permy(jx,jy,jz)
!!!              end if
!!!              if (permx(jx,jy,jz) < permminX) then
!!!                permminX = permx(jx,jy,jz)
!!!              end if
!!!              if (permy(jx,jy,jz) < permminY) then
!!!                permminY = permy(jx,jy,jz)
!!!                if (permminY == 0.0d0) then
!!!                   write(*,*) jx,jy
!!!                   read(*,*)
!!!                end if
!!!              end if
!!!
!!!            end do
!!!          end do
!!!        end do
!!!
!!!
!!!        write(*,*) ' MinMax in X: ', permminX,permmaxX
!!!        write(*,*) ' MinMax in Y: ', permminY,permmaxY
!!!        write(*,*)

      perminx = permx
      perminy = permy
      perminz = permz

      GO TO 8004

      8001   WRITE(*,*) ' Error opening X permeability file: STOP'
      READ(*,*)
      STOP
      8002   WRITE(*,*) ' Error opening Y permeability file: STOP'
      READ(*,*)
      STOP
      8003   WRITE(*,*) ' Error opening Z permeability file: STOP'
      READ(*,*)
      8005   WRITE(*,*) ' Error opening 3D permeability file: STOP'
      READ(*,*)

8004    CONTINUE

      CONTINUE

    END IF

!!!!      ELSE

      ALLOCATE(permzonex(0:mperm))
      ALLOCATE(permzoney(0:mperm))
      ALLOCATE(permzonez(0:mperm))

      ALLOCATE(jxxpermx_lo(mperm))
      ALLOCATE(jxxpermx_hi(mperm))
      ALLOCATE(jyypermx_lo(mperm))
      ALLOCATE(jyypermx_hi(mperm))
      ALLOCATE(jzzpermx_lo(mperm))
      ALLOCATE(jzzpermx_hi(mperm))

      ALLOCATE(jxxpermy_lo(mperm))
      ALLOCATE(jxxpermy_hi(mperm))
      ALLOCATE(jyypermy_lo(mperm))
      ALLOCATE(jyypermy_hi(mperm))
      ALLOCATE(jzzpermy_lo(mperm))
      ALLOCATE(jzzpermy_hi(mperm))

      ALLOCATE(jxxpermz_lo(mperm))
      ALLOCATE(jxxpermz_hi(mperm))
      ALLOCATE(jyypermz_lo(mperm))
      ALLOCATE(jyypermz_hi(mperm))
      ALLOCATE(jzzpermz_lo(mperm))
      ALLOCATE(jzzpermz_hi(mperm))

!  Reading permeability distribution directly from input file


      CALL read_permx(nout,nx,ny,nz,npermx)
      readpermx = .false.
       CALL read_permxfile(nout,nx,ny,nz,permxfile,lfile,readpermx,permxFileFormat)

       IF (readpermx) then

        permxfile(1:lfile) = permxfile(1:lfile)
        INQUIRE(FILE=permxfile,EXIST=ext)
        IF (.NOT. ext) THEN
          WRITE(*,*)
          WRITE(*,*) ' permx file not found: ',permxfile(1:lfile)
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        OPEN(UNIT=23,FILE=permxfile,STATUS='old')
        FileTemp = permxfile
        CALL stringlen(FileTemp,FileNameLength)

        IF (permxFileFormat == 'ContinuousRead') THEN
          READ(23,*,END=1020) (((VG_alpha(jx,jy,jz),jx=0,nx+1),jy=1,ny),jz=1,nz)
        ELSE IF (permxFileFormat == 'SingleColumn') THEN
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) permx(jx,jy,jz)
              END DO
            END DO
          END DO
        ELSE IF (permxFileFormat == 'FullForm') THEN
          IF (ny > 1 .AND. nz > 1) THEN
            DO jz = 1,nz
              DO jy = 1,ny
                DO jx= 1,nx
                  READ(23,*,END=1020) xdum,ydum,zdum,permx(jx,jy,jz)
                END DO
              END DO
            END DO
          ELSE IF (ny > 1 .AND. nz == 1) THEN
            jz = 1
            DO jy = 1,ny
              DO jx= 1,nx
                READ(23,*,END=1020) xdum,ydum,permx(jx,jy,jz)
              END DO
            END DO
          ELSE
          jz = 1
          jy = 1
          DO jx= 1,nx
            READ(23,*,END=1020) xdum,permx(jx,jy,jz)
          END DO
          END IF
        ELSE IF (permxFileFormat == 'Unformatted') THEN
        READ(23,END=1020) permx
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' VGalpha file format not recognized'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
      else
      IF (permzonex(0) == 0.0 .and. .NOT. ReadPerm) THEN
        WRITE(*,*)
        WRITE(*,*) ' No default X permeability given'
        WRITE(*,*) ' X permeability should be followed by "default" or blank string'
        WRITE(*,*) ' Or, turn off flow calculation (calculate_flow = false)'
        WRITE(*,*)
        STOP
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Background X permeability = ',permzonex(0)
        WRITE(*,*)
      END IF

! First, initialize X permeability to default permeability (permzonex(0))


      IF (.NOT. ReadPerm) THEN
        perminx = permzonex(0)
      END IF

!  Next, initialize permeability from various zones

      IF (npermx > 0) THEN
        DO l = 1,npermx
          DO jz = jzzpermx_lo(l),jzzpermx_hi(l)
            DO jy = jyypermx_lo(l),jyypermx_hi(l)
              DO jx = jxxpermx_lo(l),jxxpermx_hi(l)
                perminx(jx,jy,jz) = permzonex(l)
              END DO
            END DO
          END DO
        END DO
      END IF

      permmaxX= 0.0
      permx = perminx
      permmaxX = MAXVAL(DABS(permx))

      DEALLOCATE(permzonex)
      DEALLOCATE(jxxpermx_lo)
      DEALLOCATE(jxxpermx_hi)
      DEALLOCATE(jyypermx_lo)
      DEALLOCATE(jyypermx_hi)
      DEALLOCATE(jzzpermx_lo)
      DEALLOCATE(jzzpermx_hi)

    endif

      IF (ny == 1) THEN
        permy = 0.0
        perminy = 0.0
        DEALLOCATE(permzoney)
        DEALLOCATE(jxxpermy_lo)
        DEALLOCATE(jxxpermy_hi)
        DEALLOCATE(jyypermy_lo)
        DEALLOCATE(jyypermy_hi)
        DEALLOCATE(jzzpermy_lo)
        DEALLOCATE(jzzpermy_hi)
      ELSE

        readpermy = .false.
        CALL read_permyfile(nout,nx,ny,nz,permyfile,lfile,readpermy,permyFileFormat)

        IF (readpermy) then

         permyfile(1:lfile) = permyfile(1:lfile)
         INQUIRE(FILE=permyfile,EXIST=ext)
         IF (.NOT. ext) THEN
           WRITE(*,*)
           WRITE(*,*) ' permy file not found: ',permyfile(1:lfile)
           WRITE(*,*)
           READ(*,*)
           STOP
         END IF
         OPEN(UNIT=23,FILE=permyfile,STATUS='old')
         FileTemp = permyfile
         CALL stringlen(FileTemp,FileNameLength)
             jz = 1
             DO jy = 0,ny+1
               DO jx= 1,nx
                 READ(23,*,END=1020) xdum,ydum,dum1
                 if (dum1<-20) then
                permy(jx,jy,1:nz)=0
                 else
                 permy(jx,jy,1:nz)=10**dum1
                 endif
               END DO
             END DO
       else


        CALL read_permy(nout,nx,ny,nz,npermy)
        IF (permzoney(0) == 0.0 .and. .NOT.ReadPerm) THEN
          WRITE(*,*)
          WRITE(*,*) ' No default Y permeability given'
          WRITE(*,*) ' Y permeability should be followed by "default" or blank string'
          WRITE(*,*) ' Or, turn off flow calculation (calculate_flow = false)'
          WRITE(*,*)
          STOP
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Background Y permeability = ',permzoney(0)
          WRITE(*,*)
        END IF

! Initialize Y permeability to default permeability (permzoney(0))

!!        DO jy = 1,ny
!!          perminy(:,jy,:) = permzoney(0)
!!        END DO

      IF (.NOT. ReadPerm) THEN
        perminy = permzoney(0)
      END IF


!  Next, initialize permeability from various zones

        IF (npermy > 0) THEN
          DO l = 1,npermy
            DO jz = jzzpermy_lo(l),jzzpermy_hi(l)
              DO jy = jyypermy_lo(l),jyypermy_hi(l)
                DO jx = jxxpermy_lo(l),jxxpermy_hi(l)
                  perminy(jx,jy,jz) = permzoney(l)
                END DO
              END DO
            END DO
          END DO
        END IF

        permmaxy = 0.0
        permy = perminy
        permmaxy = MAXVAL(DABS(permy))

        DEALLOCATE(permzoney)
        DEALLOCATE(jxxpermy_lo)
        DEALLOCATE(jxxpermy_hi)
        DEALLOCATE(jyypermy_lo)
        DEALLOCATE(jyypermy_hi)
        DEALLOCATE(jzzpermy_lo)
        DEALLOCATE(jzzpermy_hi)

      END IF
      END IF

      IF (nz == 1) THEN

        IF (nmmLogical .and. FractureNetwork) THEN

          IF (CriticalZone) THEN

            do jy = 1,ny
              do jx = 1,nx

                if (jinit(jx,jy,1) == 2) then
                  perminx(jx,jy,1) = 1.0D-12
                  permx(jx,jy,1) = 1.0D-12
                  perminy(jx,jy,1) = 1.0D-12
                  permy(jx,jy,1) = 1.0D-12
                end if

              end do
            end do

            do jy = 1,2
              do jx = 1,nx          !!! Soil layer 2 grid cells deep

                  perminx(jx,jy,1) = 1.0D-12
                  permx(jx,jy,1) = 1.0D-12
                  perminy(jx,jy,1) = 1.0D-12
                  permy(jx,jy,1) = 1.0D-12
                  jinit(jx,jy,1) = 3

              end do
            end do

          ELSE

            do jy = 1,ny
              do jx = 1,nx

                if (jinit(jx,jy,1) == 2) then
                  perminx(jx,jy,1) = 1.0D-11
                  permx(jx,jy,1) = 1.0D-11
                  perminy(jx,jy,1) = 1.0D-11
                  permy(jx,jy,1) = 1.0D-11
                end if

              end do
            end do

          END IF


        END IF

        permz = 0.0
        perminz = 0.0
        DEALLOCATE(permzonez)
        DEALLOCATE(jxxpermz_lo)
        DEALLOCATE(jxxpermz_hi)
        DEALLOCATE(jyypermz_lo)
        DEALLOCATE(jyypermz_hi)
        DEALLOCATE(jzzpermz_lo)
        DEALLOCATE(jzzpermz_hi)

      ELSE


        CALL read_permz(nout,nx,ny,nz,npermz)
        IF (permzonez(0) == 0.0 .and. .NOT.ReadPerm) THEN
          WRITE(*,*)
          WRITE(*,*) ' No default Z permeability given'
          WRITE(*,*) ' Z permeability should be followed by "default" or blank string'
          WRITE(*,*) ' Or, turn off flow calculation (calculate_flow = false)'
          WRITE(*,*)
          STOP
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Background Z permeability = ',permzonez(0)
          WRITE(*,*)
        END IF

! Initialize Y permeability to default permeability (permzonez(0))

!!        DO jz = 1,nz
!!          perminz(:,:,jz) = permzonez(0)
!!        END DO

      IF (.NOT. ReadPerm) THEN
        perminz = permzonez(0)
        permz = permzonez(0)
      END IF


!  Next, initialize permeability from various zones

        IF (npermz > 0) THEN
          DO l = 1,npermz
            DO jz = jzzpermz_lo(l),jzzpermz_hi(l)
              DO jy = jyypermz_lo(l),jyypermz_hi(l)
                DO jx = jxxpermz_lo(l),jxxpermz_hi(l)
                  perminz(jx,jy,jz) = permzonez(l)
                END DO
              END DO
            END DO
          END DO
        END IF

        permmaxz = 0.0
        permz = perminz
        permmaxz = MAXVAL(DABS(permz))

        DEALLOCATE(permzonez)
        DEALLOCATE(jxxpermz_lo)
        DEALLOCATE(jxxpermz_hi)
        DEALLOCATE(jyypermz_lo)
        DEALLOCATE(jyypermz_hi)
        DEALLOCATE(jzzpermz_lo)
        DEALLOCATE(jzzpermz_hi)

      END IF
    
    ! ********************************************
    ! Edit by Toshiyuki Bandai, 2023 May
    Toshi_permeability: IF (Richards) THEN
      ! allocate permeability at faces
      IF (ALLOCATED(K_faces_x)) THEN
        DEALLOCATE(K_faces_x)
        ALLOCATE(K_faces_x(0:nx, ny, nz))
      ELSE
        ALLOCATE(K_faces_x(0:nx, ny, nz))
      END IF

      K_faces_x(0, 1, 1) = permx(1, 1, 1)
      K_faces_x(nx, 1, 1) = permx(nx, 1, 1)

      ! compute face permeability
      DO i = 1, nx - 1
        IF (ABS(permx(i, 1, 1) - permx(i + 1, 1, 1)) < 1.0d-20) THEN
          K_faces_x(i, 1, 1) = permx(i, 1, 1)
        ELSE
          numerator = permx(i, 1, 1) * permx(i+1, 1, 1) * (x(i+1) - x(i))
          denominator = 0.5d0 * permx(i, 1, 1) * dxx(i) + 0.5d0 * permx(i+1, 1, 1) * dxx(i+1)
          K_faces_x(i, 1, 1) = numerator / denominator
        END IF
      END DO
    END IF Toshi_permeability
    
    ! Read lower and upper boundary conditions for steady-state problem
    Toshi_boundary_conditions: IF (Richards) THEN
      IF (Richards_steady) THEN
        ! lower boundary condition
        BC_location = 0
        ! the arguments upper_BC_file, lfile, tslength are not used for the steady-state problem
        CALL read_boundary_condition_Richards(nout, Richards_steady, BC_location, lower_BC_type_steady, upper_BC_file, value_lower_BC_steady, lfile, upper_constant_BC_steady, tslength)
      
        ! unit conversion
        SELECT CASE (lower_BC_type_steady)
        CASE ('constant_dirichlet')
          value_lower_BC_steady = value_lower_BC_steady/dist_scale
          
          IF (.NOT. psi_is_head) THEN
            value_lower_BC_steady = (value_lower_BC_steady - pressure_air)/(rho_water*9.80665d0)
          END IF
          
        CASE ('constant_neumann')
          WRITE(*,*)
          WRITE(*,*) ' The upper boundary condition type ', upper_BC_type, ' is not supported for the upper boundary condition for the steady-state Richards equation. '
          WRITE(*,*)
          READ(*,*)
          STOP
          !CONTINUE ! no unit conversion
        CASE ('constant_flux')
          value_lower_BC_steady = value_lower_BC_steady/(dist_scale * time_scale)
        CASE DEFAULT
          WRITE(*,*)
          WRITE(*,*) ' The lower boundary condition type ', lower_BC_type_steady, ' is not supported for the steady-state Richards equation. '
          WRITE(*,*)
          READ(*,*)
          STOP
        END SELECT
      
        ! upper boundary condition
        BC_location = 1
        ! the arguments upper_BC_file, lfile, upper_constant_BC, tslength are not used for the steady-state problem
        CALL read_boundary_condition_Richards(nout, Richards_steady, BC_location, upper_BC_type_steady, upper_BC_file, value_upper_BC_steady, lfile, upper_constant_BC_steady, tslength)
      
        ! unit conversion
        SELECT CASE (upper_BC_type_steady)
        CASE ('constant_dirichlet')
          value_upper_BC_steady = value_upper_BC_steady/dist_scale
          IF (.NOT. psi_is_head) THEN
            value_upper_BC_steady = (value_upper_BC_steady - pressure_air)/(rho_water*9.80665d0)
          END IF
          
        CASE ('constant_neumann')
          CONTINUE ! no unit conversion
        CASE ('constant_flux')
          value_upper_BC_steady = value_upper_BC_steady/(dist_scale * time_scale)
        CASE DEFAULT
          WRITE(*,*)
          WRITE(*,*) ' The upper boundary condition type ', upper_BC_type_steady, ' is not supported for the steady-state Richards equation. '
          WRITE(*,*)
          READ(*,*)
          STOP
        END SELECT
      
      END IF
    
      ! read boundary conditions for transient problem
      ! read lower boundary condition
      BC_location = 0
      lower_constant_BC = .TRUE.
    
      CALL read_boundary_condition_Richards(nout, .FALSE., BC_location, lower_BC_type, lower_BC_file, value_lower_BC, lfile, lower_constant_BC, tslength)
        
      ! unit conversion and import time series data if the boundary condition is variable
      SELECT CASE (lower_BC_type)
      CASE ('constant_dirichlet')
        value_lower_BC = value_lower_BC/dist_scale
        IF (.NOT. psi_is_head) THEN
            value_lower_BC = (value_lower_BC - pressure_air)/(rho_water*9.80665d0)
        END IF
      CASE ('constant_neumann')
        CONTINUE ! no unit conversion
      CASE ('constant_flux')
        value_lower_BC = value_lower_BC/(dist_scale * time_scale)
      CASE ('variable_dirichlet')
      
        IF (ALLOCATED(t_lower_BC)) THEN
          DEALLOCATE(t_lower_BC)
        END IF
        IF (ALLOCATED(values_lower_BC)) THEN
          DEALLOCATE(values_lower_BC)
        END IF
        ALLOCATE(t_lower_BC(tslength))
        ALLOCATE(values_lower_BC(tslength))
        CALL read_timeseries(nout, nx, ny, nz, t_lower_BC, values_lower_BC, lfile, lower_BC_file, tslength)
      
        values_lower_BC = values_lower_BC/dist_scale
        
        IF (.NOT. psi_is_head) THEN
            values_lower_BC = (values_lower_BC - pressure_air)/(rho_water*9.80665d0)
        END IF
        
        ! unit conversion for the time for the variable boundary condition
        t_lower_BC = t_lower_BC*time_scale
      CASE ('variable_neumann')
      
        IF (ALLOCATED(t_lower_BC)) THEN
          DEALLOCATE(t_lower_BC)
        END IF
        IF (ALLOCATED(values_lower_BC)) THEN
          DEALLOCATE(values_lower_BC)
        END IF
        ALLOCATE(t_lower_BC(tslength))
        ALLOCATE(values_lower_BC(tslength))
        CALL read_timeseries(nout, nx, ny, nz, t_lower_BC, values_lower_BC, lfile, lower_BC_file, tslength)
      
        ! unit conversion for the time for the variable boundary condition
        t_lower_BC = t_lower_BC*time_scale
        CONTINUE ! no unit conversion
      CASE ('variable_flux')
      
        IF (ALLOCATED(t_lower_BC)) THEN
          DEALLOCATE(t_lower_BC)
        END IF
        IF (ALLOCATED(values_lower_BC)) THEN
          DEALLOCATE(values_lower_BC)
        END IF
        ALLOCATE(t_lower_BC(tslength))
        ALLOCATE(values_lower_BC(tslength))
        CALL read_timeseries(nout, nx, ny, nz, t_lower_BC, values_lower_BC, lfile, lower_BC_file, tslength)
      
        values_lower_BC = values_lower_BC/(dist_scale * time_scale)
        ! unit conversion for the time for the variable boundary condition
        t_lower_BC = t_lower_BC*time_scale
      CASE DEFAULT
        WRITE(*,*)
        WRITE(*,*) ' The lower boundary condition type ', lower_BC_type, ' is not supported. '
        WRITE(*,*)
        READ(*,*)
        STOP
      END SELECT
    
      ! read upper boundary condition
      upper_constant_BC = .TRUE.
      BC_location = 1
      CALL read_boundary_condition_Richards(nout, .FALSE., BC_location, upper_BC_type, upper_BC_file, value_upper_BC, lfile, upper_constant_BC, tslength)
      
      ! unit conversion and import time series for upper boundary condition if the boundary condition is time-dependent (variable)
      SELECT CASE (upper_BC_type)
      CASE ('constant_dirichlet')
        value_upper_BC = value_upper_BC/dist_scale
        IF (.NOT. psi_is_head) THEN
            value_upper_BC = (value_upper_BC - pressure_air)/(rho_water*9.80665d0)
        END IF
      CASE ('constant_neumann')
        WRITE(*,*)
        WRITE(*,*) ' The upper boundary condition type ', upper_BC_type, ' is not supported for the upper boundary condition. '
        WRITE(*,*)
        READ(*,*)
        STOP
        !CONTINUE ! no unit conversion
      CASE ('constant_flux')
        value_upper_BC = value_upper_BC/(dist_scale * time_scale)
      CASE ('variable_dirichlet')
        ! import time series for upper boundary condition
        IF (ALLOCATED(t_upper_BC)) THEN
          DEALLOCATE(t_upper_BC)
        END IF
        IF (ALLOCATED(values_upper_BC)) THEN
          DEALLOCATE(values_upper_BC)
        END IF
        ALLOCATE(t_upper_BC(tslength))
        ALLOCATE(values_upper_BC(tslength))
        CALL read_timeseries(nout, nx, ny, nz, t_upper_BC, values_upper_BC, lfile, upper_BC_file, tslength)
      
        values_upper_BC = values_upper_BC/dist_scale
        IF (.NOT. psi_is_head) THEN
            values_upper_BC = (values_upper_BC - pressure_air)/(rho_water*9.80665d0)
        END IF
        ! unit conversion for the time for the variable boundary condition
        t_upper_BC = t_upper_BC*time_scale
      CASE ('variable_neumann')
        WRITE(*,*)
        WRITE(*,*) ' The upper boundary condition type ', upper_BC_type, ' is not supported for the upper boundary condition. '
        WRITE(*,*)
        READ(*,*)
        STOP
        !CONTINUE ! no unit conversion
      CASE ('variable_flux')
        ! import time series for upper boundary condition
        IF (ALLOCATED(t_upper_BC)) THEN
          DEALLOCATE(t_upper_BC)
        END IF
        IF (ALLOCATED(values_upper_BC)) THEN
          DEALLOCATE(values_upper_BC)
        END IF
        ALLOCATE(t_upper_BC(tslength))
        ALLOCATE(values_upper_BC(tslength))
        CALL read_timeseries(nout, nx, ny, nz, t_upper_BC, values_upper_BC, lfile, upper_BC_file, tslength)
      
        values_upper_BC = values_upper_BC/(dist_scale * time_scale)
        ! unit conversion for the time for the variable boundary condition
        t_upper_BC = t_upper_BC*time_scale
    
      CASE ('environmental_forcing')
        ! import infiltration, evaporation, and transpiration data
      
        ! bool for whether the time series is provided or not
        infiltration_timeseries = .FALSE.
        evapotimeseries = .FALSE.
        transpitimeseries = .FALSE.
      
        ! bool for whether the forcing is constant or not
        infiltration_fix = .FALSE.
        evapofix = .FALSE.
        transpifix = .FALSE.
      
        ! infiltration
        CALL read_infiltration_2(nout,nx,ny,nz,infiltration_file,infiltration_rate,lfile,infiltration_fix,infiltration_timeseries,tslength)
        IF (infiltration_timeseries) THEN
          IF (ALLOCATED(t_infiltration)) THEN
            DEALLOCATE(t_infiltration)
          END IF
          IF (ALLOCATED(qt_infiltration)) THEN
            DEALLOCATE(qt_infiltration)
          END IF
          ALLOCATE(t_infiltration(tslength))
          ALLOCATE(qt_infiltration(tslength))
          CALL read_timeseries(nout,nx,ny,nz,t_infiltration,qt_infiltration,lfile,infiltration_file,tslength)
        
          ! unit conversion for the time for the variable boundary condition
          qt_infiltration = qt_infiltration/(dist_scale * time_scale)
          t_infiltration = t_infiltration*time_scale
        ENDIF
      
        IF (infiltration_fix) THEN
          infiltration_rate = infiltration_rate/(dist_scale * time_scale)
        END IF
      
        ! evaporation
        CALL read_evaporation(nout,nx,ny,nz,evapofile,evaporate,lfile,evapofix,evapotimeseries,tslength)
        IF (evapotimeseries) THEN
          IF (ALLOCATED(t_evapo)) THEN
            DEALLOCATE(t_evapo)
          END IF
          IF (ALLOCATED(qt_evapo)) THEN
            DEALLOCATE(qt_evapo)
          END IF
          ALLOCATE(t_evapo(tslength))
          ALLOCATE(qt_evapo(tslength))
          CALL read_timeseries(nout,nx,ny,nz,t_evapo,qt_evapo,lfile,evapofile,tslength)
        
          ! unit conversion for the time for the variable boundary condition
          qt_evapo = qt_evapo/(dist_scale * time_scale)
          t_evapo = t_evapo*time_scale
        
          ! evaporate = qt_evapo(1)
          ! STOP
        ENDIF
      
        IF (evapofix) THEN
          evaporate = evaporate/(dist_scale * time_scale)
        END IF


      
        ! transpiration
        IF (ALLOCATED(transpirate_cell)) THEN
          DEALLOCATE(transpirate_cell)
        END IF
        ALLOCATE(transpirate_cell(nx))
        transpirate_cell = 0
        
        CALL read_transpiration(nout,nx,ny,nz,transpifile,transpirate,lfile,transpifix,transpitimeseries,transpicells,tslength)
        IF (transpitimeseries) THEN
          IF (ALLOCATED(t_transpi)) THEN
            DEALLOCATE(t_transpi)
          END IF
          IF (ALLOCATED(qt_transpi)) THEN
            DEALLOCATE(qt_transpi)
          END IF
          ALLOCATE(t_transpi(tslength))
          ALLOCATE(qt_transpi(tslength))
          CALL read_timeseries(nout,nx,ny,nz,t_transpi,qt_transpi,lfile,transpifile,tslength)
          qt_transpi = qt_transpi/transpicells
        
          qt_transpi = qt_transpi/(dist_scale * time_scale)
          t_transpi = t_transpi*time_scale
          !transpirate = qt_transpi(1)
          !WRITE(*,*) transpirate
         ! STOP
        ENDIF
      
        IF (transpifix) THEN
          transpirate = transpirate / transpicells
          transpirate = transpirate/(dist_scale * time_scale)
          WRITE(*,*) 'stop'
        END IF

        
      
      CASE DEFAULT
        WRITE(*,*)
        WRITE(*,*) ' The upper boundary condition type ', upper_BC_type, ' is not supported. '
        WRITE(*,*)
        READ(*,*)
        STOP
      END SELECT
    
      ! Read initial condition for steady-state or transient problem
      parchar = 'read_richards_ic_file'
      parfind = ' '
      Richards_IC_File = ' '
      CALL readFileName(nout,lchar,parchar,parfind,dumstring,section,Richards_IC_FileFormat)
      IF (parfind == ' ') THEN
        WRITE(*,*) ' The initial condition file was not found. Set to zero water potential at all cells. '
        psi = 0.0d0
      ELSE
        Richards_IC_File = dumstring
        CALL stringlen(Richards_IC_File,ls)
        WRITE(*,*) ' Reading the initial condition for the Richards equation from file: ',Richards_IC_File(1:ls)
      END IF
    
    END IF Toshi_boundary_conditions
    
    Toshi_initial_conditions: IF (Richards) THEN
      read_ic_Rihcards: IF (Richards_IC_File /= ' ') THEN
        INQUIRE(FILE=Richards_IC_File,EXIST=ext)
      IF (.NOT. ext) THEN
        CALL stringlen(Richards_IC_File,ls)
        WRITE(*,*)
        WRITE(*,*) ' Initial condition file not found: ', Richards_IC_File(1:ls)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
        OPEN(UNIT=52,FILE=Richards_IC_File,STATUS='OLD',ERR=6001)
        FileTemp = Richards_IC_File
        CALL stringlen(FileTemp,FileNameLength)
        IF (Richards_IC_FileFormat == 'SingleColumn') THEN
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx= 1,nx
                READ(52,*,END=1020) psi(jx,jy,jz)
              END DO
            END DO
          END DO
          ! convert unit
          psi = psi / dist_scale
          
          ! the input value is in pressure [Pa], convert to pressure head [m]
          IF (.NOT. psi_is_head) THEN
            psi = (psi - pressure_air)/(rho_water*9.80665d0)
          END IF
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Richards Initial condition file format not recognized'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
        CLOSE(UNIT=52)
      END IF read_ic_Rihcards
    END IF Toshi_initial_conditions
    ! End of edit by Toshiyuki Bandai, 2023 May
    ! ********************************************
      
!!!  End of section where choosing between perm file read and zone specification      END IF

    CALL read_gravity(nout)

    IF (ALLOCATED(activecellPressure)) THEN
      DEALLOCATE(activecellPressure)
      ALLOCATE(activecellPressure(0:nx+1,0:ny+1,0:nz+1))
    ELSE
      ALLOCATE(activecellPressure(0:nx+1,0:ny+1,0:nz+1))
    END IF

    activecellPressure = 1
    
    CALL read_pressureAlternative(nout,nx,ny,nz,npressure)

    pres = PressureZone(0)

!  Next, initialize pressure from various zones
      IF (npressure > 0) THEN
        DO l = 1,npressure
          DO jz = jzzPressure_lo(l),jzzPressure_hi(l)
            DO jy = jyyPressure_lo(l),jyyPressure_hi(l)
              DO jx = jxxPressure_lo(l),jxxPressure_hi(l)
                pres(jx,jy,jz) = PressureZone(l)
                activecellPressure(jx,jy,jz) = PressureFix(l)
              END DO
            END DO
          END DO
        END DO
      END IF

    DEALLOCATE(PressureZone)
    DEALLOCATE(PressureFix)
    DEALLOCATE(jxxPressure_lo)
    DEALLOCATE(jxxPressure_hi)
    DEALLOCATE(jyyPressure_lo)
    DEALLOCATE(jyyPressure_hi)
    DEALLOCATE(jzzPressure_lo)
    DEALLOCATE(jzzPressure_hi)



    watertabletimeseries = .FALSE.
    CALL read_watertablefile(nout,nx,ny,nz,watertablefile,lfile,watertabletimeseries,WatertableFileFormat)
    IF (watertabletimeseries) THEN
    CALL  read_watertable_timeseries(nout,nx,ny,nz,lfile,watertablefile,WatertableFileFormat)
    !!else
    !!  CALL read_pump(nout,nx,ny,nz,nchem)
    ENDIF

    parchar = 'initialize_hydrostatic'
    parfind = ' '
    InitializeHydrostatic = .FALSE.
    CALL read_logical(nout,lchar,parchar,parfind,InitializeHydrostatic)
    IF (gimrt) THEN
      WRITE(*,*)
      WRITE(*,*) ' --> Initializing flow field to be hydrostatic '
      WRITE(*,*)
    ELSE
      CONTINUE
    END IF

  END IF flow_if  ! End of block within which flow calculation parameters are read

END IF

IF (constant_gasflow) THEN
  WRITE(*,*)
  WRITE(*,*) ' Constant gas flow specified'
  readgasvelocity = .FALSE.
  WRITE(*,*)
ELSE

!  No constant gas flow field specified, so look for file read

  readgasvelocity = .false.
  CALL read_gasflowfile(nout,nx,ny,nz,constant_gasflow,  &
      qxgasinit,qygasinit,qzgasinit,gasvelocityfile,lfile,GasVelocityFileFormat)

  IF (gasvelocityfile /= ' ') THEN
    readgasvelocity = .true.
!!      WRITE(*,*)
!!      WRITE(*,*) ' Gas velocities to be read from file ',gasvelocityfile(1:lfile)
!!      WRITE(*,*)

  END IF

END IF

IF (CalculateFlow) THEN
  
  IF (CalciteCreep) THEN
    DO jz = 1,nz
      DO jy = 1,ny
         DO jx = 1,nx
          IF (jinit(jx,jy,jz) == 1) THEN
            permX(jx,jy,jz) =   1.0E-10
            permY(jx,jy,jz) =   1.0E-10
            perminX(jx,jy,jz) = 1.0E-10
            perminY(jx,jy,jz) = 1.0E-10
 !!!           porin(jx,jy,jz) = 0.999
 !!!           por(jx,jy,jz) =   0.999
          ELSE IF (jinit(jx,jy,jz) == 3) THEN
            permX(jx,jy,jz) =   1.0E-13
            permY(jx,jy,jz) =   1.0E-13
            perminX(jx,jy,jz) = 1.0E-13
            perminY(jx,jy,jz) = 1.0E-13
 !!!           porin(jx,jy,jz) = 0.999
 !!!           por(jx,jy,jz) =   0.999
          ELSE
            permX(jx,jy,jz) = 1.0E-18
            permY(jx,jy,jz) = 1.0E-18
            perminX(jx,jy,jz) = 1.0E-18
            perminY(jx,jy,jz) = 1.0E-18
 !!!           porin(jx,jy,jz) = 0.0001
 !!!           por(jx,jy,jz)   = 0.0001
          ENDIF
        END DO
      END DO
    END DO
  ENDIF
  
  permxOld = permx
  permyOld = permy
  permzOld = permz
END IF

ELSE
WRITE(*,*)
WRITE(*,*) ' No flow block found'
readvelocity = .FALSE.
readgasvelocity = .FALSE.
WRITE(*,*) ' Assuming Darcy fluxes = 0'
WRITE(*,*)
qxinit = 0.0
qyinit = 0.0
qzinit = 0.0
qxgasinit = 0.0
qygasinit = 0.0
qzgasinit = 0.0
constant_flow = .TRUE.
constant_gasflow = .TRUE.
END IF

!  If constant_flow = true, then use constant velocities and skip file read

IF (constant_flow) THEN
WRITE(*,*)
WRITE(*,*) '  Constant velocities have been specified'
WRITE(*,*)
WRITE(*,*) '  X Darcy velocity = ',qxinit
WRITE(*,*) '  Y Darcy velocity = ',qyinit
WRITE(*,*) '  Z Darcy velocity = ',qzinit
WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,*) ' Constant velocities have been specified'
WRITE(iunit2,*)
WRITE(iunit2,*) '  X Darcy velocity = ',qxinit
WRITE(iunit2,*) '  Y Darcy velocity = ',qyinit
WRITE(iunit2,*)

!  Convert units if necessary (converting to years)

qxinit = qxinit/(time_scale*dist_scale)
qyinit = qyinit/(time_scale*dist_scale)
qzinit = qzinit/(time_scale*dist_scale)

qxmax = qxinit
qx = qxinit

qymax = qyinit
qy = qyinit

qzmax = qzinit
qz = qzinit

END IF

IF (constant_gasflow) THEN
WRITE(*,*)
WRITE(*,*) '  Constant gas velocities have been specified'
WRITE(*,*)
WRITE(*,*) '  X gas flux  = ',qxgasinit
WRITE(*,*) '  Y gas flux  = ',qygasinit
WRITE(*,*) '  Z gas flux  = ',qzgasinit
WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,*) ' Constant gas velocities have been specified'
WRITE(iunit2,*)
WRITE(iunit2,*) '  X gas flux  = ',qxgasinit
WRITE(iunit2,*) '  Y gas flux  = ',qygasinit
WRITE(iunit2,*) '  Z gas flux  = ',qzgasinit
WRITE(iunit2,*)

!  Convert units if necessary (converting to years)

qxgasinit = qxgasinit/(time_scale*dist_scale)
qygasinit = qygasinit/(time_scale*dist_scale)
qzgasinit = qzgasinit/(time_scale*dist_scale)

qxgas = qxgasinit
qygas = qygasinit
qzgas = qzgasinit

END IF


! Constant_flow = false, so read from file

IF (readvelocity) THEN
CALL stringlen(velocityfile,lfile)
WRITE(*,*)
WRITE(*,*) '  Reading velocities from file: ',velocityfile(1:lfile)
WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,*) '  Reading velocities from file: ',velocityfile(1:lfile)
WRITE(iunit2,*)

vxfile(1:lfile) = velocityfile(1:lfile)
vxfile(lfile+1:lfile+3) = '.vx'
INQUIRE(FILE=vxfile,EXIST=ext)
IF (.NOT. ext) THEN
  WRITE(*,*)
  WRITE(*,*) ' X velocity file not found: ',vxfile(1:lfile+3)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
OPEN(UNIT=23,FILE=vxfile,STATUS='old',ERR=5001)

FileTemp = vxfile
CALL stringlen(FileTemp,FileNameLength)
IF (VelocityFileFormat == 'ContinuousRead') THEN
  READ(23,*,END=1020) (((qx(jx,jy,jz),jx=0,nx),jy=1,ny),jz=1,nz)
ELSE IF (VelocityFileFormat == 'SingleColumn') THEN
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx= 0,nx
        READ(23,*,END=1020) qx(jx,jy,jz)
      END DO
    END DO
  END DO
ELSE IF (VelocityFileFormat == 'FullForm') THEN
!!    IF (ny > 1 .AND. nz > 1) THEN
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx= 0,nx
          READ(23,*,END=1020) xdum,ydum,zdum,qx(jx,jy,jz)
        END DO
      END DO
    END DO
!!    ELSE IF (ny > 1 .AND. nz == 1) THEN
!!      jz = 1
!!      DO jy = 1,ny
!!        DO jx= 0,nx
!!          READ(23,*,END=1020) xdum,ydum,qx(jx,jy,jz)
!!        END DO
!!      END DO
!!    ELSE
!!      jz = 1
!!      jy = 1
!!      DO jx= 0,nx
!!        READ(23,*,END=1020) xdum,qx(jx,jy,jz)
!!      END DO
!!    END IF
ELSE IF (VelocityFileFormat == 'Unformatted') THEN
  READ(23,END=1020) qx
ELSE
  WRITE(*,*)
  WRITE(*,*) ' X Velocity file format not recognized'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

CLOSE(UNIT=23,STATUS='keep')

IF (ny /= 1) THEN
  vyfile(1:lfile) = velocityfile(1:lfile)
  vyfile(lfile+1:lfile+3) = '.vy'
  INQUIRE(FILE=vyfile,EXIST=ext)
  IF (.NOT. ext) THEN
    WRITE(*,*)
    WRITE(*,*) ' Y velocity file not found: ',vyfile(1:lfile+3)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  OPEN(UNIT=23,FILE=vyfile,STATUS='old',ERR=5002)
  FileTemp = vyfile
  CALL stringlen(FileTemp,FileNameLength)
  IF (VelocityFileFormat == 'ContinuousRead') THEN
    READ(23,*,END=1020) (((qy(jx,jy,jz),jx=1,nx),jy=0,ny),jz=1,nz)
  ELSE IF (VelocityFileFormat == 'SingleColumn') THEN
    DO jz = 1,nz
      DO jy = 0,ny
        DO jx= 1,nx
          READ(23,*,END=1020) qy(jx,jy,jz)
        END DO
      END DO
    END DO
  ELSE IF (VelocityFileFormat == 'FullForm') THEN
    IF (ny > 1 .AND. nz > 1) THEN
      DO jz = 1,nz
        DO jy = 0,ny
          DO jx= 1,nx
            READ(23,*,END=1020) xdum,ydum,zdum,qy(jx,jy,jz)
          END DO
        END DO
      END DO
    ELSE IF (ny > 1 .AND. nz == 1) THEN
      jz = 1
      DO jy = 0,ny
        DO jx= 1,nx
          READ(23,*,END=1020) xdum,ydum,qy(jx,jy,jz)
        END DO
      END DO
    ELSE
      CONTINUE
    END IF
  ELSE IF (VelocityFileFormat == 'Unformatted') THEN
    READ(23,END=1020) qy
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Y Velocity file format not recognized'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  CLOSE(UNIT=23,STATUS='keep')
END IF     !!  End of IF block for reading Y velocities

IF (nz /= 1) THEN
  vzfile(1:lfile) = velocityfile(1:lfile)
  vzfile(lfile+1:lfile+3) = '.vz'
  INQUIRE(FILE=vzfile,EXIST=ext)
  IF (.NOT. ext) THEN
    WRITE(*,*)
    WRITE(*,*) ' Z velocity file not found: ',vzfile(1:lfile+3)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  OPEN(UNIT=23,FILE=vzfile,STATUS='old',ERR=5003)

  FileTemp = vzfile
  CALL stringlen(FileTemp,FileNameLength)
  IF (VelocityFileFormat == 'ContinuousRead') THEN
    READ(23,*,END=1020) (((qz(jx,jy,jz),jx=1,nx),jy=1,ny),jz=0,nz)
  ELSE IF (VelocityFileFormat == 'SingleColumn') THEN
    DO jz = 0,nz
      DO jy = 1,ny
        DO jx= 1,nx
          READ(23,*,END=1020) qz(jx,jy,jz)
        END DO
      END DO
    END DO
  ELSE IF (VelocityFileFormat == 'FullForm') THEN
    IF (ny > 1 .AND. nz > 1) THEN
      DO jz = 0,nz
        DO jy = 1,ny
          DO jx= 1,nx
            READ(23,*,END=1020) xdum,ydum,zdum,qz(jx,jy,jz)
          END DO
        END DO
      END DO
    ELSE
      CONTINUE
    END IF
  ELSE IF (VelocityFileFormat == 'Unformatted') THEN
    READ(23,END=1020) qz
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Z Velocity file format not recognized'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  CLOSE(UNIT=23,STATUS='keep')
END IF   !!  End of IF block for Z velocity

GO TO 5004

5001   WRITE(*,*) ' Error opening X velocity file: STOP'
READ(*,*)
STOP
5002   WRITE(*,*) ' Error opening Y velocity file: STOP'
READ(*,*)
STOP
5003   WRITE(*,*) ' Error opening Z velocity file: STOP'
READ(*,*)
STOP

5004   CONTINUE

qx = qx/(time_scale*dist_scale)
qxmax = MAXVAL(qx)

qy = qy/(time_scale*dist_scale)
qymax = MAXVAL(qy)

qz = qz/(time_scale*dist_scale)
qzmax = MAXVAL(qz)


END IF   ! End of block for file read of velocities


WRITE(*,*)
!!WRITE(*,*) '  Maximum X velocity (m/yr) = ',qxmax
!!WRITE(*,*) '  Maximum Y velocity (m/yr) = ',qymax
!!WRITE(*,*) '  Maximum Z velocity (m/yr) = ',qzmax
WRITE(*,*)
WRITE(iunit2,*)
!!WRITE(iunit2,*) '  Maximum X velocity (m/yr) = ',qxmax
!!WRITE(iunit2,*) '  Maximum Y velocity (m/yr) = ',qymax
!!WRITE(iunit2,*) '  Maximum Z velocity (m/yr) = ',qzmax
WRITE(iunit2,*)


! Constant_gasflow = false, so read from file

IF (readgasvelocity) THEN
call stringlen(gasvelocityfile,lfile)
WRITE(*,*)
WRITE(*,*) '  Reading gas velocities from file: ',gasvelocityfile(1:lfile)
WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,*) '  Reading gas velocities from file: ',gasvelocityfile(1:lfile)
WRITE(iunit2,*)

vxfile = ' '
vyfile = ' '
vzfile = ' '

vxfile(1:lfile) = gasvelocityfile(1:lfile)
vxfile(lfile+1:lfile+3) = '.vx'
INQUIRE(FILE=vxfile,EXIST=ext)
IF (.NOT. ext) THEN
  WRITE(*,*)
  WRITE(*,*) ' X gas velocity file not found: ',vxfile(1:lfile+3)
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
OPEN(UNIT=23,FILE=vxfile,STATUS='old',ERR=7001)

FileTemp = vxfile
CALL stringlen(FileTemp,FileNameLength)
IF (GasVelocityFileFormat == 'ContinuousRead') THEN
  READ(23,*,END=1020) (((qxgas(jx,jy,jz),jx=0,nx),jy=1,ny),jz=1,nz)
ELSE IF (GasVelocityFileFormat == 'SingleColumn') THEN
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx= 0,nx
        READ(23,*,END=1020) qxgas(jx,jy,jz)
      END DO
    END DO
  END DO
ELSE IF (GasVelocityFileFormat == 'FullForm') THEN
  IF (ny > 1 .AND. nz > 1) THEN
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx= 0,nx
          READ(23,*,END=1020) xdum,ydum,zdum,qxgas(jx,jy,jz)
        END DO
      END DO
    END DO
  ELSE IF (ny > 1 .AND. nz == 1) THEN
    jz = 1
    DO jy = 1,ny
      DO jx= 0,nx
        READ(23,*,END=1020) xdum,ydum,qxgas(jx,jy,jz)
      END DO
    END DO
  ELSE
    jz = 1
    jy = 1
    DO jx= 0,nx
      READ(23,*,END=1020) xdum,qxgas(jx,jy,jz)
    END DO
  END IF
ELSE IF (GasVelocityFileFormat == 'Unformatted') THEN
  READ(23,END=1020) qxgas
ELSE
  WRITE(*,*)
  WRITE(*,*) ' X Gas Velocity file format not recognized'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF

CLOSE(UNIT=23,STATUS='keep')

IF (ny /= 1) THEN
  vyfile(1:lfile) = gasvelocityfile(1:lfile)
  vyfile(lfile+1:lfile+3) = '.vy'
  INQUIRE(FILE=vyfile,EXIST=ext)
  IF (.NOT. ext) THEN
    WRITE(*,*)
    WRITE(*,*) ' Y gas velocity file not found: ',vyfile(1:lfile+3)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  OPEN(UNIT=23,FILE=vyfile,STATUS='old',ERR=7002)
  FileTemp = vyfile
  CALL stringlen(FileTemp,FileNameLength)
  IF (GasVelocityFileFormat == 'ContinuousRead') THEN
    READ(23,*,END=1020) (((qygas(jx,jy,jz),jx=1,nx),jy=0,ny),jz=1,nz)
  ELSE IF (GasVelocityFileFormat == 'SingleColumn') THEN
    DO jz = 1,nz
      DO jy = 0,ny
        DO jx= 1,nx
          READ(23,*,END=1020) qygas(jx,jy,jz)
        END DO
      END DO
    END DO
  ELSE IF (GasVelocityFileFormat == 'FullForm') THEN
    IF (ny > 1 .AND. nz > 1) THEN
      DO jz = 1,nz
        DO jy = 0,ny
          DO jx= 1,nx
            READ(23,*,END=1020) xdum,ydum,zdum,qygas(jx,jy,jz)
          END DO
        END DO
      END DO
    ELSE IF (ny > 1 .AND. nz == 1) THEN
      jz = 1
      DO jy = 0,ny
        DO jx= 1,nx
          READ(23,*,END=1020) xdum,ydum,qygas(jx,jy,jz)
        END DO
      END DO
    ELSE
      CONTINUE
    END IF
  ELSE IF (GasVelocityFileFormat == 'Unformatted') THEN
    READ(23,END=1020) qygas
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Y Velocity file format not recognized'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  CLOSE(UNIT=23,STATUS='keep')
END IF     !!  End of IF block for reading Y velocities

IF (nz /= 1) THEN
  vzfile(1:lfile) = gasvelocityfile(1:lfile)
  vzfile(lfile+1:lfile+3) = '.vz'
  INQUIRE(FILE=vzfile,EXIST=ext)
  IF (.NOT. ext) THEN
    WRITE(*,*)
    WRITE(*,*) ' Z velocity file not found: ',vzfile(1:lfile+3)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  OPEN(UNIT=23,FILE=vzfile,STATUS='old',ERR=7003)

  FileTemp = vzfile
  CALL stringlen(FileTemp,FileNameLength)
  IF (GasVelocityFileFormat == 'ContinuousRead') THEN
    READ(23,*,END=1020) (((qzgas(jx,jy,jz),jx=1,nx),jy=1,ny),jz=0,nz)
  ELSE IF (GasVelocityFileFormat == 'SingleColumn') THEN
    DO jz = 0,nz
      DO jy = 1,ny
        DO jx= 1,nx
          READ(23,*,END=1020) qzgas(jx,jy,jz)
        END DO
      END DO
    END DO
  ELSE IF (GasVelocityFileFormat == 'FullForm') THEN
    IF (ny > 1 .AND. nz > 1) THEN
      DO jz = 0,nz
        DO jy = 1,ny
          DO jx= 1,nx
            READ(23,*,END=1020) xdum,ydum,zdum,qzgas(jx,jy,jz)
          END DO
        END DO
      END DO
    ELSE
      CONTINUE
    END IF
  ELSE IF (GasVelocityFileFormat == 'Unformatted') THEN
    READ(23,END=1020) qzgas
  ELSE
    WRITE(*,*)
    WRITE(*,*) ' Z gas velocity file format not recognized'
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  CLOSE(UNIT=23,STATUS='keep')

END IF   !!  End of IF block for Z velocity

GO TO 7004

7001   WRITE(*,*) ' Error opening X gas velocity file: STOP'
READ(*,*)
STOP
7002   WRITE(*,*) ' Error opening Y gas velocity file: STOP'
READ(*,*)
STOP
7003   WRITE(*,*) ' Error opening Z gas velocity file: STOP'
READ(*,*)
STOP

7004   CONTINUE

qxgas = qxgas/(time_scale*dist_scale)
qygas = qygas/(time_scale*dist_scale)
qzgas = qzgas/(time_scale*dist_scale)

END IF   ! End of block for file read of velocities

time_scale = 1.0d0
dist_scale = 1.0d0

!  *****************EROSION/BURIAL BLOCK***********************

ierode = 0
erodex = 0.0d0
erodey = 0.0d0

IF (ALLOCATED(FluidBuryX)) THEN
DEALLOCATE(FluidBuryX)
ALLOCATE(FluidBuryX(0:nx))
ELSE
ALLOCATE(FluidBuryX(0:nx))
END IF

IF (ALLOCATED(FluidBuryY)) THEN
DEALLOCATE(FluidBuryY)
ALLOCATE(FluidBuryY(0:ny))
ELSE
ALLOCATE(FluidBuryY(0:ny))
END IF

IF (ALLOCATED(SolidBuryX)) THEN
DEALLOCATE(SolidBuryX)
ALLOCATE(SolidBuryX(0:nx))
ELSE
ALLOCATE(SolidBuryX(0:nx))
END IF

IF (ALLOCATED(SolidBuryY)) THEN
DEALLOCATE(SolidBuryY)
ALLOCATE(SolidBuryY(0:ny))
ELSE
ALLOCATE(SolidBuryY(0:ny))
END IF

FluidBuryX = 0.0
FluidBuryY = 0.0
SolidBuryX = 0.0
SolidBuryY = 0.0

section = 'erosion/burial'
CALL readblock(nin,nout,section,found,ncount)

IF (found) THEN
WRITE(*,*)
WRITE(*,*) ' Erosion/burial block found'
WRITE(*,*)


CALL units_time(nout,section,time_scale)
CALL units_distance(nout,section,dist_scale)

!  Look for file containing fluid burial rates

BurialFile = ' '
readburial = .false.
CALL read_burialfile(nout,nx,ny,nz,BurialFile,lfile,ierode,BurialFileFormat)

IF (BurialFile /= ' ') THEN
  readburial = .true.
END IF

IF (.NOT. readburial) THEN

  CALL read_erosion(nout,erodex,erodey)

  erodex = erodex/(time_scale*dist_scale)
  erodey = erodey/(time_scale*dist_scale)


  SolidBuryX = erodex
  SolidBuryY = erodey

  NoFluidBury = .TRUE.

  IF (NoFluidBury) THEN
    FluidBuryX = 0.0D0
    FluidBuryY = 0.0D0
  ELSE
    FluidBuryX = erodex
    FluidBuryY = erodey
  END IF

  IF (erodex == 0.0 .AND. erodey == 0.0) THEN
    ierode = 0
  ELSE
    ierode = 1
  END IF

ELSE

  ierode = 1

  WRITE(*,*)
  WRITE(*,*) '  Reading fluid and solid burial rates from file: ',BurialFile(1:lfile)
  WRITE(*,*)
  WRITE(iunit2,*)
  WRITE(iunit2,*) '  Reading fluid and solid burial rates from file: ',BurialFile(1:lfile)
  WRITE(iunit2,*)

  INQUIRE(FILE=BurialFile,EXIST=ext)
  IF (.NOT. ext) THEN
    CALL stringlen(BurialFile,ls)
    WRITE(*,*)
    WRITE(*,*) ' Solid/fluid burial file not found: ', BurialFile(1:ls)
    WRITE(*,*)
    READ(*,*)
    STOP
  END IF
  OPEN(UNIT=52,FILE=BurialFile,STATUS='old',ERR=5015)
  FileTemp = BurialFile
  CALL stringlen(FileTemp,FileNameLength)
!!   ********************************************************
IF (BurialFileFormat == 'ContinuousRead') THEN
  WRITE(*,*)
  WRITE(*,*) ' Only column format allowed for burial/erosion rate format'
  WRITE(*,*)
  READ(*,*)
  STOP
ELSE IF (BurialFileFormat == 'SingleColumn') THEN
  DO jx= 0,nx
    READ(52,*,END=1020) FluidBuryX(jx),SolidBuryX(jx)
  END DO
ELSE IF (BurialFileFormat == 'FullForm') THEN
  DO jx= 0,nx
    READ(52,*,END=1020) xdum,FluidBuryX(jx),SolidBuryX(jx)
  END DO
ELSE IF (BurialFileFormat == 'Unformatted') THEN
  READ(52,END=1020) FluidBuryX
  READ(52,END=1020) SolidBuryX
ELSE
  WRITE(*,*)
  WRITE(*,*) ' Burial file format not recognized'
  WRITE(*,*)
  READ(*,*)
  STOP
END IF
!!   ********************************************************

  CLOSE(UNIT=52,STATUS='keep')

!!  For now, set burial/erosion rates in Y direction to 0

  FluidBuryY = 0.0
  SolidBuryY = 0.0

  FluidBuryX = FluidBuryX/(time_scale*dist_scale)
  FluidBuryY = FluidBuryY/(time_scale*dist_scale)
  SolidBuryX = SolidBuryX/(time_scale*dist_scale)
  SolidBuryY = SolidBuryY/(time_scale*dist_scale)

  GOTO 5016

  5015   WRITE(*,*) ' Error opening burial/erosion file: STOP'
  READ(*,*)
  STOP

  5016 CONTINUE

ENDIF

dist_scale = 1.0d0
time_scale = 1.0d0

END IF

!!!  For erosion and burial, transport mineral/solid properties

!!!IF (gimrt .AND. ierode == 1) THEN

IF (ALLOCATED(specificByGrid)) THEN
  DEALLOCATE(specificByGrid)
END IF
ALLOCATE(specificByGrid(nrct,0:nx+1,0:ny+1,0:nz+1))

IF (ALLOCATED(areainByGrid)) THEN
  DEALLOCATE(areainByGrid)
END IF
ALLOCATE(areainByGrid(nrct,0:nx+1,0:ny+1,0:nz+1))

IF (ALLOCATED(volinByGrid)) THEN
  DEALLOCATE(volinByGrid)
END IF
ALLOCATE(volinByGrid(nrct,0:nx+1,0:ny+1,0:nz+1))

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      !!    Map mineral parameters by Condition to grid
      DO k = 1,nrct
        jpoint = jinit(jx,jy,jz)
        specificByGrid(k,jx,jy,jz) = specific(k,jpoint)
        areainByGrid(k,jx,jy,jz) = areain(k,jpoint)
        volinByGrid(k,jx,jy,jz) = volin(k,jpoint)
      END DO

    END DO
  END DO
END DO

!*****************
!Stolze Lucien: overwrite volin and areain if minerals are imported from dat file
if (readmineral) then
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      DO k = 1,nrct
        areainByGrid(k,jx,jy,jz) = area(k,jx,jy,jz)
        volinByGrid(k,jx,jy,jz) = volfx(k,jx,jy,jz)
      END DO
    END DO
  END DO
END DO
!*****************

ENDIF

!!  Now the boundaries, if NX > 1
IF (nx > 1 .AND. jinit(0,1,1) /= 0 .AND. jinit(nx+1,1,1) /=0) THEN

  DO k = 1,nrct

    jPoint = jinit(0,1,1)
    specificByGrid(k,0,1,1) = specific(k,jPoint)
    areainByGrid(k,0,1,1) = areain(k,jpoint)
    volinByGrid(k,0,jy,jz) = volin(k,jpoint)
    area(k,0,1,1) = areain(k,jpoint)

    jPoint = jinit(nx+1,1,1)
    specificByGrid(k,nx+1,1,1) = specific(k,jPoint)
    areainByGrid(k,nx+1,1,1) = areain(k,jpoint)
    volinByGrid(k,nx+1,jy,jz) = volin(k,jpoint)
    area(k,nx+1,1,1) = areain(k,jpoint)

  END DO

END IF


!!!END IF

!  *****************TRANSPORT BLOCK***********************

section = 'transport'
CALL readblock(nin,nout,section,found,ncount)

!*****************
!Stolze Lucien, June 23, load parameters specific to eastriver simulations
east_river = .FALSE.                  
parchar = 'east_river'
parfind = ' '
CALL read_eastriver(nin,nx,ny,nz)
!***********************
!Stolze Lucien, June 23, choose if solutes are extracted via transpiration (not used for now)
transpisoluteflux = .FALSE.
parchar = 'transpisoluteflux'
parfind = ' '
CALL read_logical(nout,lchar,parchar,parfind,transpisoluteflux)
!*****************

formation = 1.0d0
uli = 1.0d0
dcoeff = 0
dgas = 0
anisotropyY = 1.0d0
!!!anisotropyZ = 1.0d0

UseThresholdPorosity = .FALSE.
  MillingtonQuirk = .TRUE.
TortuosityOption = 'none'

IF (ALLOCATED(tortuosity)) THEN
DEALLOCATE(tortuosity)
ALLOCATE(tortuosity(nx,ny,nz))
ELSE
ALLOCATE(tortuosity(nx,ny,nz))
END IF

tortuosity = 1.0d0

IF (ALLOCATED(anisotropyZ)) THEN
DEALLOCATE(anisotropyZ)
ALLOCATE(anisotropyZ(nx,ny,nz))
ELSE
ALLOCATE(anisotropyZ(nx,ny,nz))
END IF

anisotropyZ = 1.0d0

IF (MontTerri) THEN

  DO jz = 1,nz
    DO jy = 1,ny
      DO jx= 1,nx
             
        IF (jinit(jx,jy,jz) == 7) THEN    
          anisotropyZ(jx,jy,jz) = 3.0   
          tortuosity(jx,jy,jz) =  0.07  
        ELSE IF (jinit(jx,jy,jz) == 6) THEN    
          anisotropyZ(jx,jy,jz) = 3.0   
          tortuosity(jx,jy,jz) =  0.14   
        ELSE IF (jinit(jx,jy,jz) == 5 ) THEN
          anisotropyZ(jx,jy,jz) = 3.0
          tortuosity(jx,jy,jz) =  0.14
        ELSE IF (jinit(jx,jy,jz) == 4 ) THEN
          anisotropyZ(jx,jy,jz) = 1.0
          tortuosity(jx,jy,jz) =  0.05
        ELSE IF (jinit(jx,jy,jz) == 3 ) THEN
          anisotropyZ(jx,jy,jz) = 1.0
          tortuosity(jx,jy,jz) =  0.05
        ELSE IF (jinit(jx,jy,jz) == 2 ) THEN
          anisotropyZ(jx,jy,jz) = 1.0
          tortuosity(jx,jy,jz) =  1.0E-10
        ELSE
          anisotropyZ(jx,jy,jz) = 1.0   !!! CT:
          tortuosity(jx,jy,jz) =  10.0  !!! CT:
        END IF

      END DO
    END DO
  END DO
  
  jz = 0
  DO jy = 1,ny
    DO jx= 1,nx             
      tortuosity(jx,jy,jz) =  0.0
    END DO
  END DO
  
  jz = nz+1
  DO jy = 1,ny
    DO jx= 1,nx             
      tortuosity(jx,jy,jz) =  0.0
    END DO
  END DO
  
  jy = 0
  DO jz = 1,nz
    DO jx= 1,nx             
      tortuosity(jx,jy,jz) =  0.0
    END DO
  END DO
  
  jy = ny+1
  DO jz = 1,nz
    DO jx= 1,nx             
      tortuosity(jx,jy,jz) =  0.0
    END DO
  END DO
  
  jx = 0
  DO jz = 1,nz
    DO jy= 1,ny             
      tortuosity(jx,jy,jz) =  0.0
    END DO
  END DO
  
  jx = nx+1
  DO jz = 1,nz
    DO jy= 1,ny             
      tortuosity(jx,jy,jz) =  0.0
    END DO
  END DO

  
  

  
!!! 1: Chamber
!!! 2: Not in Chamber
!!! 3: Cement
!!! 4: Skin
!!! 5: Outer disturbed zone
!!! 6: Inner disturbed zone
!!! 7: Opalinus Clay (OPA)    
  
END IF

IF (found) THEN

WRITE(*,*)
WRITE(*,*) ' Transport block found'
WRITE(*,*)

IF (nx == 1 .OR. nx == 2) THEN
  WRITE(*,*)
  WRITE(*,*) ' Need at least 3 grid cells for transport'
  WRITE(*,*)
  STOP
END IF

CALL units_time(nout,section,time_scale)
CALL units_distance(nout,section,dist_scale)

CALL read_diffusion(nout,nx,ny,nz)

DO ik = 1,ncomp+nspec
  IF (idiffus == 0) THEN
    d_sp(ik) = dzero
  ELSE
    d_sp(ik) = dcoeff
  END IF
END DO

CALL read_gasdiffusion(nout,nx,ny,nz)

CALL read_speciesdiffusion(nout,ncomp,nspec,ndiff)

IF (ndiff >= 1) THEN
  species_diffusion = .TRUE.
ELSE
  species_diffusion = .FALSE.
END IF

IF (species_diffusion .AND. nsurf > 0) THEN
  WRITE(*,*)
  WRITE(*,*) ' --Combined use of a surface complexation model (not electrically balanced) and '
  WRITE(*,*) ' --electrochemical migration (species-specific diffusion) not recommended.'
  WRITE(*,*) ' --Aqueous phase needs to be electrically neutral'
  WRITE(*,*) ' --ABORTING RUN'
  WRITE(*,*)
!!!    READ(*,*)
!!!    STOP
END IF

CALL read_dispersion(nout,nx,ny,nz,alfl,alft)

alfl = alfl/dist_scale
alft = alft/dist_scale

!  Convert units if necessary (converting to years and meters)

dcoeff = dcoeff/(time_scale*dist_scale*dist_scale)
dzero = dzero/(time_scale*dist_scale*dist_scale)
dgas = dgas/(time_scale*dist_scale*dist_scale)

DO ik = 1,ncomp+nspec
  d_sp(ik) = d_sp(ik)/(time_scale*dist_scale*dist_scale)
END DO

time_scale = 1.0d0
dist_scale = 1.0d0

!!  Section on tortuosity

constant_tortuosity = .FALSE.                  !!  Treat this as false for default, since setting of this to true will disable other "tortuosity" options
call read_ConstantTortuosity(nout,nx,ny,nz,constant_tortuosity,TortuosityConstant,TortuosityOption)

IF (constant_tortuosity) THEN
  WRITE(*,*)
  WRITE(*,*) ' Constant tortuosity option specified'
  WRITE(*,*)
  MillingtonQuirk = .TRUE.
  IF (TortuosityOption /= 'none') THEN
    CALL stringlen(TortuosityOption,ls)
    WRITE(*,*)
    WRITE(*,*) ' Tortuosity will be calculated using: ', TortuosityOption(1:ls)
    WRITE(*,*)
  END IF
  tortuosity = TortuosityConstant
ELSE

!   No constant tortuosity specified, so look for file read or for tortuosity set by zones

  TortuosityFile = ' '
  ReadTortuosity = .FALSE.
  CALL read_TortuosityFile(nout,nx,ny,nz,constant_tortuosity,TortuosityFile,lfile,TortuosityFileFormat)

  IF (TortuosityFile == ' ') THEN
    ReadTortuosity = .FALSE.
  ELSE
    ReadTortuosity = .TRUE.
  END IF

!   Reading tortuosity zones directly from input file
  IF (.NOT. ReadTortuosity) THEN

    ALLOCATE(TortuosityZone(0:mperm))

    TortuosityZone = 0.0d0

    ALLOCATE(jxxTortuosity_lo(mperm))
    ALLOCATE(jxxTortuosity_hi(mperm))
    ALLOCATE(jyyTortuosity_lo(mperm))
    ALLOCATE(jyyTortuosity_hi(mperm))
    ALLOCATE(jzzTortuosity_lo(mperm))
    ALLOCATE(jzzTortuosity_hi(mperm))

    CALL read_TortuosityByZone(nout,nx,ny,nz)

     IF (TortuosityZone(0) == 0.0d0 .AND. nTortuosityZone==0) THEN

!!        WRITE(*,*)
!!        WRITE(*,*) ' No default tortuosity given'
!!        WRITE(*,*) ' Tortuosity should be followed by "default" or blank string'
!!        WRITE(*,*)
!!        STOP

    ELSE
  MillingtonQuirk = .TRUE.
      WRITE(*,*)
      WRITE(*,*) ' Default tortuosity = ',TortuosityZone(0)
      WRITE(*,*)
    END IF

! First, initialize the tortuosity to default tortuosity (TortuosityZone(0))

    IF (TortuosityZone(0) > 0.0d0 .OR. nTortuosityZone > 0) THEN
      MillingtonQuirk = .TRUE.
      Tortuosity = TortuosityZone(0)

!       Next, initialize tortuosity from various zones

      IF (nTortuosityZone > 0) THEN
        DO l = 1,nTortuosityZone
          DO jz = jzzTortuosity_lo(l),jzzTortuosity_hi(l)
            DO jy = jyyTortuosity_lo(l),jyyTortuosity_hi(l)
              DO jx = jxxTortuosity_lo(l),jxxTortuosity_hi(l)
                Tortuosity(jx,jy,jz) = TortuosityZone(l)
              END DO
            END DO
          END DO
        END DO
      END IF

!!      Check to see if any of the nodes are uninitialized with a non-zero value

      CheckSum = MINVAL(Tortuosity)

      IF (checkSum < eps) THEN
        WRITE(*,*)
        WRITE(*,*) ' Tortuosity is not initialized to a non-zero value everywhere'
        WRITE(*,*)
        STOP
      END IF

    ELSE
  MillingtonQuirk = .TRUE.
    END IF

    DEALLOCATE(TortuosityZone)
    DEALLOCATE(jxxTortuosity_lo)
    DEALLOCATE(jxxTortuosity_hi)
    DEALLOCATE(jyyTortuosity_lo)
    DEALLOCATE(jyyTortuosity_hi)
    DEALLOCATE(jzzTortuosity_lo)
    DEALLOCATE(jzzTortuosity_hi)

  ELSE                             !!  Otherwise, read from tortuosity file

    IF (TortuosityFile /= ' ') THEN
      INQUIRE(FILE=TortuosityFile,EXIST=ext)
      IF (.NOT. ext) THEN
        CALL stringlen(TortuosityFile,ls)
        WRITE(*,*)
        WRITE(*,*) ' Tortuosity file not found: ', TortuosityFile(1:ls)
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF
  MillingtonQuirk = .TRUE.
      OPEN(UNIT=52,FILE=TortuosityFile,STATUS='OLD',ERR=6002)
      FileTemp = TortuosityFile
      CALL stringlen(FileTemp,FileNameLength)
      IF (TortuosityFileFormat == 'ContinuousRead') THEN
        READ(52,*,END=1020) (((tortuosity(jx,jy,jz),jx=1,nx),jy=1,ny),jz=1,nz)
        OPEN(UNIT=53,FILE='SynchrotronStructure.dat',STATUS='UNKNOWN')

        jz = 1
          DO jy = 1,ny
            DO jx= 1,nx
              IF (tortuosity(jx,jy,jz) == 1.0d0) THEN   !! pore space
                WRITE(53,*) 'porespace  ',jx,jy,jz, '  fix'
              ELSE                                      !! matrix
                WRITE(53,*) 'matrix     ',jx,jy,jz, '  fix'
              END IF
            END DO
          END DO

        DO jz = 2,nz
          DO jy = 1,ny
            DO jx= 1,nx
              IF (tortuosity(jx,jy,jz) == 1.0d0) THEN   !! pore space
                WRITE(53,*) 'porespace  ',jx,jy,jz
              ELSE                                      !! matrix
                WRITE(53,*) 'matrix     ',jx,jy,jz
              END IF
            END DO
          END DO
        END DO
        CLOSE(UNIT=53,STATUS='keep')

      ELSE IF (TortuosityFileFormat == 'SingleColumn') THEN
        DO jz = 1,nz
          DO jy = 1,ny
            DO jx= 1,nx
              READ(52,*,END=1020) tortuosity(jx,jy,jz)
            END DO
          END DO
        END DO
        OPEN(UNIT=53,FILE='SynchrotronStructure.dat',STATUS='UNKNOWN')
        jz = 1
          DO jy = 1,ny
            DO jx= 1,nx
              IF (tortuosity(jx,jy,jz) == 1.0d0) THEN   !! pore space
                WRITE(53,*) 'porespace  ',jx,jy,jz, '  fix'
              ELSE                                      !! matrix
                WRITE(53,*) 'matrix     ',jx,jy,jz, '  fix'
              END IF
            END DO
          END DO

        DO jz = 2,nz
          DO jy = 1,ny
            DO jx= 1,nx
              IF (tortuosity(jx,jy,jz) == 1.0d0) THEN   !! pore space
                WRITE(53,*) 'porespace  ',jx,jy,jz
              ELSE                                      !! matrix
                WRITE(53,*) 'matrix     ',jx,jy,jz
              END IF
            END DO
          END DO
        END DO
        CLOSE(UNIT=53,STATUS='keep')

      ELSE IF (TortuosityFileFormat == 'FullForm') THEN
        IF (ny > 1 .AND. nz > 1) THEN
          DO jz = 1,nz
            DO jy = 1,ny
              DO jx= 1,nx
                READ(52,*,END=1020) xdum,ydum,zdum,tortuosity(jx,jy,jz)
              END DO
            END DO
          END DO
        OPEN(UNIT=53,FILE='SynchrotronStructure.dat',STATUS='UNKNOWN')
         jz = 1
          DO jy = 1,ny
            DO jx= 1,nx
              IF (tortuosity(jx,jy,jz) == 1.0d0) THEN   !! pore space
                WRITE(53,*) 'porespace  ',jx,jy,jz, '  fix'
              ELSE                                      !! matrix
                WRITE(53,*) 'matrix     ',jx,jy,jz, '  fix'
              END IF
            END DO
          END DO

        DO jz = 2,nz
          DO jy = 1,ny
            DO jx= 1,nx
              IF (tortuosity(jx,jy,jz) == 1.0d0) THEN   !! pore space
                WRITE(53,*) 'porespace  ',jx,jy,jz
              ELSE                                      !! matrix
                WRITE(53,*) 'matrix     ',jx,jy,jz
              END IF
            END DO
          END DO
        END DO
        CLOSE(UNIT=53,STATUS='keep')

        ELSE IF (ny > 1 .AND. nz == 1) THEN
          jz = 1
          DO jy = 1,ny
            DO jx= 1,nx
              READ(52,*,END=1020) xdum,ydum,tortuosity(jx,jy,jz)
            END DO
          END DO
        ELSE
          jz = 1
          jy = 1
          DO jx= 1,nx
            READ(52,*,END=1020) xdum,tortuosity(jx,jy,jz)
          END DO
        END IF
      ELSE IF (TortuosityFileFormat == 'Unformatted') THEN
        READ(52,END=1020) tortuosity
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Tortuosity file format not recognized'
        WRITE(*,*)
        READ(*,*)
        STOP
      END IF

      CLOSE(UNIT=52,STATUS='KEEP')
    END IF
  END IF
END IF

parchar = 'anisotropy_ratioY'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter "anistropy_ratioY" not found
  anisotropyY = 1.000            ! Use default
ELSE
  anisotropyY = realjunk
END IF

parchar = 'anisotropy_ratioZ'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter "anistropy_ratioZ" not found
  anisotropyZ = 1.000            ! Use default
ELSE
  anisotropyZ = realjunk
END IF

ThresholdPorosity = 0.0d0
TortuosityBelowThreshold = 1.0d0
TortuosityAboveThreshold = 1.0d0

parchar = 'threshold_porosity'
parfind = ' '
realjunk = 0.0
CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
IF (parfind == ' ') THEN  ! Parameter "threshold_porosity" not found
  ThresholdPorosity = 0.0d0            ! Use default
ELSE
  ThresholdPorosity = realjunk
END IF

!!  ThresholdPorosity :: real
!!  UseThresholdPorosity :: logical
!!  TortuosityBelowThreshold:  real
!!  TortuosityAboveThreshold:  real

IF (ThresholdPorosity == 0.0d0) THEN
  UseThresholdPorosity = .FALSE.
ELSE
  UseThresholdPorosity = .TRUE.
  MillingtonQuirk = .TRUE.
  parchar = 'tortuosity_below'
  parfind = ' '
  realjunk = 0.0
  CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
  IF (parfind == ' ') THEN  ! Parameter "tortuosity_below" not found
!!      TortuosityBelowThreshold = 1.0d0            ! Use default
    WRITE(*,*)
    WRITE(*,*) ' Non-zero threshold porosity has been specified'
    WRITE(*,*) ' Tortuosity below the threshold, "tortuosity_below", must be given'
    WRITE(*,*)
    STOP
  ELSE
    TortuosityBelowThreshold = realjunk
  END IF
  parchar = 'tortuosity_above'
  parfind = ' '
  realjunk = 0.0
  CALL read_par(nout,lchar,parchar,parfind,realjunk,section)
  IF (parfind == ' ') THEN  ! Parameter "tortuosity_above" not found
!!      TortuosityAboveThreshold = 1.0d0            ! Use default
    WRITE(*,*)
    WRITE(*,*) ' Non-zero threshold porosity has been specified'
    WRITE(*,*) ' Tortuosity above the threshold, "tortuosity_above", must be given'
    WRITE(*,*)
    STOP
  ELSE
    TortuosityAboveThreshold = realjunk
  END IF

END IF

!!  Check for specification of how to calculate mean diffusivity between grid cells

MeanDiffusion = 0
parchar = 'meandiffusion'
parfind = ' '
CALL read_string(nout,lchar,parchar,parfind,dumstring,section)
IF (parfind == ' ') THEN
  MeanDiffusion = 2
ELSE
  IF (dumstring == 'GeometricMean' .OR. dumstring == 'geometricmean' .OR. dumstring == 'geometric') THEN
    MeanDiffusion = 0
  ELSE IF (dumstring == 'ArithmeticMean' .OR. dumstring == 'arithmeticmean' .OR. dumstring == 'arithmetic') THEN
    MeanDiffusion = 1
  ELSE IF (dumstring == 'HarmonicMean' .OR. dumstring == 'harmonicmean' .OR. dumstring == 'harmonic') THEN
    MeanDiffusion = 2
  ELSE
    MeanDiffusion = 2
  END IF
END IF

ELSE
WRITE(*,*)
WRITE(*,*) ' No transport block found'
WRITE(*,*) ' Assuming dispersivity and diffusion = 0'
WRITE(*,*)
alfl = 0.0d0
alft = 0.0d0
dcoeff = 0.0d0
dgas = 0.0d0
dzero = 0.0d0
idiffus = 1
anisotropyY = 1.0d0
anisotropyZ = 1.0d0
END IF

IF (idiffus == 0) THEN
WRITE(*,*)
WRITE(*,*) ' Calculating a temperature-dependent diffusion coefficient'
WRITE(*,*)
ELSE
WRITE(*,*)
WRITE(*,*) '  Using constant diffusion coefficient'
WRITE(*,1110) dcoeff
WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,*) '  Using constant diffusion coefficient'
WRITE(iunit2,1110) dcoeff
WRITE(iunit2,*)
END IF

WRITE(*,*)
WRITE(*,1111) alfl
WRITE(*,1112) alft
WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,1111) alfl
WRITE(iunit2,1112) alft
WRITE(iunit2,*)

!!!CALL dispersivity(nx,ny,nz)
dspx = 0.0
dspy = 0.0
dspz = 0.0

IF (ReadGeochemicalConditions) THEN

!!     do jy = 1,ny
!!       do jx = 1,nx

!!         IF (qx(jx,jy,1) == 0.0d0 .AND. qy(jx,jy-1,1) == 0.0d0 .and. qx(jx-1,jy,1) /= 0.0d0 .and. qy(jx,jy,1) /= 0.0d0) THEN  ! NW face
!!               jinit(jx,jy,1) = 1
!!               volfx(1,jx,jy,1) = 0.20
!!         ELSE IF (qx(jx,jy,1) == 0.0d0 .AND. qy(jx,jy,1) == 0.0d0 .and. qx(jx-1,jy,1) /= 0.0d0 .and. qy(jx,jy-1,1) /= 0.0d0) THEN  ! SW face
!!               jinit(jx,jy,1) = 2
!!               volfx(1,jx,jy,1) = 0.20
!!         ELSE IF (qx(jx-1,jy,1) == 0.0d0 .AND. qy(jx,jy-1,1) == 0.0d0 .and. qx(jx,jy,1) /= 0.0d0 .and. qy(jx,jy,1) /= 0.0d0) THEN  ! NE face
!!               jinit(jx,jy,1) = 1
!!               volfx(1,jx,jy,1) = 0.20
!!         ELSE IF (qx(jx-1,jy,1) == 0.0d0 .AND. qy(jx,jy,1) == 0.0d0 .and. qx(jx,jy,1) /= 0.0d0 .and. qy(jx,jy-1,1) /= 0.0d0) THEN  ! SE face
!!               jinit(jx,jy,1) = 2
!!               volfx(1,jx,jy,1) = 0.20

!!         IF (qx(jx,jy,1) == 0.0d0 .and. qx(jx-1,jy,1) /= 0.0d0) THEN  ! W face
!!               jinit(jx,jy,1) = 1
!!               volfx(1,jx,jy,1) = 0.20
!!         ELSE IF (qy(jx,jy,1) == 0.0d0 .and. qy(jx,jy-1,1) /= 0.0d0 .AND. jy /=1 .AND. jy /= ny) THEN  ! S face
!!               jinit(jx,jy,1) = 2
!!               volfx(1,jx,jy,1) = 0.20
!!         ELSE IF (qy(jx,jy-1,1) == 0.0d0 .and. qy(jx,jy,1) /= 0.0d0 .AND. jy /=1 .AND. jy /= ny) THEN  ! N face
!!               jinit(jx,jy,1) = 1
!!               volfx(1,jx,jy,1) = 0.20
!!         ELSE IF (qx(jx-1,jy,1) == 0.0d0 .AND. qx(jx,jy,1) /= 0.0d0) THEN  ! E face
!!               jinit(jx,jy,1) = 2
!!               volfx(1,jx,jy,1) = 0.20
!!         ELSE IF (qx(jx-1,jy,1) == 0.0d0 .AND. qy(jx,jy,1) == 0.0d0 .and. qx(jx,jy,1) == 0.0d0 .and. qy(jx,jy-1,1) == 0.0d0) THEN  ! Crystal interior
!!             volfx(1,jx,jy,1) = 0.0d0
!!             jinit(jx,jy,1) = 4
!!             activecell(jx,jy,1) = 0
!!             tortuosity(jx,jy,1) = 0.0d0
!!         ELSE
!!             continue
!!         END IF

!!           if (qy(jx,jy,1) /= 0.0d0 .OR. qy(jx,jy-1,1) /= 0.0d0 .OR. qx(jx+1,jy,1) /= 0.0d0 .OR. qx(jx-1,jy,1) /= 0.0d0) then
!!             volfx(1,jx,jy,1) = 0.20
!!             if (jy>(ny/2)) then
!!               jinit(jx,jy,1) = 1
!!             activecell(jx,jy,1) = 0
!!             else
!!               jinit(jx,jy,1) = 2
!!             activecell(jx,jy,1) = 0
!!             end if
!!           else
!!             volfx(1,jx,jy,1) = 0.0d0
!!             jinit(jx,jy,1) = 4
!!             activecell(jx,jy,1) = 0
!!             tortuosity(jx,jy,1) = 0.0d0
!!           end if
!!         end if


!!       end do
!!     end do
END IF


!      write(*,*)
!      do jx = 0,nx
!        write(*,1009) (dspx(jx,jy,1),jy=1,ny)
!      end do
!      write(*,*)
!      pause
!      write(*,*)
!      do jx = 1,nx
!        write(*,1009) (dspy(jx,jy,1),jy=0,ny)
!      end do
!      write(*,*)
!      pause

!! Evapotranspiration

!evapofix = .FALSE.
!evapotimeseries = .FALSE.
!transpifix = .FALSE.
!transpitimeseries = .FALSE.
!CALL read_evaporation(nout,nx,ny,nz,evapofile,evaporate,lfile,evapofix,evapotimeseries,tslength)
!IF (evapotimeseries) THEN
!  IF (ALLOCATED(t_evapo)) THEN
!    DEALLOCATE(t_evapo)
!  END IF
!  IF (ALLOCATED(qt_evapo)) THEN
!    DEALLOCATE(qt_evapo)
!  END IF
!  ALLOCATE(t_evapo(tslength))
!  ALLOCATE(qt_evapo(tslength))
!  CALL read_timeseries(nout,nx,ny,nz,t_evapo,qt_evapo,lfile,evapofile,tslength)
!  evaporate = qt_evapo(1)
!  !STOP
!ENDIF
!CALL read_transpiration(nout,nx,ny,nz,transpifile,transpirate,lfile,transpifix,transpitimeseries,transpicells,tslength)
!IF (transpitimeseries) THEN
!  IF (ALLOCATED(t_transpi)) THEN
!    DEALLOCATE(t_transpi)
!  END IF
!  IF (ALLOCATED(qt_transpi)) THEN
!    DEALLOCATE(qt_transpi)
!  END IF
!  ALLOCATE(t_transpi(tslength))
!  ALLOCATE(qt_transpi(tslength))
!  CALL read_timeseries(nout,nx,ny,nz,t_transpi,qt_transpi,lfile,transpifile,tslength)
!  qt_transpi=qt_transpi/transpicells 
!  transpirate = qt_transpi(1)
!  !WRITE(*,*) transpirate
! ! STOP
!ENDIF




!!!   ******************  NMM Coupling  ****************************************************


!!!   ******************  NMM Coupling  ****************************************************!!!

call BreakthroughInitialize(ncomp,nspec,nkin,nrct,ngas,npot,nx,ny,nz,nseries,  &
                      nexchange,nexch_sec,nsurf,nsurf_sec,ikpH,nplotsurface,nplotexchange )

call MineralBreakthroughInitialize(ncomp,nkin,nrct,nx,ny,nz,minseries)

call IsotopeBreakthroughInitialize(ncomp,nspec,nkin,nrct,nseries,nx,ny,nz)

call AqueousFluxInitialize(ncomp,nspec,nx,ny,nz )

call FluxWeightedConcentrationInitialize(ncomp,nspec,nx,ny,nz )

IF (ALLOCATED(CumulativeXflux)) THEN
DEALLOCATE(CumulativeXflux)
END IF
ALLOCATE(CumulativeXflux(nplotAqueousFlux,nAqueousFluxSeriesFile))
CumulativeXflux = 0.0d0

IF (ALLOCATED(InstantaneousXflux)) THEN
DEALLOCATE(InstantaneousXflux)
END IF
ALLOCATE(InstantaneousXflux(nplotAqueousFlux))
InstantaneousXflux = 0.0d0

IF (ALLOCATED(XfluxWeightedConcentration)) THEN
DEALLOCATE(XfluxWeightedConcentration)
END IF
ALLOCATE(XfluxWeightedConcentration(nplotFluxWeightedConcentration))
XfluxWeightedConcentration = 0.0d0

IF (ALLOCATED(YfluxWeightedConcentration)) THEN
DEALLOCATE(YfluxWeightedConcentration)
END IF
ALLOCATE(YfluxWeightedConcentration(nplotFluxWeightedConcentration))
YfluxWeightedConcentration = 0.0d0

IF (ALLOCATED(ZfluxWeightedConcentration)) THEN
DEALLOCATE(ZfluxWeightedConcentration)
END IF
ALLOCATE(ZfluxWeightedConcentration(nplotFluxWeightedConcentration))
ZfluxWeightedConcentration = 0.0d0


3001 FORMAT('VARIABLES = "Time (yrs)"',                   100(', "',A16,'"'))
3006 FORMAT('  Time          ',                                       100(A17) )
3007 FORMAT('  Yrs           ',                                        100(A17) )
3008 FORMAT('  Days          ',                                        100(A17) )
3009 FORMAT('  Hrs           ',                                        100(A17) )
3010 FORMAT('  Min           ',                                        100(A17) )
3011 FORMAT('  Sec           ',                                        100(A17) )
3002 FORMAT('VARIABLES = "Time (day)"',                   100(', "',A16,'"'))
3012 FORMAT('VARIABLES = "Time (day)" ',11(', "',A18,'" ')               ,2x,     &
                                       '"delCaAqueous     ",   &
                                        "delCaMineral     ",   &
                                        "del34S_Sulfate   ",   &
                                        "del34S_Sulfide   ",   &
                                        "del34S_Mineral   "'     )
3003 FORMAT('VARIABLES = "Time (hrs)"',                   100(', "',A16,'"'))
3004 FORMAT('VARIABLES = "Time (min)"',                   100(', "',A16,'"'))
3005 FORMAT('VARIABLES = "Time (sec)"',                   100(', "',A16,'"'))
3701 FORMAT('  Time(yrs)',4x,150(1X,a16))
3702 FORMAT('  Time(day)',4x,150(1X,a16))
3712 FORMAT('  Time(day)',2x,11(3X,a20),3x,'delCaAqueous        ',3x,     &
                                        'delCaMineral        ',3x,     &
                                        'del34S_Sulfate      ',3x,     &
                                        'del34S_Sulfide      ',3x,     &
                                        'del34S_Mineral      '     )
3036 FORMAT('  Time     ', 2x,11(3X,a20),3x,'delCaAqueous        ',3x,     &
                                        'delCaMineral        ',3x,     &
                                        'del34S_Sulfate      ',3x,     &
                                        'del34S_Sulfide      ',3x,     &
                                        'del34S_Mineral      '     )
3703 FORMAT('  Time(hrs)',4x,150(1X,a16))
3704 FORMAT('  Time(min)',4x,150(1X,a16))
3705 FORMAT('  Time(sec)',4x,150(1X,a16))

!DEALLOCATE(eqgas)
!DEALLOCATE(eqhom)
!DEALLOCATE(eqsurf)

#endif
! end of block to skip for ALQUIMIA
IF (ALLOCATED(stringarray)) THEN
DEALLOCATE(stringarray)
END IF
DEALLOCATE(condtitle)
#ifndef ALQUIMIA
DEALLOCATE(condlabel)
DEALLOCATE(keqmin_tmp)
DEALLOCATE(keqaq_tmp)
DEALLOCATE(keqgas_tmp)
DEALLOCATE(keqsurf_tmp)
DEALLOCATE(sptmp)
DEALLOCATE(sptmp10)
DEALLOCATE(gamtmp)
DEALLOCATE(stmp)
#endif
DEALLOCATE(nbasin)
DEALLOCATE(nbkin)
DEALLOCATE(dxxt)
DEALLOCATE(dyyt)
DEALLOCATE(dzzt)
DEALLOCATE(nvx)
DEALLOCATE(nvy)
DEALLOCATE(nvz)
#ifndef ALQUIMIA
DEALLOCATE(sexch)
DEALLOCATE(spextmp)
DEALLOCATE(spextmp10)
DEALLOCATE(totextmp)
DEALLOCATE(spgastmp)
DEALLOCATE(spgastmp10)
DEALLOCATE(sgastmp)
DEALLOCATE(spsurftmp)
DEALLOCATE(spsurftmp10)
DEALLOCATE(ssurftmp)
DEALLOCATE(dpsi)
DEALLOCATE(ctot)
DEALLOCATE(guess)
DEALLOCATE(itype)
!!!DEALLOCATE(c_surf)
DEALLOCATE(guess_surf)
#endif
DEALLOCATE(gaspp)
#ifndef ALQUIMIA
DEALLOCATE(totexch)
!!DEALLOCATE(cec)
DEALLOCATE(ncon)
#endif
DEALLOCATE(namdep_nyf)
#ifndef ALQUIMIA
DEALLOCATE(tempcond)
DEALLOCATE(SkipAdjust)
DEALLOCATE(rocond)
DEALLOCATE(porcond)
DEALLOCATE(SaturationCond)
DEALLOCATE(PressureCond)
DEALLOCATE(equilibrate)
DEALLOCATE(fsurftmp)
#endif
DEALLOCATE(constraint)
DEALLOCATE(realmult)
DEALLOCATE(pH)
DEALLOCATE(guesspH)
DEALLOCATE(ndist)
DEALLOCATE(jxxlo)
DEALLOCATE(jxxhi)
DEALLOCATE(jyylo)
DEALLOCATE(jyyhi)
DEALLOCATE(jzzlo)
DEALLOCATE(jzzhi)
DEALLOCATE(jjfix)
#ifndef ALQUIMIA
DEALLOCATE(surfcharge_init)
DEALLOCATE(LogPotential_tmp)
#endif
DEALLOCATE(unitsflag)
#ifndef ALQUIMIA
DEALLOCATE(conversion)
#endif
DEALLOCATE(OneOverMassFraction)
!!DEALLOCATE(SolidDensity)
#ifndef ALQUIMIA
DEALLOCATE(SolidSolutionRatio)
#endif
DEALLOCATE(SolidDensityFrom)
IF (Duan .OR. Duan2006) THEN
DEALLOCATE(vrInitial)
END IF
CLOSE(UNIT=8)

2001 FORMAT(2X,'Darcian flux (m**3/m**2/yr) = ',1PE11.2)
2011 FORMAT(2X,'Dispersivity (m) =            ',f11.3)
2002 FORMAT(2X,'Temperature (C) at J = 1      ',f11.2)
2012 FORMAT(2X,'Temperature gradient (C/m) =  ',f11.2)
2003 FORMAT(2X,'Diffusion coeff (m**2/sec) =  ',1PE11.3)
2004 FORMAT(2X,'Formation factor =            ',f11.2)
2031 FORMAT(2X,'Porosity dependence in D* =   ',f11.2)

!**********************
1  FORMAT(1(/))
2  FORMAT(2(/))
3  FORMAT(3(/))
4  FORMAT(4(/))
5  FORMAT(5(/))
6  FORMAT(6(/))
7  FORMAT(7(/))
8  FORMAT(8(/))
9  FORMAT(9(/))
10  FORMAT(10(/))
11  FORMAT(11(/))
12  FORMAT(12(/))
13  FORMAT(13(/))

705 FORMAT(1X,1PE12.5,150(1X,1PE12.4))
!706 FORMAT('#   Time (yrs) ',150(1X,a12))
706 FORMAT('  Time      ',3x,150(1X,a12))
707 FORMAT(100(1X,1PE12.4))
1013 FORMAT(1X,'Graphics files printed at (yrs)')
1014 FORMAT(20(1X,1PE12.4,/))
1009 FORMAT(3X,'Parameters for this run')
1010 FORMAT(2X,a72)
1011 FORMAT(2X,'Maximum timestep (yrs) = ',1PE12.3)
1015 FORMAT(2X,'Initial timestep (yrs) = ',1PE12.3)
1012 FORMAT(2X,'Maximum concentration correction = ',1PE12.4)
1016 FORMAT(1X,'Number of zones in X direction ',i2)
1026 FORMAT(1X,'Number of zones in Y direction ',i2)
1036 FORMAT(1X,'Number of zones in Z direction ',i2)
1018 FORMAT(2X,'Zone',3X,'Width (m)')
1017 FORMAT(4X,i2,1X,1PE12.4)
2101 FORMAT(2X, ' ---> Using database: ',a100)
1023 FORMAT(1X,a18,1X,1PE12.4)
1022 FORMAT(1X,a18,' on ',a18)
598  FORMAT(2X,a18,2X,f10.5,7X,1PE10.3)
599  FORMAT(2X,'Porosity for this chemical condition = ',f8.4)
597  FORMAT(2X,' Mineral          ',2X,'Volume fraction',2X, 'Area (m2/m3)')
595  FORMAT(2X,a18,2X,i2,2X,1PE12.4,2X,1PE12.4,5X,a18)
596  FORMAT(2X,'Primary species',4X,'itype',2X,'   guess ',2X,  &
  '   Total conc.',2X,'Constraint phase')
802 FORMAT(2X,1PE12.4,2X,1PE12.3)
1101 FORMAT(2X,a18)
1102 FORMAT(2X,'Number of parallel reactions = ',i2)
1103 FORMAT(2X,'Sat. state dep.: ',1X,'Sat1 = ',1PE9.2,1X, 'Sat2 = ',1PE9.2)
1104 FORMAT(2X,'Nucleation threshold (log Q/K) = ',1PE10.2)
1105 FORMAT(4X,'Reaction ',i2,' includes ',i2,  &
  ' far from equilibrium species')
1106 FORMAT(6X,'Rate constant (mol/m**2/s) ',1PE12.4)
1108 FORMAT(6X,'Activation energy (kcal/mol) ',1PE12.4)
1107 FORMAT(6X,10(a18,1X,1PE10.2))
5010 FORMAT(1X,' Fixing porosity at ',f7.4)
5012 FORMAT(' Using temperature of ',f6.2)
5013 FORMAT(' Using temperature gradient of ',f6.2)
1110 FORMAT(2X,'Diffusion coefficient (m**2/yr) =     ',1PE12.4)
1111 FORMAT(2X,'Longitudinal dispersivity (m) =   ',1PE12.4)
1112 FORMAT(2X,'Transverse dispersivity (m) =     ',1PE12.4)
1113 FORMAT(2X,'Porosity dependence of diff. coeff.= ',1PE12.4)

RETURN

700 WRITE(*,*)
WRITE(*,*) ' Error opening temperature file'
WRITE(*,*)
READ(*,*)
STOP
701 WRITE(*,*)
WRITE(*,*) ' Error reading temperature file'
WRITE(*,*)
READ(*,*)
STOP

703 WRITE(*,*)
WRITE(*,*) ' Error opening input file'
WRITE(*,*) ' Looking for: ', filename
WRITE(*,*)
READ(*,*)
STOP
704 WRITE(*,*)
WRITE(*,*) ' Error opening ',FileOutput
WRITE(*,*)
READ(*,*)
STOP
708 WRITE(*,*)
WRITE(*,*) ' Error opening PEST control file (PestControl.ant)'
WRITE(*,*)
READ(*,*)
STOP
3000 WRITE(*,*)
WRITE(*,*) ' Error reading MODFLOW stress period information'
WRITE(*,*) ' Number of stress periods: ', nstress
WRITE(*,*) ' Looking for information on period length, no. of time steps, and multiplier '
WRITE(*,*)
READ(*,*)
STOP
6001 WRITE(*,*)
CALL stringlen(PorosityFile,ls)
WRITE(*,*) ' Error opening porosity file: ',PorosityFile(1:ls)
WRITE(*,*)
READ(*,*)
STOP
6002 WRITE(*,*)
CALL stringlen(SaturationFile,ls)
WRITE(*,*) ' Error opening liquid saturation file: ',SaturationFile(1:ls)
WRITE(*,*)
READ(*,*)
STOP

1020  WRITE(*,*) ' End of file during read'
WRITE(*,*) ' Trying to read the file: ', FileTemp(1:FileNameLength)
READ(*,*)
STOP

END SUBROUTINE StartTope
