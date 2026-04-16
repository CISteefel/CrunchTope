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

SUBROUTINE CrunchTope(NumInputFiles,InputFileCounter,NewInput)
USE crunchtype
USE params
USE runtime
USE crunch_interface
USE concentration
USE mineral
USE solver
USE medium
USE transport
USE flow
USE temperature
USE io
USE ReadFlow
USE modflowModule
USE CrunchFunctions
USE Richards_module
!!USE fparser

#include <petsc/finclude/petsc.h>
      use petsc

IMPLICIT NONE

INTERFACE
  SUBROUTINE bdgas(ncomp,nspec,nrct,ngas,nbnd,scorr)
  USE crunchtype
  IMPLICIT NONE
  INTEGER(I4B), INTENT(IN)                      :: ncomp
  INTEGER(I4B), INTENT(IN)                      :: nspec
  INTEGER(I4B), INTENT(IN)                      :: nrct
  INTEGER(I4B), INTENT(IN)                      :: ngas
  INTEGER(I4B), INTENT(IN)                      :: nbnd
  REAL(DP), DIMENSION(:), INTENT(OUT)           :: scorr
  END SUBROUTINE bdgas
END INTERFACE

LOGICAL(LGT), INTENT(IN OUT)                                      :: NewInput
INTEGER(I4B), INTENT(IN OUT)                                      :: NumInputFiles
INTEGER(I4B), INTENT(IN OUT)                                      :: InputFileCounter

EXTERNAL AssembleGlobal
EXTERNAL SolveDiffuse
EXTERNAL dgetrf
EXTERNAL dgetrs
EXTERNAL CrunchPETScInitializePressure
EXTERNAL CrunchPETScInitializeChemistry
EXTERNAL CrunchPETScInitializeDiffusion
EXTERNAL CrunchPETScFinalizeSolver
EXTERNAL CrunchPETScTolerances

! *******************begin PETSc declarations of f90 variables***********
INTEGER(I4B)             ::numprocs
INTEGER(I4B)             ::ierr
INTEGER(I4B)             ::irank
INTEGER(I4B)             ::linefil
INTEGER(I4B)             ::itsiterate
INTEGER(I4B), PARAMETER  ::maxitsksp=200
REAL(DP), PARAMETER      ::zero=0.0
!!REAL(DP), PARAMETER      ::rtolksp=1.0d-09
!!REAL(DP), PARAMETER      ::atolksp=1.0d-50
!!REAL(DP), PARAMETER      ::dtolksp=1.0d05
REAL(DP)                 :: rtolksp
REAL(DP)                 :: atolksp
REAL(DP)                 :: dtolksp

!*********************end PETSc declarations ******************************

CHARACTER (LEN=12)                                         :: dumm1
CHARACTER (LEN=12)                                         :: dumm2
CHARACTER (LEN=12)                                         :: dumm3
INTEGER(I4B), DIMENSION(8)                                 :: curr_time

INTEGER(I4B)                                               :: str_mon
INTEGER(I4B)                                               :: str_day
INTEGER(I4B)                                               :: str_hr
INTEGER(I4B)                                               :: str_min
INTEGER(I4B)                                               :: str_sec
INTEGER(I4B)                                               :: str_millisec
INTEGER(I4B)                                               :: end_mon
INTEGER(I4B)                                               :: end_day
INTEGER(I4B)                                               :: end_hr
INTEGER(I4B)                                               :: end_min
INTEGER(I4B)                                               :: end_sec
INTEGER(I4B)                                               :: end_millisec
REAL(DP)                                                   :: StartSeconds
REAL(DP)                                                   :: EndSeconds
REAL(DP)                                                   :: MilliSeconds
INTEGER(I4B), DIMENSION(1)                                              :: MaxFx

REAL(DP), DIMENSION(:), ALLOCATABLE                        :: tdata
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: break
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: c0
REAL(DP), DIMENSION(:), ALLOCATABLE                        :: csave
REAL(DP), DIMENSION(:,:), ALLOCATABLE                      :: aaaTemp
REAL(DP), PARAMETER                                        :: atol=1.e-09
REAL(DP), PARAMETER                                        :: rtol=1.e-06
REAL(DP), PARAMETER                                        :: correx=1.0
REAL(DP), PARAMETER                                        :: tiny=1.0E-13

INTEGER(I4B), PARAMETER                                    :: idetail=0
INTEGER(I4B), PARAMETER                                    :: nitmax=100
INTEGER(I4B), PARAMETER                                    :: north=1
INTEGER(I4B), PARAMETER                                    :: nend=999999999
INTEGER(I4B), PARAMETER                                    :: maxiterat=1
INTEGER(I4B), PARAMETER                                    :: newton=50
INTEGER(I4B), PARAMETER                                    :: iprint3=1

LOGICAL(LGT)                                               :: solve_flag
LOGICAL(LGT)                                               :: ext

INTEGER(I4B)                                               :: ncomp
INTEGER(I4B)                                               :: nspec
INTEGER(I4B)                                               :: nkin
INTEGER(I4B)                                               :: nrct
INTEGER(I4B)                                               :: ngas
INTEGER(I4B)                                               :: nx
INTEGER(I4B)                                               :: ny
INTEGER(I4B)                                               :: nz
INTEGER(I4B)                                               :: ipath
INTEGER(I4B)                                               :: igamma
INTEGER(I4B)                                               :: ikmast
INTEGER(I4B)                                               :: ikph
INTEGER(I4B)                                               :: ikO2
INTEGER(I4B)                                               :: jpor
INTEGER(I4B)                                               :: ikin
INTEGER(I4B)                                               :: nstop
INTEGER(I4B)                                               :: nstopsave
INTEGER(I4B)                                               :: nseries
INTEGER(I4B)                                               :: minseries
INTEGER(I4B)                                               :: nexchange
INTEGER(I4B)                                               :: nexch_sec
INTEGER(I4B)                                               :: nsurf
INTEGER(I4B)                                               :: nsurf_sec
INTEGER(I4B)                                               :: npot
INTEGER(I4B)                                               :: ndecay
INTEGER(I4B)                                               :: ijunk1
INTEGER(I4B)                                               :: iures

CHARACTER (LEN=mls)                                        :: data1
CHARACTER (LEN=2*mls)                                      :: ltitle
CHARACTER (LEN=mls)                                        :: dumstring
CHARACTER (LEN=mls)                                        :: file
CHARACTER (LEN=mls)                                        :: prefix

REAL(DP)                                                   :: tstep
REAL(DP)                                                   :: delt
REAL(DP)                                                   :: deltmin
REAL(DP)                                                   :: ttol
REAL(DP)                                                   :: corrmax
REAL(DP)                                                   :: time
REAL(DP)                                                   :: timeflow
REAL(DP)                                                   :: dtflow
REAL(DP)                                                   :: MaxDivergence
REAL(DP)                                                   :: MaximumCorrection
REAL(DP)                                                   :: TotalMass
REAL(DP)                                                   :: InitialTotalMass
REAL(DP)                                                   :: dt_GIMRT
REAL(DP)                                                   :: t_rich
REAL(DP)                                                   :: PrintSeconds
REAL(DP)                                                   :: ExtraSeconds
REAL(DP)                                                   :: StartSimulation

INTEGER(I4B)                                               :: reason
INTEGER(I4B)                                               :: nflow
INTEGER(I4B)                                               :: nt
INTEGER(I4B)                                               :: neqn
INTEGER(I4B)                                               :: iprnt
INTEGER(I4B)                                               :: nint
INTEGER(I4B)                                               :: nxyz
INTEGER(I4B)                                               :: np
INTEGER(I4B)                                               :: iteration_tot
INTEGER(I4B)                                               :: nn
INTEGER(I4B)                                               :: i
INTEGER(I4B)                                               :: jx
INTEGER(I4B)                                               :: jy
INTEGER(I4B)                                               :: jz
INTEGER(I4B)                                               :: j
INTEGER(I4B)                                               :: jydum
INTEGER(I4B)                                               :: newtmax
INTEGER(I4B)                                               :: ne
INTEGER(I4B)                                               :: icvg
!!INTEGER(I4B)                                               :: iterat
INTEGER(I4B)                                               :: i2
INTEGER(I4B)                                               :: ind
INTEGER(I4B)                                               :: ier
INTEGER(I4B)                                               :: ix
INTEGER(I4B)                                               :: nex
INTEGER(I4B)                                               :: is
INTEGER(I4B)                                               :: itskip
INTEGER(I4B)                                               :: ik
INTEGER(I4B)                                               :: ineg
INTEGER(I4B)                                               :: k
INTEGER(I4B)                                               :: jxph
INTEGER(I4B)                                               :: jyph
INTEGER(I4B)                                               :: jzph
INTEGER(I4B)                                               :: ll
INTEGER(I4B)                                               :: intfile
INTEGER(I4B)                                               :: isimu_hr
INTEGER(I4B)                                               :: isimu_min
INTEGER(I4B)                                               :: isimu_sec

INTEGER(I4B)                                               :: ndata
INTEGER(I4B)                                               :: nchar
INTEGER(I4B)                                               :: n3
INTEGER(I4B)                                               :: ntt
INTEGER(I4B)                                               :: npt
INTEGER(I4B)                                               :: ls
INTEGER(I4B)                                               :: ncount
INTEGER(I4B)                                               :: idummy
INTEGER(I4B)                                               :: jxmax
INTEGER(I4B)                                               :: jymax
INTEGER(I4B)                                               :: jydum2
INTEGER(I4B)                                               :: jyCheck

REAL(DP)                                                   :: dtold
REAL(DP)                                                   :: tempc
REAL(DP)                                                   :: det
REAL(DP)                                                   :: errmax
REAL(DP)                                                   :: tolmax
REAL(DP)                                                   :: ddtold
REAL(DP)                                                   :: waterold
REAL(DP)                                                   :: waternew
REAL(DP)                                                   :: ratio
REAL(DP)                                                   :: dtnewest
REAL(DP)                                                   :: rinv
REAL(DP)                                                   :: dxe
REAL(DP)                                                   :: dxw
REAL(DP)                                                   :: phmax
REAL(DP)                                                   :: phchg
REAL(DP)                                                   :: phwrite
REAL(DP)                                                   :: simu_t
REAL(DP)                                                   :: slope
REAL(DP)                                                   :: dtmin
REAL(DP)                                                   :: dtemp
REAL(DP)                                                   :: check
REAL(DP)                                                   :: checkN
REAL(DP)                                                   :: checkS
REAL(DP)                                                   :: checkE
REAL(DP)                                                   :: checkW
REAL(DP)                                                   :: checkPlus
REAL(DP)                                                   :: checkMinus
REAL(DP)                                                   :: RealSum
REAL(DP)                                                   :: aascale
REAL(DP)                                                   :: dtmax
REAL(DP)                                                   :: dtmaxcour
REAL(DP)                                                   :: totcharge
REAL(DP)                                                   :: PrintTime
REAL(DP)                                                   :: DummyReal
REAL(DP)                                                   :: distCheck
REAL(DP)                                                   :: distCheck3
REAL(DP)                                                   :: distCheck2
REAL(DP), PARAMETER                                        :: eps=1.D-12

REAL(DP)                                                   :: minsat
REAL(DP)                                                   :: dtnuft
REAL(DP)                                                   :: porsatro
REAL(DP)                                                   :: dt
REAL(DP)                                                   :: RoAveLeft
REAL(DP)                                                   :: RoAveRight
REAL(DP)                                                   :: qxSum

INTEGER(I4B)                                               :: jyy

CHARACTER (LEN=1)                                           :: Coordinate

INTEGER(I4B)                                               :: phloc(3)
INTEGER(I4B)                                               :: DummyInteger
INTEGER(I4B)                                               :: nplotsurface
INTEGER(I4B)                                               :: nplotexchange

LOGICAL(LGT)                                               :: SteadyFlow
LOGICAL(LGT)                                               :: ActiveFlow
LOGICAL(LGT)                                               :: TrueFalse
LOGICAL(LGT)                                               :: steady
LOGICAL(LGT)                                               :: scale
LOGICAL(LGT)                                               :: FlowConverged
LOGICAL(LGT)                                               :: FirstCall

REAL(DP), DIMENSION(:), ALLOCATABLE                        :: tempreal
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: InitialMass

REAL(DP)                                                   :: satchange
REAL(DP)                                                   :: sumvertical
REAL(DP)                                                   :: sumhorizontal
REAL(DP)                                                   :: SumCheck
REAL(DP)                                                   :: massbalance
REAL(DP)                                                   :: RealDummy
REAL(DP)                                                   :: AqueousToBulk
REAL(DP)                                                   :: adjust
REAL(DP)                                                   :: pi
REAL(DP)                                                   :: SumSlice
REAL(DP)                                                   :: SumQx
REAL(DP)                                                   :: MaxPermeabilityX
REAL(DP)                                                   :: MaxPermeabilityY
REAL(DP)                                                   :: MinPermeabilityX
REAL(DP)                                                   :: MinPermeabilityY
REAL(DP)                                                   :: ConcentrationNormalized
REAL(DP)                                                   :: Nr
REAL(DP)                                                   :: ScaledNr
REAL(DP)                                                   :: CPU_unit
REAL(DP)                                                   :: TimeHours
REAL(DP)                                                   :: DelCaAqueous
REAL(DP)                                                   :: delCaMineralInstant
REAL(DP)                                                   :: delCaMineralBulk
REAL(DP)                                                   :: del34S_sulfate
REAL(DP)                                                   :: del34S_sulfide
REAL(DP)                                                   :: del34S_SulfideMineral
REAL(DP)                                                   :: rateFeS
REAL(DP)                                                   :: rateS
REAL(DP)                                                   :: CumulativeSulfate
REAL(DP)                                                   :: CumulativeFe
REAL(DP)                                                   :: CumulativeO2
REAL(DP)                                                   :: CumulativeCO2

CHARACTER (LEN=1)                                          :: trans

INTEGER(I4B)                                               :: info
INTEGER(I4B)                                               :: npz
INTEGER(I4B)                                               :: PetscLogPrintSummary
INTEGER(I4B)                                               :: jxslice
INTEGER(I4B)                                               :: nbnd

INTEGER(I4B), PARAMETER                                    :: ione=1

REAL(DP)                                                   :: Tk
REAL(DP)                                                   :: MassBalanceError
REAL(DP)                                                   :: AqueousMassBalanceError
REAL(DP)                                                   :: RateMassBalanceError
REAL(DP)                                                   :: DynamicMassBalanceError
REAL(DP)                                                   :: aq_accum
REAL(DP)                                                   :: sumCalcite
REAL(DP)                                                   :: sumCO2
REAL(DP)                                                   :: sumPlagioclaseArea
REAL(DP)                                                   :: denominator
REAL(DP)                                                   :: totPor
REAL(DP)                                                   :: totVol
REAL(DP)                                                   :: totCheck
REAL(DP)                                                   :: totChange
REAL(DP)                                                   :: totCheckInitial
REAL(DP)                                                   :: CheckMass1
REAL(DP)                                                   :: CheckMass2
REAL(DP)                                                   :: pumpterm

INTEGER(I4B)                                               :: nBoundaryConditionZone
INTEGER(I4B)                                               :: nco
INTEGER(I4B)                                               :: i_substep
INTEGER(I4B)                                               :: kk

! transient pump time series (Lucien Stolze)
REAL(DP)                                                    :: time_norm
REAL(DP)                                                    :: qgdum
REAL(DP)                                                    :: time_dum
! transient temperature (Lucien Stolze)
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: temp_dum
! transient water table (Lucien Stolze)
REAL(DP)                                                    :: wattab
REAL(DP), DIMENSION(:), ALLOCATABLE                         :: depth
INTEGER(I4B)                                                :: depthwattab

CHARACTER (LEN=3)                                           :: ulabPrint
REAL(DP)                                                    :: sionPrint

REAL(DP), DIMENSION(:), ALLOCATABLE                          :: GasFlux_FaceWest
REAL(DP), DIMENSION(:), ALLOCATABLE                          :: GasFlux_FaceEast
REAL(DP), DIMENSION(:), ALLOCATABLE                          :: GasFlux_FaceSouth
REAL(DP), DIMENSION(:), ALLOCATABLE                          :: GasFlux_FaceNorth
REAL(DP), DIMENSION(:), ALLOCATABLE                          :: AqueousFlux_FaceWest
REAL(DP), DIMENSION(:), ALLOCATABLE                          :: AqueousFlux_FaceEast
REAL(DP), DIMENSION(:), ALLOCATABLE                          :: AqueousFlux_FaceSouth
REAL(DP), DIMENSION(:), ALLOCATABLE                          :: AqueousFlux_FaceNorth

!*************************************************************************
! Edit by Toshiyuki Bandai, 2024 Oct.
INTEGER(I4B)                                                  :: lowX
INTEGER(I4B)                                                  :: lowY
INTEGER(I4B)                                                  :: lowZ
INTEGER(I4B)                                                  :: highX
INTEGER(I4B)                                                  :: highY
INTEGER(I4B)                                                  :: highZ
CHARACTER (LEN=15)                                            :: text
! End of Edit by Toshiyuki Bandai, 2024 Oct. 
!*************************************************************************

!*************************************************************************
! Added by Toshiyuki Bandai, 2024 Aug. to measure walltime for time stepping
!character(10) :: time_WallTime1
!character(10) :: time_WallTime2
!*************************************************************************


! ******************** PETSC declarations ********************************
PetscFortranAddr     userC(6),userD(6),userP(6),user(6)
Mat                  amatpetsc,amatD,amatP
Vec                  bvec,xvec,bvecD,xvecD,bvecP,xvecP
!!SLES                 sles
PC                   pc
KSP                  ksp
!!Scalar               zeroPetsc
! ************************end PETSc declarations of PETSc variables ******
  
Switcheroo = .FALSE.

MassBalanceError = 0.0d0
AqueousMassBalanceError = 0.0d0
RateMassBalanceError = 0.0d0
CumulativeSulfate = 0.0d0
CumulativeFe = 0.0d0
CumulativeO2 = 0.0d0
CumulativeCO2 = 0.0d0

INQUIRE(FILE='PestControl.ant',EXIST=ext)
IF (EXT) THEN
  RunningPest = .TRUE.
ELSE
  RunningPest = .FALSE.
END IF

WRITE(*,*)
WRITE(*,*) '   ************************** CrunchFlow ******************************'
WRITE(*,*)
WRITE(*,*) '                 Authors:  C.I. STEEFEL, T. BANDAI, S. MOLINS'

WRITE(*,*) '                      *** Copyright Notice ***          '
WRITE(*,*) '     Copyright (c) 2016, The Regents of the University of California, '
WRITE(*,*) '            through Lawrence Berkeley National Laboratory             '
WRITE(*,*) '                        All Rights Reserved                           '
WRITE(*,*)
WRITE(*,*) '                  For full license agreement with disclaimers:        '
WRITE(*,*) '                  See https://github.com/CISteefel/CrunchTope         '
WRITE(*,*) '   ********************************************************************'

IF (.NOT. RunningPest) THEN

!!  Go here when running subsequent input files

  9999 CONTINUE

  call CPU_TIME(PrintSeconds)
  StartSimulation = PrintSeconds
  StartSimulation = 0.0d0

END IF    !  Above is not used when running PEST

FirstCall = .TRUE.

!!  Initialize
trans = 'N'
solve_flag = .FALSE.
steady = .FALSE.
OS3Dpetsc = .FALSE.
pi = DACOS(-1.0D0)
ncomp = 0
ngas = 0
nspec = 0
ndecay = 0
nsurf = 0
nexchange = 0
nexch_sec = 0
nsurf_sec = 0
nrct = 0
nkin = 0
ikin = 0
NumSourceTerms = 0
isaturate = 0
dtmaxcour = 0.0
iprnt = 0
ncounter = 0

CALL StartTope(ncomp,nspec,nkin,nrct,ngas,npot,nx,ny,nz,data1,ipath,igamma,               &
    ikmast,ikph,iko2,ltitle,tstep,delt,deltmin,ttol,jpor,ikin,nstop,                      &
    corrmax,nseries,minseries,nexchange,nexch_sec,nsurf,nsurf_sec,ndecay,str_mon,         &
    str_day,str_hr,str_min,str_sec,NumInputFiles,InputFileCounter,nBoundaryConditionZone)

str_mon = 0
str_day = 0
str_hr  = 0
str_min = 0
str_sec = 0
str_millisec = 0

IF (ALLOCATED(GasFlux_FaceWest)) THEN
  DEALLOCATE(GasFlux_FaceWest)
END IF
ALLOCATE(GasFlux_FaceWest(ngas))
GasFlux_FaceWest = 0.0

IF (ALLOCATED(GasFlux_FaceEast)) THEN
  DEALLOCATE(GasFlux_FaceEast)
END IF
ALLOCATE(GasFlux_FaceEast(ngas))
GasFlux_FaceEast = 0.0

IF (ALLOCATED(GasFlux_FaceSouth)) THEN
  DEALLOCATE(GasFlux_FaceSouth)
END IF
ALLOCATE(GasFlux_FaceSouth(ngas))
GasFlux_FaceSouth = 0.0

IF (ALLOCATED(GasFlux_FaceNorth)) THEN
  DEALLOCATE(GasFlux_FaceNorth)
END IF
ALLOCATE(GasFlux_FaceNorth(ngas))
GasFlux_FaceNorth = 0.0

!!!!!!!!!!!!!!!!!!!!
IF (ALLOCATED(AqueousFlux_FaceWest)) THEN
  DEALLOCATE(AqueousFlux_FaceWest)
END IF
ALLOCATE(AqueousFlux_FaceWest(ncomp+nspec))
AqueousFlux_FaceWest = 0.0

IF (ALLOCATED(AqueousFlux_FaceEast)) THEN
  DEALLOCATE(AqueousFlux_FaceEast)
END IF
ALLOCATE(AqueousFlux_FaceEast(ncomp+nspec))
AqueousFlux_FaceEast = 0.0

IF (ALLOCATED(AqueousFlux_FaceSouth)) THEN
  DEALLOCATE(AqueousFlux_FaceSouth)
END IF
ALLOCATE(AqueousFlux_FaceSouth(ncomp+nspec))
AqueousFlux_FaceSouth = 0.0

IF (ALLOCATED(AqueousFlux_FaceNorth)) THEN
  DEALLOCATE(AqueousFlux_FaceNorth)
END IF
ALLOCATE(AqueousFlux_FaceNorth(ncomp+nspec))
GasFlux_FaceNorth = 0.0

! ************ Initialize PETSc stuff ***************************************
IF ( InputFileCounter == 1) THEN
  
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call MPI_COMM_SIZE(PETSC_COMM_WORLD,numprocs,ierr)
  IF (numprocs /= 1) THEN
    call MPI_COMM_RANK(PETSC_COMM_WORLD,irank,ierr)
    if (irank == 0) THEN
      write(*,*) ' This is the 1-processor case'
    endif
    call PetscFinalize(ierr)
    READ(*,*)
    STOP
    
  END IF
END IF
!  ************** End PETSc initialize ***************************************

!!! Last two unknowns are ionic strength (sion[jx,jy,jz] and activity of water (gammawater[jx,jy,jz] )
neqn = ncomp + nexchange + nsurf + npot + 1 + 1

IF (nstop == 0) THEN
  
  WRITE(*,*)
  WRITE(*,*) '       DISTRIBUTION OF SPECIES COMPLETED'
  WRITE(*,*) '          NO TIME STEPPING SPECIFIED'
  WRITE(*,*)
!*******************petsc addition to cleanup *******************************************
  call PetscFinalize(ierr)
!******************end petsc addition for cleanup on error ******************************
  READ(*,*)
  STOP
  
END IF

dtmax = tstep
nint = 1
OutputTimeCounter = 1
nxyz = nx*ny*nz
np = nxyz*neqn
dtold = delt
time = 0.0

WRITE(*,*)
WRITE(*,*) '        INITIALIZATION COMPLETED '
WRITE(*,*) '        --> STARTING TIMESTEPPING '
WRITE(*,*)
WRITE(iunit2,*)
WRITE(iunit2,*) '        INITIALIZATION COMPLETED '
WRITE(iunit2,*) '        --> STARTING TIMESTEPPING '
WRITE(iunit2,*)

CLOSE(iunit2,STATUS='keep')

IF (GIMRT) THEN
  
  WRITE(*,*)
  WRITE(*,*) ' ---> RUNNING IN:  GLOBAL IMPLICIT REACTION AND TRANSPORT (GIMRT) MODE'
  WRITE(*,*)
  IF (nz > 1) THEN
    WRITE(*,*)
    WRITE(*,*) ' No Z discretization allowed in GLOBAL IMPLICIT mode at present time'
    WRITE(*,*)
!*******************petsc addition to cleanup *******************************************
    call PetscFinalize(ierr)
!******************end petsc addition for cleanup on error ******************************
    READ(*,*)
    STOP
	
  END IF
  
ELSE
  
  WRITE(*,*)
  WRITE(*,*) ' ---> RUNNING IN:  TIME SPLITTING REACTION AND TRANSPORT (OS3D) MODE'
  WRITE(*,*) ' Except this needs fixing, so check back after ionic strength and'      
  WRITE(*,*) '   activity of water are updated there'
  WRITE(*,*)
  STOP
  
END IF

WRITE(*,*)
WRITE(*,*) ' Grid cells in X direction (NX): ', nx
WRITE(*,*) ' Grid cells in Y direction (NY): ', ny
WRITE(*,*) ' Grid cells in Z direction (NZ): ', nz
WRITE(*,*)

!!!  *******************************************************
!!!  ***************  FLOW CALCULATION  ********************

IF (CalculateFlow) THEN

  
  SteadyFlow = .FALSE.
  
  CALL CrunchPETScInitializePressure(nx,ny,nz,userP,ierr,xvecP,bvecP,amatP)
  
  dtflow = delt

  harx = 0.0
  hary = 0.0
  harz = 0.0

  CALL GlobalDensity(nx,ny,nz)

!!  Initialize ghost cells (to be eliminated once ghost cells are removed)
  ro(0,:,:) = ro(1,:,:)
  ro(nx+1,:,:) = ro(nx,:,:)
  ro(:,0,:) = ro(:,1,:)
  ro(:,ny+1,:) = ro(:,ny,:)
  ro(:,:,0) = ro(:,:,1)
  ro(:,:,nz+1) = ro(:,:,nz)

  ! compute face permeability (harx etc.) based on distance-weighted harmonic mean
  CALL harmonic(nx,ny,nz)

  MaxPermeabilityX = MAXVAL(harx)
  MaxPermeabilityY = MAXVAL(hary)
  WRITE(*,*) ' Maximum X permeability: ',MaxPermeabilityX
  WRITE(*,*) ' Maximum Y permeability: ',MaxPermeabilityY
  WRITE(*,*)

  MinPermeabilityX = MINVAL(harx)
  MinPermeabilityY = MINVAL(hary)
  WRITE(*,*) ' Minimum X permeability: ',MinPermeabilityX
  WRITE(*,*) ' Minimum Y permeability: ',MinPermeabilityY
  WRITE(*,*)

  IF (InitializeHydrostatic) THEN
    WRITE(*,*) ' Initializing to hydrostatic'
    DO jz = 1,nz
      DO jy = 1,ny
		    DO jx = 1,nx
          harx(jx,jy,jz) = (1.0d0 - COSD(x_angle))*1.0D-22 + COSD(x_angle)*harx(jx,jy,jz)
          hary(jx,jy,jz) = (1.0d0 - COSD(y_angle))*1.0D-22 + COSD(y_angle)*hary(jx,jy,jz)
          harz(jx,jy,jz) = (1.0d0 - COSD(z_angle))*1.0D-22 + COSD(z_angle)*harz(jx,jy,jz)
        END DO
      END DO
    END DO
  END IF
  
 ! Edit by Toshiyuki Bandai 2024 Oct
 ! Because the 1D Richards solver by Toshiyuki Bandai does not use PETSc, we need to diverge here
 
 initial_flow_solver_if: IF (Richards) THEN
 ! ******************************************************************
 ! Steady-state Richards solver by Toshiyuki Bandai, 2023 May

   steady_Richards: IF (Richards_Options%is_steady) THEN
   ! solve the 1D state-state Richards equation
     WRITE(*,*) ' Solves the steady-state Richards equation to obtain the the initial condition. '
     Richards_BCs_pointer => Richards_BCs_steady
     CALL RichardsSolve(nx, ny, nz, delt)
     Richards_Options%is_steady = .FALSE.
     Richards_BCs_pointer => Richards_BCs

   ELSE steady_Richards
     
     WRITE(*,*) ' Steady-state Richards equation was not used to obtain the initial condition. '
     ! compute water flux from the initial condition and the boundary conditions at t = 0
     CALL RichardsFlux(nx, ny, nz)
     
   END IF steady_Richards

 ! End of edit by Toshiyuki Bandai, 2024 Oct
 ! ******************************************************************
 ELSE initial_flow_solver_if

  atolksp = 1.D-50
  rtolksp = GIMRTRTOLKSP
  rtolksp = 1.0D-25
  dtolksp = 1.0D-30

  pc%v = userP(5)
  ksp%v = userP(6)

  CALL KSPSetOperators(ksp,amatP,amatP,ierr)
  
  WRITE(*,*)
  WRITE(*,*) ' Running flow field to steady state prior to chemistry'
  WRITE(*,*)
  
  IF (NavierStokes) THEN
    CALL pressureNS(nx,ny,nz,dtflow,amatP,SteadyFlow)
  ELSE
    CALL pressure(nx,ny,nz,dtflow,amatP,SteadyFlow)
  END IF

  CALL CrunchPETScTolerances(userP,rtolksp,atolksp,dtolksp,maxitsksp,ierr)

  CALL KSPSolve(ksp,BvecP,XvecP,ierr)
  CALL KSPGetIterationNumber(ksp,itsiterate,ierr)

  WRITE(*,*)
  WRITE(*,*) ' Number of iterations for initial flow calculation = ', itsiterate
  WRITE(*,*)

  CALL KSPGetConvergedReason(ksp,reason,ierr)
  IF (reason == 2) THEN
    WRITE(*,*) ' Converged on relative tolerance in linear solver'
  END IF
  IF (reason == 3) THEN
    WRITE(*,*) ' Converged on absolute tolerance in linear solver'
  END IF
  IF (reason == 4) THEN
    WRITE(*,*) ' Converged based on iterations in linear solver'
  END IF
  IF (reason < 0) THEN
    WRITE(*,*)
    WRITE(*,*) ' Steady state flow failed to converge, but continue on '
    ierr = 0
  END IF

 
!     ***** PETSc closeout*******

  IF (ierr /= 0) then
    WRITE(*,*)
    WRITE(*,*) ' Error solving pressure equation in KSPSolve', ierr

!   ***** PETSc closeout**************
    IF (petscon) then
      call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
    END IF
    IF (CalculateFlow) then
      call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
    END IF
    IF (OS3Dpetsc) then
      call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
    END IF
    call PetscFinalize(ierr)
!   ****** PETSc closeout finished *********
    READ(*,*)
    STOP
  END IF

  FORALL (jx=1:nx, jy=1:ny, jz=1:nz)
    pres(jx,jy,jz) = XvecCrunchP((jz-1)*nx*ny + (jy-1)*nx + jx - 1)
  END FORALL
  

  IF (NavierStokes) THEN
      CALL velocalcNS(nx,ny,nz,dtflow)
  ELSE
      CALL velocalc(nx,ny,nz)
  END IF
  
  END If initial_flow_solver_if
 ! End of If construct for solvers needing PETSc or not
   
  ! **********************************************
  ! Edit by Toshiyuki Bandai, 2024 Oct.
  ! calculate saturation from volumetric water content
  IF (Richards) THEN
    
    ! get the initial water content
    
    
    DO jz = 1, nz
      DO jy = 1, ny
        DO jx = 1, nx
            satliq(jx,jy,jz) = Richards_State%theta(jx,jy,jz)/por(jx,jy,jz)
        END DO
      END DO
    END DO
           
    ! fill ghost points by zero-order extrapolation
    text = 'Liquid_Saturation'
    satliq(0,:,:) = satliq(1,:,:)
    lowX = LBOUND(satliq,1)
    lowY = LBOUND(satliq,2)
    lowZ = LBOUND(satliq,3)
    highX = UBOUND(satliq,1)
    highY = UBOUND(satliq,2)
    highZ = UBOUND(satliq,3)
    call GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,satliq,TEXT)
    
    satliqold = satliq
  
  END IF
  ! End of Edit by Toshiyuki Bandai, 2024 Oct.

!!  Check divergence of flow field

  MaxDivergence = 0.00
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        Coordinate = 'X'
        call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
        checkw = RoAveLeft*qx(jx-1,jy,jz)*dyy(jy)*dzz(jx,jy,jz)
        checke = RoAveRight*qx(jx,jy,jz)*dyy(jy)*dzz(jx,jy,jz)
        Coordinate = 'Y'
        call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
        checkn = RoAveRight*qy(jx,jy,jz)*dxx(jx)*dzz(jx,jy,jz)
        checks = RoAveLeft*qy(jx,jy-1,jz)*dxx(jx)*dzz(jx,jy,jz)
        Coordinate = 'Z'
        call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
        checkPlus = RoAveRight*qz(jx,jy,jz)*dxx(jx)*dyy(jy)
        checkMinus = RoAveLeft*qz(jx,jy,jz-1)*dxx(jx)*dyy(jy)

        pumpterm = 0.0d0

        !! 1) pump time series
        IF (pumptimeseries) THEN
          IF (npump(jx,jy,jz)>0) THEN
            qg(1,jx,jy,jz)=qgdum
          ELSE
            qg(1,jx,jy,jz)=0
          END IF

          pumpterm = pumpterm + qg(1,jx,jy,jz)

          !! 2) normal pump
          ELSEIF (wells) THEN

            DO npz = 1,npump(jx,jy,jz)
              pumpterm = pumpterm + qg(npz,jx,jy,jz)
            END DO

          END IF
        RealSum = ro(jx,jy,jz)* pumpterm + checkw+checks+checkMinus-checkn-checke-CheckPlus

        IF (DABS(RealSum) > MaxDivergence) THEN
            jxmax = jx
            jymax = jy
        END IF
        MaxDivergence = DMAX1(MaxDivergence,DABS(RealSum))
      END DO
    END DO
  END DO

  WRITE(*,*) ' Maximum divergence in flow field = ', MaxDivergence
  write(*,*) ' At grid cells: ',jxmax,jymax
  WRITE(*,*)

  SteadyFlow = .FALSE.

END IF   !!  End of steady state flow block

!!!  **********  END OF FLOW BLOCK  ************************
!!!  *******************************************************


!!  Initial mass in system
IF (ALLOCATED(InitialMass)) THEN
  DEALLOCATE(InitialMass)
END IF
ALLOCATE(InitialMass(nx,ny,nz))
    
InitialTotalMass = 0.0d0
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      
      DO ik = 1,ncomp+nspec
        ulabPrint = ulab(ik)
        IF (ulabPrint(1:3) /= 'H2O' .and. ulabPrint(1:3) /= 'HHO') THEN
         InitialTotalMass = InitialTotalMass + sp10(ik,jx,jy,jz)*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
        ENDIF
      END DO
      
    END DO
  END DO
END DO

WRITE(*,*)
WRITE(*,*) ' Initial aqueous mass in system = ', InitialTotalMass
WRITE(*,*)
WRITE(*,*)


!!  Initial allocation for OS3D and GIMRT modes

IF (OS3D) then
  call AllocateOS3D(ncomp,nspec,ngas,nrct,nexchange,nsurf,nsurf_sec,npot,neqn,nx,ny,nz)
END IF

IF (GIMRT) THEN
  call AllocateGIMRT(ncomp,nspec,ngas,nrct,nexchange,nsurf,nsurf_sec,npot,neqn,nx,ny,nz)
END IF

call AllocateALL(ncomp,nspec,ngas,nrct,nexchange,nsurf,nsurf_sec,npot,neqn,nx,ny,nz)

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
      
      DO npt = 1,npot
        LogPotential(npt,jx,jy,jz) = LogPotentialInit(npt,jinit(jx,jy,jz))
      END DO
      
      IF (Duan .OR. Duan2006) GasPressureTotal(jx,jy,jz) = GasPressureTotalInit(jinit(jx,jy,jz))
      
    END DO
  END DO
END DO


!!  Plotting for exchange and surface complexation
nplotexchange = 0
nplotsurface = 0
IF (nexchange > 0) THEN
  nplotexchange = ncomp
END IF
IF (nsurf > 0) THEN
  nplotsurface = ncomp
END IF

!***************** Initialize PETSc solver, vectors, matrix ************************
IF (GIMRT .AND. (nxyz > 1) .AND. petscon) THEN  ! set up PETSc vectors, matrix, and solver

  CALL CrunchPETScInitializeChemistry(nx,ny,nz,neqn,xvec,bvec,amatpetsc,userC,ierr)

endif   !  PETSc setup
!*****************************end PETSc setup **************************************

IF (irestart == 1) THEN

  CALL restart(time,nn,nint,nexchange,nsurf,nrct,nx,ny,nz,nstop,nstopsave,          &
               delt,dtold,tstep,deltmin,dtmaxcour,dtmax,userC,userD,userP,user,     &
               amatpetsc,amatD,amatP,bvec,xvec,bvecD,xvecD,bvecP,xvecP)

END IF


!!  Check to see if the problem is fully saturated initially and changes to unsaturated
!!  with the call to "ReadFlowField".

IF (isaturate == 1) THEN            !!  Already treated as a unsaturated problem
  AlreadyUnsaturated = .TRUE.
END IF

IF (GIMRT .AND. isaturate == 1 .AND. .NOT. AlreadyUnsaturated) THEN
  CALL AllocateGasesGIMRT(nx,ny,nz,ncomp)
END IF

IF (OS3D .AND. isaturate == 1) THEN         !!  NOTE:  No previous allocation of gases for OS3D
  CALL AllocateGasesOS3D(nx,ny,nz,ncomp)
END IF

!  **********************  START OS3D BLOCK  *********************************

IF (OS3D) THEN

    call FindMaxFlow(nx,ny,nz)
    IF (xflow .OR. yflow .OR. zflow) THEN
      call CourantStepAlt(nx,ny,nz,dtmaxcour)
    END IF

END IF

!  **********************  END OS3D BLOCK  *********************************


!  ************  OS3D diffusion (OS3Dpetsc flag)  ****************************

IF (OS3D .AND. nxyz > 1) THEN
  
  !!! If dispersion or diffusion /= 0.0, then allocate gas transport coefficients
  !!! QUESTION: Why is this using parameters for aqueous diffusion or dispersion??
  IF (alfL > 0.0 .OR. alfT > 0.0 .OR. dcoeff > 0.0 .OR. dzero > 0.0) THEN
    
    OS3Dpetsc = .TRUE.
    IF (isaturate == 1) THEN
      IF (ALLOCATED(ag)) THEN
        DEALLOCATE(ag)
      END IF
      ALLOCATE(ag(1:nx+1,ny,nz))
      IF (ALLOCATED(bg)) THEN
        DEALLOCATE(bg)
      END IF
      ALLOCATE(bg(nx,ny,nz))
      IF (ALLOCATED(cg)) THEN
        DEALLOCATE(cg)
      END IF
      ALLOCATE(cg(nx,ny,nz))
      IF (ALLOCATED(dg)) THEN
        DEALLOCATE(dg)
      END IF
      ALLOCATE(dg(nx,ny,nz))
      IF (ALLOCATED(eg)) THEN
        DEALLOCATE(eg)
      END IF
      ALLOCATE(eg(nx,ny,nz))
      IF (ALLOCATED(fg)) THEN
        DEALLOCATE(fg)
      END IF
      ALLOCATE(fg(0:nx+1,0:ny+1,nz))

      ag = 0.0
      bg = 0.0
      cg = 0.0
      dg = 0.0
      eg = 0.0
      fg = 0.0
    END IF
    
  END IF
END IF

IF ( OS3Dpetsc) THEN

  call CrunchPETScInitializeDiffusion(nx,ny,nz,xvecD,bvecD,amatD,userD,ierr)

END IF

iteration_tot = 0
nn = 0
!**********************************************
! record initial state by Toshiyuki Bandai 2024, Oct.
IF (Richards) THEN
  OPEN(unit = 10, file = 'initial_condition_Richards.tec', ACCESS='sequential',STATUS='unknown')
  WRITE(10,*) 'TITLE = "Initial condition for Richards solver" '
  WRITE(10,*) 'VARIABLES = "X" "Y" "Z" "Water Content" "Water Potential" "Liquid Saturation"'
  WRITE(10,*) 'ZONE I=', nx,  ', J=',ny, ', K=',nz, ' F=POINT'
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        WRITE(10, '(6(1X,1PE16.7))') x(jx)*OutputDistanceScale,y(jy)*OutputDistanceScale, &
              z(jz)*OutputDistanceScale,Richards_State%theta(jx,jy,jz), Richards_State%psi(jx,jy,jz), satliq(jx,jy,jz)
      END DO
    END DO
  END DO
END IF
!**********************************************

i_substep = 0
DO WHILE (nn <= nend)

  IF (nn == 0 .AND. dtflow > 1.0E-15 ) THEN
    delt = dtflow
  END IF

  i_substep = i_substep + 1

  nn = nn + 1

!!!  Call timestep-dependent re-speciation routine to check whether initial state should be updated  FLASH

  jx = 0
  jy = 1
  jz = 1
    
  roOld = ro
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        CALL density(jx,jy,jz)
      END DO
    END DO
  END DO

  IF (nexchange > 0) THEN
    CALL UpdateExchanger(nx,ny,nz,nexchange)
  END IF

!!! Deleted code from Lucien (Lucien-Delete1.txt)

  IF (CalculateFlow) THEN

    !!! Initialize ghost cells for flow
    ro(0,:,:) = ro(1,:,:)
    ro(nx+1,:,:) = ro(nx,:,:)
    ro(:,0,:) = ro(:,1,:)
    ro(:,ny+1,:) = ro(:,ny,:)
    ro(:,:,0) = ro(:,:,1)
    ro(:,:,nz+1) = ro(:,:,nz)

    IF (jpor == 1 .OR. jpor == 3) THEN
      IF (.not. CubicLaw) THEN
        CALL porperm(nx,ny,nz)
      END IF
    END IF
        
    ! Edit by Toshiyuki Bandai 2023 May
    ! Because the 1D Richards solver by Toshiyuki Bandai does not use PETSc, we need to diverge here
    flow_solver_if_time: IF (Richards) THEN
        ! ******************************************************************
          IF (Richards_Options%is_print) THEN
            WRITE(*,*) ' Solves the time-dependent Richards equation at t = ', time + delt ! get the solution at t = time + delt
          END IF
               
          ! update fluid property
          CALL RichardsUpdateFluid(t, Richards_State%xi)
          
          ! update permeability at faces
          CALL harmonic(nx,ny,nz)
          
          ! store the previous time step water content
          Richards_State%theta_prev = Richards_State%theta         
          
          ! update transient boundary condition
          DO i = 1, Richards_Base%n_bfaces
            IF (Richards_BCs_pointer(i)%is_variable) THEN
              Richards_Variable_BC_ptr => Richards_Variable_BC
              
              DO j = 1, Richards_BCs_pointer(i)%variable_BC_index - 1
                Richards_Variable_BC_ptr => Richards_Variable_BC%p
              END DO
              
              CALL interp(time + delt, Richards_Variable_BC_ptr%BC_time, Richards_Variable_BC_ptr%BC_values, Richards_BCs_pointer(i)%BC_value, size(Richards_Variable_BC_ptr%BC_time))
            END IF
          END DO
             
          ! solve the Richards equation
          CALL RichardsSolve(nx, ny, nz, delt)
          
          ! update saturation
          satliqold = satliq
          
          DO jz = 1, nz
            DO jy = 1, ny
              DO jx = 1, nx
                  satliq(jx,jy,jz) = Richards_State%theta(jx,jy,jz)/por(jx,jy,jz)
              END DO
            END DO
          END DO
    
       
          ! fill ghost points by zero-order extrapolation
          text = 'Liquid_Saturation'
          lowX = LBOUND(satliq,1)
          lowY = LBOUND(satliq,2)
          lowZ = LBOUND(satliq,3)
          highX = UBOUND(satliq,1)
          highY = UBOUND(satliq,2)
          highZ = UBOUND(satliq,3)
          call GhostCells(nx,ny,nz,lowX,lowY,lowZ,highX,highY,highZ,satliq,TEXT)
          
          ! the velocity at the boundary is forced to zero when the vector goes outward
          ! not to consider chemcial transport via evaporation

          IF (Richards_Options%evaporation_boundary) THEN
            IF (nx > 1 .AND. ny == 1 .AND. nz == 1) THEN ! one-dimensional problem
              jy = 1
              jz = 1
              
              DO i = 1, Richards_Base%n_bfaces
                IF (Richards_BCs_pointer(i)%is_atmosphere) THEN
                  IF (i == 1) THEN
                  ! left boundary
                    jx = 0
                    IF (qx(jx, jy, jz) < 0.0d0) THEN
                      qx(jx, jy, jz) = 0.0d0
                    END IF
                  ELSE
                  ! right boundary
                    jx = nx
                    IF (qx(jx, jy, jz) > 0.0d0) THEN
                      qx(jx, jy, jz) = 0.0d0
                    END IF
                  END IF
                END IF
                
              END DO
            ELSE IF (nx > 1 .AND. ny > 1 .AND. nz == 1) THEN ! two-dimensional problem
              jz = 1
              IF (Richards_Base%spatial_domain == 'regular') THEN  
                DO i = 1, Richards_Base%n_bfaces
                IF (Richards_BCs_pointer(i)%is_atmosphere) THEN
                  IF (i <= nx) THEN
                  ! bottom boundary
                    jx = i
                    jy = 0
                    IF (qy(jx, jy, jz) < 0.0d0) THEN
                      qy(jx, jy, jz) = 0.0d0
                    END IF
                  ELSE IF (i > nx .AND. i <= nx+ny) THEN
                    ! right boundary
                    jx = nx
                    jy = i - nx
                    IF (qx(jx, jy, jz) > 0.0d0) THEN
                      qx(jx, jy, jz) = 0.0d0
                    END IF
                  ELSE IF (i > nx + ny .AND. i <= 2*nx+ny) THEN  
                    ! top boundary
                    jx = 2*nx + ny - i + 1
                    jy = ny
                    IF (qy(jx, jy, jz) > 0.0d0) THEN
                      qy(jx, jy, jz) = 0.0d0
                    END IF
                  ELSE
                  ! left boundary
                    jx = 0
                    jy = 2*nx + 2*ny - i + 1
                    IF (qx(jx, jy, jz) < 0.0d0) THEN
                      qx(jx, jy, jz) = 0.0d0
                    END IF
                  END IF
                END IF
                
              END DO
                
              ELSE
                WRITE(*,*)
                WRITE(*,*) ' Currently, two-dimensional Richards solver does not support the shape ', Richards_Base%spatial_domain
                WRITE(*,*)
                READ(*,*)
                STOP
              END IF
            ELSE IF (nx > 1 .AND. ny > 1 .AND. nz > 1) THEN
              WRITE(*,*)
              WRITE(*,*) ' Currently, three-dimensional Richards solver is supported.'
              WRITE(*,*)
              READ(*,*)
              STOP
            END IF
          END IF
        ! End of edit by Toshiyuki Bandai, 2024 Oct
          
        ! ******************************************************************
        !!! Original version using PETSc solve for fully saturated flow
        ELSE flow_solver_if_time
   
          atolksp = 1.D-50
          rtolksp = GIMRTRTOLKSP
          rtolksp = 1.0D-25
          dtolksp = 1.0D-30

          pc%v = userP(5)
          ksp%v = userP(6)

          SteadyFlow = .FALSE.
          CALL KSPSetOperators(ksp,amatP,amatP,ierr)
          CALL harmonic(nx,ny,nz)

          IF (NavierStokes) THEN
            CALL pressureNS(nx,ny,nz,delt,amatP,SteadyFlow)
          ELSE
            CALL pressure(nx,ny,nz,delt,amatP,SteadyFlow)
          END IF

          CALL CrunchPETScTolerances(userP,rtolksp,atolksp,dtolksp,maxitsksp,ierr)

!!!!  To invoke direct solve
!!!!    call KSPGetPC(ksp,pc,ierr)
!!!!    call PCSetType(pc,PCLU,ierr)

        IF (SolverMethod /= 'direct') THEN
          CALL KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
        END IF

    !!!    CALL MatView(amatP, PETSC_VIEWER_STDOUT_SELF,ierr)
    !!!    CALL VecView(BvecP, PETSC_VIEWER_STDOUT_SELF,ierr)
        CALL KSPSolve(ksp,BvecP,XvecP,ierr)
    !!!    CALL KSPGetIterationNumber(ksp,itsiterate,ierr)

    !!!    CALL VecView(XvecP, PETSC_VIEWER_STDOUT_SELF,ierr)
  
        IF (ierr /= 0) then
          WRITE(*,*)
          WRITE(*,*) ' Error solving pressure equation in KSPSolve', ierr

    !     ***** PETSc closeout**************
          IF (petscon) then
            call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
          END IF
          IF (CalculateFlow) then
            call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
          END IF
          IF (OS3Dpetsc) then
            call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
          END IF
          call PetscFinalize(ierr)
    !     ****** PETSc closeout finished *********
          READ(*,*)
          STOP
        END IF
        
        FORALL (jx=1:nx, jy=1:ny, jz=1:nz)
          pres(jx,jy,jz) = XvecCrunchP((jz-1)*nx*ny + (jy-1)*nx + jx - 1)
        END FORALL

        IF (NavierStokes) THEN
            CALL velocalcNS(nx,ny,nz,delt)
        ELSE
            CALL velocalc(nx,ny,nz)
        END IF
        
        END If flow_solver_if_time
          
  END IF

!! Return here to restart time step after failure

  4000 CONTINUE

!  *********  End NUFT block within time stepping  ****************

  IF (GIMRT) THEN         !  Update dispersivity
    CALL dispersivity(nx,ny,nz)
  END IF

!  **************  OS3D BLOCK    **********************
  IF (OS3D) THEN

    maxQx = MAXVAL(DABS(qx(0:nx,1:ny,1:nz)))
    IF (maxQx > 0.0) THEN
      xflow = .TRUE.
    ELSE
      xflow = .FALSE.
    END IF

    maxQy = MAXVAL(DABS(qy(1:nx,0:ny,1:nz)))
    IF (maxQy > 0.0) THEN
      yflow = .TRUE.
    ELSE
      yflow = .FALSE.
    END IF

    maxQz = MAXVAL(DABS(qz(1:nx,1:ny,0:nz)))
    IF (maxQz > 0.0) THEN
      zflow = .TRUE.
    ELSE
      zflow = .FALSE.
    END IF

    IF (xflow .OR. yflow .OR. zflow) THEN
      call CourantStepAlt(nx,ny,nz,dtmaxcour)
      CALL dispersivity(nx,ny,nz)
    END IF

    IF (.NOT. CalculateFlow) THEN
      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            CALL density(jx,jy,jz)
          END DO
        END DO
      END DO
    END IF

    IF (OS3Dpetsc) THEN
      
      IF (spherical) THEN
        CALL coeffSphericalNew(nx,ny,nz)
      ELSE if (cylindrical) THEN
        CALL coeffCylinderDiffuse(nx,ny,nz)
      ELSE
        CALL coeffDiffuse(nx,ny,nz)
        IF (isaturate == 1) THEN
          CALL gasdiff(nx,ny,nz)
        END IF
      END IF
      
    END IF

    IF (delt > dtmaxcour .AND. dtmaxcour /= 0.0) THEN
      delt = dtmaxcour
    END IF

    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          
          CALL oldcon(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)
          CALL oldkd(ncomp,jx,jy,jz)
          IF (isaturate == 1) THEN
            CALL oldcongas(ncomp,ngas,jx,jy,jz)
          END IF
          CALL oldsurf(ncomp,nsurf,nsurf_sec,jx,jy,jz)
          
        END DO
      END DO
    END DO

    CALL xmass(nx,ny,nz,ncomp,nspec)
    xgramOld = xgram

!!! Now do the boundaries for the TVD used by OS3D
    DO i = 1,ncomp

      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = -1,0
            ctvd(jx,jy,jz) = sbnd(i,1)
          END DO
        END DO
      END DO
      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = nx+1,nx+2
            ctvd(jx,jy,jz) = sbnd(i,2)
          END DO
        END DO
      END DO
      DO jz = 1,nz
        DO jy = -1,0
          DO jx = 1,nx
            ctvd(jx,jy,jz) = sbnd(i,3)
          END DO
        END DO
      END DO
      DO jz = 1,nz
        DO jy = ny+1,ny+2
          DO jx = 1,nx
            ctvd(jx,jy,jz) = sbnd(i,4)
          END DO
        END DO
      END DO
      DO jz = -1,0
        DO jy = 1,ny
          DO jx = 1,nx
            ctvd(jx,jy,jz) = sbnd(i,5)
          END DO
        END DO
      END DO
      DO jz = nz+1,nz+2
        DO jy = 1,ny
          DO jx = 1,nx
            ctvd(jx,jy,jz) = sbnd(i,6)
          END DO
        END DO
      END DO

      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            ctvd(jx,jy,jz) = sn(i,jx,jy,jz)
          END DO
        END DO
      END DO


      CALL tvd(nx,ny,nz,delt,i)
      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            IF (activecell(jx,jy,jz) == 0) THEN
              CONTINUE
            ELSE
              sn(i,jx,jy,jz) = ctvd(jx,jy,jz)
            END IF
          END DO
        END DO
      END DO

      IF (OS3Dpetsc) THEN

        atolksp = 1.D-50
        rtolksp = 1.D-09
        dtolksp = 1.D+05

        CALL CrunchPETScTolerances(userD,rtolksp,atolksp,dtolksp,maxitsksp,ierr)

        pc%v = userD(5)
        ksp%v = userD(6)

        IF (spherical) THEN
          CALL SolveDiffuseSpherical(nx,ny,nz,nn,i,delt,userD,amatD)
        ELSE
          CALL SolveDiffuse(nx,ny,nz,nn,i,delt,userD,amatD)
        END IF

        DO jz = 1,nz
          DO jy = 1,ny
            DO jx = 1,nx
              
              j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
              IF (activecell(jx,jy,jz) == 0) THEN
                CONTINUE
              ELSE
                XvecCrunchD(j) = sn(i,jx,jy,jz)
              END IF
                
            END DO
          END DO
        END DO

        CALL KSPSetOperators(ksp,amatD,amatD,ierr)
        CALL KSPSolve(ksp,BvecD,XvecD,ierr)
        CALL KSPGetIterationNumber(ksp,itsiterate,ierr)

        IF (ierr /= 0) then
          WRITE(*,*)
          WRITE(*,*) ' Error solving diffusion equation in KSPSolve', ierr

!         ***** PETSc closeout**************
          IF (petscon) then
            call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
          END IF
          IF (CalculateFlow) then
            call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
          END IF
          IF (OS3Dpetsc) then
            call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
          END IF
          call PetscFinalize(ierr)
!       ****** PETSc closeout finished *********
          READ(*,*)
          STOP
        END IF

!         Now, update the SN array

        DO jz = 1,nz
          DO jy = 1,ny
            DO jx = 1,nx
                
              j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
              IF (activecell(jx,jy,jz) == 0) THEN
                CONTINUE
              ELSE
                sn(i,jx,jy,jz) = XvecCrunchD(j)
              END IF
                
            END DO
          END DO
        END DO

      END IF
      
    END DO   ! End of species

!!  Start of gas diffusion

    IF (OS3Dpetsc .AND. isaturate==1) THEN
      atolksp = 1.D-50
      rtolksp = 1.D-09
      dtolksp = 1.D+05

      CALL CrunchPETScTolerances(userD,rtolksp,atolksp,dtolksp,maxitsksp,ierr)

      pc%v = userD(5)
      ksp%v = userD(6)

      DO i = 1,ncomp

        CALL SolveGasDiffuse(nx,ny,nz,nn,i,delt,userD,amatD)

        DO jz = 1,nz
          DO jy = 1,ny
            DO jx = 1,nx
              j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
              IF (activecell(jx,jy,jz) == 0) THEN
                CONTINUE
              ELSE
                XvecCrunchD(j) = sgasn(i,jx,jy,jz)
              END IF
            END DO
          END DO
        END DO

        CALL KSPSetOperators(ksp,amatD,amatD,ierr)
        CALL KSPSolve(ksp,BvecD,XvecD,ierr)
        CALL KSPGetIterationNumber(ksp,itsiterate,ierr)

        IF (ierr /= 0) then
          WRITE(*,*)
          WRITE(*,*) ' Error solving diffusion equation in KSPSolve', ierr

!     ***** PETSc closeout**************
          IF (petscon) then
            CALL CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
          END IF
          IF (CalculateFlow) then
            CALL CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
          END IF
          IF (OS3Dpetsc) then
            CALL CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
          END IF
          CALL PetscFinalize(ierr)
!     ****** PETSc closeout finished *********
          READ(*,*)
          STOP
        END IF

!       Now, update the SGASN array

        DO jz = 1,nz
          DO jy = 1,ny
            DO jx = 1,nx
              
              j = (jz-1)*nx*ny + (jy-1)*nx + jx - 1
              IF (activecell(jx,jy,jz) == 0) THEN
                CONTINUE
              ELSE
                sgasn(i,jx,jy,jz) = XvecCrunchD(j)
              END IF
              
            END DO
          END DO
        END DO

      END DO   ! End of cycling through components

    END IF     ! End of gas diffusion in OS3D mode


!  Start sweep through spatial domain for reaction step

!   NOTE:  call to OLDCON not needed since total concentrations are transported
!     (SN already calculated)

    newtmax = 0
    loop4001: DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          
          CALL keqcalc2(ncomp,nrct,nspec,ngas,nsurf_sec,jx,jy,jz)
          IF (igamma == 3) THEN

            IF (Duan .OR. Duan2006) THEN
!!!              CALL gamma_co2(ncomp,nspec,ngas,jx,jy,jz)
            ELSE
!!!              CALL !!! gammaUpdated(ncomp,nspec,nsurf,nexchange,npot,jx,jy,jz,igamma)
            END IF
            CALL gammaUpdated(ncomp,nspec,nsurf,nexchange,npot,jx,jy,jz,igamma)

          END IF

          CALL AqueousToBulkConvert(jx,jy,jz,AqueousToBulk)

          IF (ncomp == 1 .AND. ulab(1) == 'Tracer') THEN
             sp(1,jx,jy,jz) = DLOG(sn(1,jx,jy,jz))
             sp10(1,jx,jy,jz) = sn(1,jx,jy,jz)
          END IF

          CALL OS3D_newton(ncomp,nspec,nkin,nrct,ngas,ikin,                  &
              nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,neqn,igamma,   &
              delt,corrmax,jx,jy,jz,iterat,icvg,nx,ny,nz,time,AqueousToBulk)

          CALL SpeciesLocal(ncomp,nspec,jx,jy,jz)
          CALL totconc(ncomp,nspec,jx,jy,jz)

          IF (iterat > newtmax) THEN
            newtmax = iterat
          END IF

          IF (icvg == 1) THEN
            EXIT loop4001
          END IF

        END DO
      END DO
    END DO loop4001


    IF (icvg == 1) THEN
      
      ddtold = delt
      delt = delt/10.0
      itskip = 1

      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            
            DO i = 1,ncomp+nspec
              sp(i,jx,jy,jz) = spold(i,jx,jy,jz)
              sp10(i,jx,jy,jz) = EXP(sp(i,jx,jy,jz))
            END DO
            DO ix = 1,nexchange
              spex(ix,jx,jy,jz) = spexold(ix,jx,jy,jz)
            END DO
            DO is = 1,nsurf+nsurf_sec
              spsurf(is,jx,jy,jz) = spsurfold(is,jx,jy,jz)
              spsurf10(is,jx,jy,jz) = EXP(spsurfold(is,jx,jy,jz))
            END DO
            CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
            
          END DO
        END DO
      END DO

      IF (igamma == 2) THEN
        DO jz = 1,nz
          DO jy = 1,ny
            DO jx = 1,nx
              
!!!            IF (Duan .OR. Duan2006) THEN
!!!              CALL gamma_co2(ncomp,nspec,ngas,jx,jy,jz)
!!!            ELSE
!!!              CALL gammaUpdated(ncomp,nspec,nsurf,nexchange,npot,jx,jy,jz,igamma)
!!!            END IF
              
            CALL gammaUpdated(ncomp,nspec,nsurf,nexchange,npot,jx,jy,jz,igamma)
            
            END DO
          END DO
        END DO
      END IF

      WRITE(*,*)
      WRITE(*,*) '***** NO CONVERGENCE OF NEWTON ITERATIONS IN OS3D *****'
      WRITE(*,*) '                REDUCING TIME STEP'
      WRITE(*,5086) ddtold*OutputTimeScale
      WRITE(*,5085) delt*OutputTimeScale
      WRITE(*,*)
      GO TO 4000   ! Loop back to start the time step over

    ELSE
      iterat = newtmax
    END IF   !  End of IF block for case where no Newton convergence


  END IF    ! End of OS3D block

!!!  *************   END OF OS3D BLOCK  **********************
!!!  ********************************************************

!!!  ********************************************************
!!!  *************   START GIMRT BLOCK  *********************

6000 IF (GIMRT) THEN
  

       
    i_substep = 1
    n_substep = 1
    dt_GIMRT = delt

    IF (i_substep == n_substep) THEN
      
        ! Invoke GIMRT calculation
        jz = 1
    !           Calculate finite difference coefficients
        CALL dispersivity(nx,ny,nz)

        IF (TortuosityOption /= 'none') THEN
          CALL CalculateTortuosity(nx,ny,nz)
        END IF

        !!  Diffusion block for GIMRT
        IF (species_diffusion) THEN 
          
          IF (spherical) THEN
            IF (xflow .OR. yflow .OR. zflow) THEN
              WRITE(*,*)
              WRITE(*,*) ' Spherical coordinates for diffusion problems only at present'
              WRITE(*,*) ' Aborting run'
              WRITE(*,*)
              READ(*,*)
              STOP
            ELSE
              CALL coeffSphericalNew_d(nx,ny,nz,ncomp,nspec)
            END IF     
          ELSE IF (cylindrical) THEN
            CALL coeffCylinder_d(nx,ny,nz,ncomp,nspec) 
          ELSE
            CALL coeff_d(nx,ny,nz,ncomp,nspec)    
          END IF
          
        ELSE
          
          IF (cylindrical) THEN
            CALL coeffCylinder(nx,ny,nz)
          ELSE IF (spherical) THEN
            CALL coeffSphericalNew(nx,ny,nz)
          ELSE
            CALL coeff(nx,ny,nz)
          END IF
          
        END IF

        IF (isaturate == 1) THEN
          IF (cylindrical) THEN
    !!        CALL gasdiffCylinder(nx,ny,nz)
            CALL gascoeffCylinder(nx,ny,nz)
          ELSE
            CALL gasdiff(nx,ny,nz)
          END IF
        END IF

        IF (ierode == 1) THEN
          CALL erosion(nx,ny,nz)
        END IF

        jz = 1
        DO jy = 1,ny
          DO jx = 1,nx
            CALL keqcalc2(ncomp,nrct,nspec,ngas,nsurf_sec,jx,jy,jz)
          END DO
        END DO

        IF (igamma == 3 .or. igamma == 2 .or. igamma == 0) THEN
          jz = 1
          DO jy = 1,ny
            DO jx = 1,nx
              
                IF (Duan .OR. Duan2006) THEN                 
!!!                  CALL gamma_co2(ncomp,nspec,ngas,jx,jy,jz)
                 ELSE
!!!                   CALL gammaUpdated(ncomp,nspec,nsurf,nexchange,npot,jx,jy,jz,igamma)
                 END IF
                  
                CALL gammaUpdated(ncomp,nspec,nsurf,nexchange,npot,jx,jy,jz,igamma)
                if (jx==8) then
                  continue
                end if
                if (jx==401) then
                  continue
                end if
            END DO
          END DO
          
        END IF
        
        jz = 1
        DO jy = 1,ny
          DO jx = 1,nx
            CALL oldcon(ncomp,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)
            IF (isaturate == 1) THEN
              CALL oldcongas(ncomp,ngas,jx,jy,jz)
            END IF
            CALL oldsurf(ncomp,nsurf,nsurf_sec,jx,jy,jz)
          END DO
        END DO

    ! **************  START NEWTON LOOP  *******************

        5000     NE = 0
        icvg = 1
        iterat = 0 ! number of Newton iterations

    newtonloop:  DO WHILE (icvg == 1 .AND. iterat <= 9)
          NE = NE + 1
          iterat = iterat + 1
          
          CALL species(ncomp,nspec,nsurf,nexchange,npot,nx,ny,nz)
		      CALL jacobian(ncomp,nspec,nx,ny,nz)								   
            
          jz = 1
          DO jy = 1,ny
            DO jx = 1,nx
                 
              IF (Duan .OR. Duan2006) THEN
!!!                    CALL gamma_co2(ncomp,nspec,ngas,jx,jy,jz)
              ELSE
!!!                    CALL gammaUpdated(ncomp,nspec,nsurf,nexchange,npot,jx,jy,jz,igamma)
              END IF
                  
              if (igamma == 2) then
                CALL gammaUpdated(ncomp,nspec,nsurf,nexchange,npot,jx,jy,jz,igamma)
              end if
 
              CALL totconc(ncomp,nspec,jx,jy,jz)
                  
            END DO
          END DO 
                          
          IF (ierode == 1) THEN
            CALL SurfaceComplex(ncomp,nsurf,nsurf_sec,nx,ny,nz)
            CALL jacsurf(ncomp,nsurf,nsurf_sec,nx,ny,nz)
          END IF

          jz = 1
          DO jy = 1,ny
            DO jx = 1,nx

              IF (isaturate == 1) THEN
                CALL gases(ncomp,ngas,jx,jy,jz)
              END IF
              IF (ierode == 1) THEN
                CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
              END IF
              IF (species_diffusion) THEN
                CALL jacobian_plus(ncomp,nspec,jx,jy,jz)
              END IF
              IF (isaturate == 1) THEN
                CALL jacgas(ncomp,ngas,jx,jy,jz)
              END IF
              IF (ierode == 1) THEN
                CALL jac_exchange(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,neqn,jx,jy,jz)
              END IF
              
            END DO
          END DO

          jz = 1
          DO jy = 1,ny
            DO jx = 1,nx
              
              IF (species_diffusion) THEN
                CALL totconc_plus(ncomp,nspec,jx,jy,jz)
              ELSE
                CALL totconc(ncomp,nspec,jx,jy,jz)
              END IF
              IF (isaturate == 1) THEN
                CALL totgas(ncomp,nspec,ngas,jx,jy,jz)
              END IF
              IF (ierode == 1) THEN
                CALL totexchange(ncomp,nexchange,nexch_sec,nsurf,nsurf_sec,jx,jy,jz)
                CALL totsurf(ncomp,nsurf,nsurf_sec,jx,jy,jz)
              END IF
              
            END DO  ! end of J loop
          END DO

          CALL xmass(nx,ny,nz,ncomp,nspec)

    !************************ Start PETSC changes for matrix fill ********************
          if(petscon) then
            call MatZeroEntries(amatpetsc,ierr)
          endif

          CALL AssembleGlobal(nx,ny,nz,ncomp,nspec,nkin,nrct,ngas,ikin,                   &
             nexchange,nexch_sec,nsurf,nsurf_sec,npot,ndecay,nn,dt_GIMRT,time, &
             userC,amatpetsc,nBoundaryConditionZone)

          if ((nn.eq.1) .and. (ne.eq.1) .and. petscon) then
    !!          call MatSetOption(amatpetsc,MAT_NO_NEW_NONZERO_LOCATIONS,ierr)
          endif
    !**************  End PETSC changes *********************************************


          IF (nxyz == 1) THEN

            bb = -fxx

    !!        CALL ludcmp90(aaa,indd,det,neqn)
    !!        CALL lubksb90(aaa,indd,bb,neqn)

            CALL dgetrf(neqn,neqn,aaa,neqn,indd,info)
            CALL dgetrs(trans,neqn,ione,aaa,neqn,indd,bb,neqn,info)

            xn = bb

          ELSE       !  One-dimensional case

            IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1 .or. Switcheroo) THEN

              IF (switcheroo) then
                write(*,*) ' Using Hindmarsh solver with PETSc constructed arrays'
              end if

              DO jx = 1,nx
                j = jx
                DO i = 1,neqn
                  ind = (j-1)*(neqn) + i
                  yh(i,jx) = -fxx(ind)
                END DO
                
              END DO

              CALL decbt90(neqn,nx,ier)
              CALL solbt90(neqn,nx)

            ELSE IF (ihindmarsh == 2) THEN

    !         ***** PETSc closeout**************
              IF (petscon) then
                call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
              END IF
              IF (CalculateFlow) then
                call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
              END IF
              IF (OS3Dpetsc) then
                call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
              END IF
              call PetscFinalize(ierr)
    !         ****** PETSc closeout finished *********
              READ(*,*)
              STOP

            ELSE         !  Run PETSc

              fxx = -fxx

    !************* Invoke PETSc solver *****************************************************
              if (petscon) then

                atolksp = 1.D-50
    !!            rtolksp = 1.D-09
                rtolksp = GIMRTRTOLKSP
                dtolksp = 1.D+05

                pc%v = userC(5)
                ksp%v = userC(6)

                call KSPSetOperators(ksp,amatpetsc,amatpetsc,ierr)
                CALL GIMRTCrunchPETScTolerances(userC,rtolksp,atolksp,dtolksp,maxitsksp,ierr)

                call KSPSolve(ksp,bvec,xvec,ierr)
                CALL KSPGetIterationNumber(ksp,itsiterate,ierr)

                if(ierr.ne.0) then
                    icvg = 1
                    write(*,*) ' KSPSolve error, reducing timestep,ierr=',ierr
                    exit newtonloop
                endif
              endif

    !      xn now contains the solution to amatpetsc * xvec = bvec
    !      unless the solver failed in which case the time step will be halved
    !******************  PETSc solver exit **************************************************

            END IF

          END IF

          IF (nxyz == nx .AND. ihindmarsh == 1 .AND. nxyz /= 1 .or. Switcheroo) THEN

            errmax = 0.0
            jz = 1
            DO jy = 1,ny
              DO jx = 1,nx
                j = (jy-1)*nx+jx
                
                DO i = 1,ncomp
                  IF (ulab(i) == 'O2(aq)') THEN
                    MaximumCorrection = corrmax
                  ELSE
                    MaximumCorrection = 2.0
                  END IF
                  IF (ABS(yh(i,jx)) > errmax) THEN
                    errmax = ABS(yh(i,jx))
                  END IF
                  IF (ABS(yh(i,jx)) > MaximumCorrection) THEN
                    yh(i,jx) = SIGN(MaximumCorrection,yh(i,jx))
                  ELSE
                    CONTINUE
                  END IF
                  sp(i,jx,jy,jz) = sp(i,jx,jy,jz) + yh(i,jx)
                  sp10(i,jx,jy,jz) = DEXP(sp(i,jx,jy,jz))
                END DO

                DO ix = 1,nexchange
                  IF (DABS(yh(ix+ncomp,jx)) > errmax) THEN
                    errmax = DABS(yh(ix+ncomp,jx))
                  END IF
                  IF (DABS(yh(ix+ncomp,jx)) > MaximumCorrection) THEN
                    yh(ix+ncomp,jx) = SIGN(MaximumCorrection,yh(ix+ncomp,jx))
                  ELSE
                    CONTINUE
                  END IF
                  spex(ix,jx,jy,jz) = spex(ix,jx,jy,jz) + yh(ix+ncomp,jx)
                  spex10(ix,jx,jy,jz) = DEXP(spex(ix,jx,jy,jz))
                END DO

                DO is = 1,nsurf
                  IF (DABS(yh(is+ncomp+nexchange,jx)) > errmax) THEN
                    errmax = DABS(yh(is+ncomp+nexchange,jx))
                  END IF
                  IF (DABS(yh(is+ncomp+nexchange,jx)) > MaximumCorrection) THEN
                    yh(is+ncomp+nexchange,jx) = SIGN(MaximumCorrection,yh(is+ncomp+nexchange,jx))
                  ELSE
                    CONTINUE
                  END IF
                  spsurf(is,jx,jy,jz) = spsurf(is,jx,jy,jz) + yh(is+ncomp+nexchange,jx)
                  spsurf10(is,jx,jy,jz) = DEXP(spsurf(is,jx,jy,jz))
                END DO

                DO npt = 1,npot
                  IF (DABS(yh(npt+ncomp+nexchange+nsurf,jx)) > 0.9d0) THEN
                    yh(npt+ncomp+nexchange+nsurf,jx) = SIGN(0.9d0,yh(npt+ncomp+nexchange+nsurf,jx))
                  ELSE
                    CONTINUE
                  END IF
                  LogPotential(npt,jx,jy,jz) = LogPotential(npt,jx,jy,jz) +   &
                      yh(npt+ncomp+nexchange+nsurf,jx)
                END DO
                
          !!!   Update ionic strength
                
                ind = ncomp + nexchange + nsurf + npot + 1 
                IF (DABS(yh(ind,jx)) > MaximumCorrection) THEN
                  yh(ind,jx) = SIGN( MaximumCorrection,yh(ind,jx) )
                ELSE
                  CONTINUE
                END IF
                
                IF ( sion(jx,jy,jz) == 0.0d0 ) THEN
                  sion(jx,jy,jz) = 1.0D-06
                ELSE
                  sion(jx,jy,jz) = EXP( LOG( sion(jx,jy,jz) ) + yh(ind,jx) )
                END IF
                
          !!!   Update lngammawater
                
                ind = ncomp + nexchange + nsurf + npot + 1 + 1
                IF (DABS(yh(ind,jx) ) > MaximumCorrection) THEN
                  yh(ind,jx) = SIGN( MaximumCorrection,yh(ind,jx) )
                ELSE
                  CONTINUE
                END IF
                
                lngammawater(jx,jy,jz) = lngammawater(jx,jy,jz) + yh(ind,jx)           
                
              END DO
            END DO

          ELSE

            jz = 1
            errmax = 0.0D0
            DO jy = 1,ny
              DO jx = 1,nx      
                j = (jy-1)*nx+jx
                
                DO i = 1,ncomp
                  IF (ulab(i) == "O2(aq)") THEN
                    MaximumCorrection = corrmax
                  ELSE
                    MaximumCorrection = 2.0d0
                  END IF
                  ind = (j-1)*(neqn) + i
                  IF (DABS(xn(ind)) > errmax) THEN
                    errmax = DABS(xn(ind))
                  END IF
                  IF (DABS(xn(ind)) > MaximumCorrection) THEN
                    xn(ind) = SIGN(MaximumCorrection,xn(ind))
                  ELSE
                    CONTINUE
                  END IF
                  sp(i,jx,jy,jz) = sp(i,jx,jy,jz) + xn(ind)
                  sp10(i,jx,jy,jz) = DEXP(sp(i,jx,jy,jz))
                END DO

                DO ix = 1,nexchange
                  ind = (j-1)*(neqn) + ix+ncomp
                  IF (DABS(xn(ind)) > errmax) THEN
                    errmax = DABS(xn(ind))
                  END IF
                  IF (DABS(xn(ind)) > MaximumCorrection) THEN
                    xn(ind) = SIGN(MaximumCorrection,xn(ind))
                  ELSE
                    CONTINUE
                  END IF
                  spex(ix,jx,jy,jz) = spex(ix,jx,jy,jz) + xn(ind)
                  spex10(ix,jx,jy,jz) = DEXP(spex(ix,jx,jy,jz))
                END DO

                DO is = 1,nsurf
                  ind = (j-1)*(neqn) + is+ncomp+nexchange
                  IF (DABS(xn(ind)) > errmax) THEN
                    errmax = DABS(xn(ind))
                  END IF
                  IF (DABS(xn(ind)) > MaximumCorrection) THEN
                    xn(ind) = SIGN(MaximumCorrection,xn(ind))
                  ELSE
                    CONTINUE
                  END IF
                  spsurf(is,jx,jy,jz) = spsurf(is,jx,jy,jz) + xn(ind)
                  spsurf10(is,jx,jy,jz) = DEXP(spsurf(is,jx,jy,jz))
                END DO

                DO npt = 1,npot
                  ind = (j-1)*(neqn) + npt+ncomp+nexchange+nsurf
                 IF (DABS(xn(ind)) > errmax) THEN
                    errmax = DABS(xn(ind))
                 END IF
                  LogPotential(npt,jx,jy,jz) = LogPotential(npt,jx,jy,jz) + xn(ind)
                END DO
                
                !!!   Update ionic strength
                
                ind = ncomp + nexchange + nsurf + npot + 1
                IF (DABS(xn(ind)) > MaximumCorrection) THEN
                  xn(ind) = SIGN(MaximumCorrection,xn(ind))
                ELSE
                  CONTINUE
                END IF
                
                IF ( sion(jx,jy,jz) == 0.d0 .AND. xn(ind) == 0.0d0) THEN
                  write(*,*)
                  write(*,*) ' Ionic strength cannot be zero '
                  write(*,*)
                  stop
                END IF
                  
                sion(jx,jy,jz) = EXP( LOG( sion(jx,jy,jz) ) + xn(ind) )
                
          !!!   Update lngammawater
                
                ind = ncomp + nexchange + nsurf + npot + 1 + 1
                IF (DABS(xn(ind)) > MaximumCorrection) THEN
                  xn(ind) = SIGN(MaximumCorrection,xn(ind))
                ELSE
                  CONTINUE
                END IF
                
                lngammawater(jx,jy,jz) = lngammawater(jx,jy,jz) + xn(ind) 
                
              END DO
            END DO

          END IF

          icvg = 0
          jz = 1
          DO jy = 1,ny
            DO jx = 1,nx
              j = (jy-1)*nx+jx
              
              DO i = 1,ncomp
                ind = (j-1)*(neqn) + i
                IF (ResidualTolerance /= 0.0d0) THEN
                  tolmax = ResidualTolerance
                ELSE
                  tolmax = atol
                END IF
                
                IF (SaltCreep) THEN
                  IF (DABS(fxx(ind)) > tolmax) THEN
                    icvg = 1
                  END IF
                ELSE
                  IF (DABS(dt_GIMRT*fxx(ind)) > tolmax) THEN
!!!                IF (DABS(fxx(ind)) > tolmax) THEN
                    icvg = 1
                  END IF
                END IF
                
              END DO
              
              DO ix = 1,nexchange
                tolmax = 1.e-10
                ind = (j-1)*(neqn) + ix+ncomp
                IF (DABS(dt_GIMRT*fxx(ind)) > tolmax) THEN
                    icvg = 1
                END IF
              END DO
              
              DO is = 1,nsurf
                tolmax = atol
                ind = (j-1)*(neqn) + is+ncomp+nexchange
                IF (DABS(dt_GIMRT*fxx(ind)) > tolmax) THEN
                    icvg = 1
                END IF
              END DO
              
              DO npt = 1,npot
                tolmax = 1.e-10
                ind = (j-1)*(neqn) + npt+ncomp+nexchange+nsurf
                IF (DABS(fxx(ind)) > tolmax) THEN
                    icvg = 1
                END IF
              END DO
              
            END DO
          END DO

          IF (ABS(errmax) < 1.e-15) THEN
            icvg = 0
          END IF

    !! If electrochemical migration is considered, always carry out at least two Newton iterations
          IF (species_diffusion .AND. iterat<=2) THEN
            icvg = 1
          ELSE IF (iterat <= 1) THEN
            icvg = 1
          ELSE
            continue
          END IF

    END DO  newtonloop    ! end of Newton iteration loop

    !  Halve the timestep if convergence not achieved

    5017 CONTINUE
        IF (icvg == 1) THEN
          MaxFx = MaxLoc(fxx)
    !!      MaxValFx = MaxVal(fxx)

          ddtold = dt_GIMRT
          dt_GIMRT = dt_GIMRT/10.0
          itskip = 1

          jz = 1
          DO jy = 1,ny
            DO jx = 1,nx
              
              DO i = 1,ncomp+nspec
                sp(i,jx,jy,jz) = spold(i,jx,jy,jz)
                sp10(i,jx,jy,jz) = DEXP(sp(i,jx,jy,jz))
              END DO
              DO ix = 1,nexchange+nexch_sec
                spex(ix,jx,jy,jz) = spexold(ix,jx,jy,jz)
              END DO
              DO nex = 1,nexch_sec
                spex10(nex+nexchange,jx,jy,jz) = DEXP(spex(nex+nexchange,jx,jy,jz))
              END DO
              DO is = 1,nsurf+nsurf_sec
                spsurf(is,jx,jy,jz) = spsurfold(is,jx,jy,jz)
                spsurf10(is,jx,jy,jz) = DEXP(spsurf(is,jx,jy,jz))
              END DO
              IF (ierode == 1) THEN
                CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
              END IF
              CALL reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,dt_GIMRT,time)
              IF (isaturate == 1) THEN
                CALL gases(ncomp,ngas,jx,jy,jz)
              END IF
              
            END DO
          END DO

          jz = 1
          DO jy = 1,ny
            DO jx = 1,nx
              spno2(jx,jy,jz) = spnno2(jx,jy,jz)
            END DO
          END DO

          WRITE(*,*)
          WRITE(*,*) '***** NO CONVERGENCE OF NEWTON ITERATIONS IN GIMRT *****'
          WRITE(*,*) '                REDUCING TIME STEP'
          WRITE(*,5086) ddtold*OutputTimeScale
          WRITE(*,5085) dt_GIMRT*OutputTimeScale
          WRITE(*,*)
          GO TO 5000   ! Loop back to start the time step over

        END IF   !  End of IF block for case where no Newton convergence

    !  Correct for change in water concentration

        IF (nn == 1) THEN
          IF (ALLOCATED(lrow)) DEALLOCATE(lrow)
          IF (ALLOCATED(levptr)) DEALLOCATE(levptr)
        END IF

        CALL species(ncomp,nspec,nsurf,nexchange,npot,nx,ny,nz)
        CALL SurfaceComplex(ncomp,nsurf,nsurf_sec,nx,ny,nz)

        jz = 1
        DO jy = 1,ny
          DO jx = 1,nx
            CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
    !!        CALL reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt)
          END DO
        END DO

        dt_GIMRT = 0.0d0
        i_substep = 0
    END IF ! end i_substep == n_substep

END IF    !  END OF GIMRT NEWTON LOOP

!!!  *********  End of NEWTON LOOP  ********************
!!!  ***************************************************
       
!! For greater accuracy, update reaction rate
DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx
!!      CALL totconc(ncomp,nspec,jx,jy,jz))
      CALL reaction(ncomp,nkin,nrct,nspec,nexchange,nsurf,ndecay,jx,jy,jz,delt,time)
    END DO
  END DO
END DO

time = time + delt

!*************************************************************
! Edit by Lucien Stolze, June 2023
! Addition of a wall time
IF (time > 1e-3) THEN
  IF (walltime) THEN
    call CPU_TIME(PrintSeconds)  
    IF ((PrintSeconds/60.0d0)>wall_t) THEN
      write(*,*)
      write(*,*) 'WALLTIME REACHED'
      write(*,*)
      STOP
    END IF
  END IF
END IF
  !*************************************************************
! end of Edit by Lucien Stolze, June 2023

  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        CALL decay(ncomp,ikin,delt,jx,jy,jz)
      END DO
    END DO
  END DO

!  Do not update mineral volume fractions in reaction path
!  mode (ipath = 1 and nxyz = 1)
  
  IF (ipath == 1 .AND. nxyz == 1) THEN
    CONTINUE
  ELSE
    CALL mineral_update(nx,ny,nz,nrct,delt,dtnewest,ineg,jpor,deltmin)
    IF (FractureNetwork .and. CubicLaw) THEN
      call rmesh51(nx,ny)
    END IF
  END IF

!  NOTE:  Will have to change some routines below for OS3D

  IF (ineg == 1) THEN
    
    itskip = 1
    time = time - delt

    delt = dtnewest

    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
          DO ik = 1,ncomp
            sp(ik,jx,jy,jz) = spold(ik,jx,jy,jz)
            sp10(ik,jx,jy,jz) = EXP(sp(ik,jx,jy,jz))
          END DO
          DO ix = 1,nexchange
            spex(ix,jx,jy,jz) = spexold(ix,jx,jy,jz)
          END DO
          DO is = 1,nsurf
            spsurf(is,jx,jy,jz) = spsurfold(is,jx,jy,jz)
          END DO
          CALL exchange(ncomp,nexchange,nexch_sec,jx,jy,jz)
          IF (isaturate == 1) THEN
            CALL gases(ncomp,ngas,jx,jy,jz)
          END IF
        END DO
      END DO
    END DO

    CALL species(ncomp,nspec,nsurf,nexchange,npot,nx,ny,nz)
    CALL SurfaceComplex(ncomp,nsurf,nsurf_sec,nx,ny,nz)

    GO TO 6000
    
  END IF

!  **********************  START GIMRT EROSION BLOCK  *********************************

  IF (GIMRT .AND. ierode == 1) THEN

!         Advect the minerals via burial or erosion
!           (positive for burial, negative for erosion)

    !!!  ******* Minerals *******************************

    rinv = 1.0/delt

    DO k=1,nkin

      jz = 1
      DO jy = 1,ny
        DO jx = -1,0
          ctvd(jx,jy,jz) = volb(k,1)
        END DO
      END DO
      jz = 1
      DO jy = 1,ny
        DO jx = nx+1,nx+2
          ctvd(jx,jy,jz) = volb(k,2)
        END DO
      END DO
      jz = 1
      DO jy = -1,0
        DO jx = 1,nx
          ctvd(jx,jy,jz) = volb(k,3)
        END DO
      END DO
      jz = 1
      DO jy = ny+1,ny+2
        DO jx = 1,nx
          ctvd(jx,jy,jz) = volb(k,4)
        END DO
      END DO

      jz = 1
      DO jy = 1,ny
        DO jx = 1,nx
          ctvd(jx,jy,jz) = volfx(k,jx,jy,jz)
        END DO
      END DO

      IF (nx > 1 .AND. ny == 1) THEN

        DO jx = 1,nx
          IF (jx == 1) THEN
            dxe = 0.5*(dxx(jx)+dxx(jx+1))
            dxw = 0.5*dxx(1)
          ELSE IF (jx == nx) THEN
            dxw = 0.5*(dxx(jx)+dxx(jx-1))
            dxe = 0.5*dxx(nx)
          ELSE
            dxe = 0.5*(dxx(jx)+dxx(jx+1))
            dxw = 0.5*(dxx(jx)+dxx(jx-1))
          END IF
          IF (SolidburyX(jx) > 0.0) THEN      ! Burial case
            aabur(jx) = -SolidburyX(jx-1)
            bbbur(jx) = SolidburyX(jx-1) + rinv*dxx(jx)
            ccbur(jx) = 0.0
            uubur(jx) = ctvd(jx,1,1)
            rrbur(jx) = ctvd(jx,1,1)*rinv*dxx(jx)
          ELSE                         ! Erosion case
            aabur(jx) = 0.0
            bbbur(jx) = -SolidburyX(jx) + rinv*dxx(jx)
            ccbur(jx) = SolidburyX(jx)
            uubur(jx) = ctvd(jx,1,1)
            rrbur(jx) = ctvd(jx,1,1)*rinv*dxx(jx)
          END IF
        END DO

        IF (SolidBuryX(1) > 0.0) THEN
          rrbur(1) = rrbur(1) - aabur(1)*volb(k,1)
        END IF
        IF (SolidBuryX(nx) < 0.0) THEN
          rrbur(nx) = rrbur(nx) - ccbur(nx)*volb(k,2)
        END IF

        CALL tridag_ser(aabur,bbbur,ccbur,rrbur,uubur)

        DO jx = 1,nx
          volfx(k,jx,1,1) = uubur(jx)
        END DO

      END IF


    END DO   ! End of Mineral loop for erosion/burial

    rinv = 1.0/delt
    jz = 1
    jy = 1

    rrbur = 0.0d0

    DO k = 1,nkin

      DO jx = 1,nx
        ctvd(jx,jy,jz) = specificByGrid(k,jx,1,1)
      END DO

      DO jx = 1,nx
        IF (jx == 1) THEN
          dxe = 0.5*(dxx(jx)+dxx(jx+1))
          dxw = 0.5*dxx(1)
        ELSE IF (jx == nx) THEN
          dxw = 0.5*(dxx(jx)+dxx(jx-1))
          dxe = 0.5*dxx(nx)
        ELSE
          dxe = 0.5*(dxx(jx)+dxx(jx+1))
          dxw = 0.5*(dxx(jx)+dxx(jx-1))
        END IF

        IF (SolidburyX(jx) > 0.0) THEN      ! Burial case

          aabur(jx) = -SolidburyX(jx-1)
          bbbur(jx) = SolidburyX(jx-1) + rinv*dxx(jx)
          ccbur(jx) = 0.0
          uubur(jx) = ctvd(jx,1,1)
          rrbur(jx) = ctvd(jx,1,1)*rinv*dxx(jx)

        ELSE                         ! Erosion case

          aabur(jx) = 0.0
          bbbur(jx) = -SolidburyX(jx) + rinv*dxx(jx)
          ccbur(jx) = SolidburyX(jx)
          uubur(jx) = ctvd(jx,1,1)
          rrbur(jx) = ctvd(jx,1,1)*rinv*dxx(jx)

        END IF
      END DO

      IF (SolidBuryX(1) > 0.0) THEN
         rrbur(1) = rrbur(1) - aabur(1)*specificByGrid(k,0,1,1)
      END IF

      IF (SolidBuryX(nx) < 0.0) THEN
        rrbur(nx) = rrbur(nx) - ccbur(nx)*specificByGrid(k,nx+1,1,1)
      END IF

      CALL tridag_ser(aabur,bbbur,ccbur,rrbur,uubur)

      DO jx = 1,nx
        specificByGrid(k,jx,1,1) = uubur(jx)
      END DO

!!   ********************************************

      DO jx = 1,nx
        ctvd(jx,jy,jz) = areainByGrid(k,jx,1,1)
      END DO

      DO jx = 1,nx

        IF (jx == 1) THEN
          dxe = 0.5*(dxx(jx)+dxx(jx+1))
          dxw = 0.5*dxx(1)
        ELSE IF (jx == nx) THEN
          dxw = 0.5*(dxx(jx)+dxx(jx-1))
          dxe = 0.5*dxx(nx)
        ELSE
          dxe = 0.5*(dxx(jx)+dxx(jx+1))
          dxw = 0.5*(dxx(jx)+dxx(jx-1))
        END IF

        IF (SolidburyX(jx) > 0.0) THEN      ! Burial case

          aabur(jx) = -SolidburyX(jx-1)
          bbbur(jx) = SolidburyX(jx-1) + rinv*dxx(jx)
          ccbur(jx) = 0.0
          uubur(jx) = ctvd(jx,1,1)
          rrbur(jx) = ctvd(jx,1,1)*rinv*dxx(jx)

        ELSE                         ! Erosion case

          aabur(jx) = 0.0
          bbbur(jx) = -SolidburyX(jx) + rinv*dxx(jx)
          ccbur(jx) = SolidburyX(jx)
          uubur(jx) = ctvd(jx,1,1)
          rrbur(jx) = ctvd(jx,1,1)*rinv*dxx(jx)

        END IF
      END DO

      IF (SolidBuryX(1) > 0.0) THEN
         rrbur(1) = rrbur(1) - aabur(1)*areab(k,1)
      END IF

      IF (SolidBuryX(nx) < 0.0) THEN
        rrbur(nx) = rrbur(nx) - ccbur(nx)*areab(k,2)
      END IF

      CALL tridag_ser(aabur,bbbur,ccbur,rrbur,uubur)

      DO jx = 1,nx
        areainByGrid(k,jx,1,1) = uubur(jx)
      END DO

!!   ********************************************

      DO jx = 1,nx
        ctvd(jx,jy,jz) = volinByGrid(k,jx,1,1)
      END DO

      DO jx = 1,nx

        IF (jx == 1) THEN
          dxe = 0.5*(dxx(jx)+dxx(jx+1))
          dxw = 0.5*dxx(1)
        ELSE IF (jx == nx) THEN
          dxw = 0.5*(dxx(jx)+dxx(jx-1))
          dxe = 0.5*dxx(nx)
        ELSE
          dxe = 0.5*(dxx(jx)+dxx(jx+1))
          dxw = 0.5*(dxx(jx)+dxx(jx-1))
        END IF

        IF (SolidburyX(jx) > 0.0) THEN      ! Burial case

          aabur(jx) = -SolidburyX(jx-1)
          bbbur(jx) = SolidburyX(jx-1) + rinv*dxx(jx)
          ccbur(jx) = 0.0
          uubur(jx) = ctvd(jx,1,1)
          rrbur(jx) = ctvd(jx,1,1)*rinv*dxx(jx)

        ELSE                         ! Erosion case

          aabur(jx) = 0.0
          bbbur(jx) = -SolidburyX(jx) + rinv*dxx(jx)
          ccbur(jx) = SolidburyX(jx)
          uubur(jx) = ctvd(jx,1,1)
          rrbur(jx) = ctvd(jx,1,1)*rinv*dxx(jx)

        END IF
      END DO


      IF (SolidBuryX(1) > 0.0) THEN
        rrbur(1) = rrbur(1) - aabur(1)*volb(k,1)
      END IF

      IF (SolidBuryX(nx) < 0.0) THEN
        rrbur(nx) = rrbur(nx) - ccbur(nx)**volb(k,2)
      END IF

      CALL tridag_ser(aabur,bbbur,ccbur,rrbur,uubur)

      DO jx = 1,nx
        volinByGrid(k,jx,1,1) = uubur(jx)
      END DO

!!   ********************************************

    END DO !  Loop through minerals

    CALL porcalc(nx,ny,nz,nkin,jpor)     ! Updates porosity and surface area


  END IF   !  End of erosion/burial block (for GIMRT only)

!  **********   END OF EROSION/BURIAL BLOCK FOR MINERALS  *****************

!!    write(*,*) ' Smectite surface area = ', area(7,nx,1,1)
!!    write(*,*) ' Smectite volume       = ', volfx(7,nx,1,1)
!!    write(*,*) ' Smectite rate         = ', dppt(7,nx,1,1)
!!    write(*,*) ' One cell in'
!!    write(*,*) ' Smectite surface area = ', area(7,nx-1,1,1)
!!    write(*,*) ' Smectite volume       = ', volfx(7,nx-1,1,1)
!!    write(*,*) ' Smectite rate         = ', dppt(7,nx-1,1,1)
!!    write(*,*)

!  **********************  END GIMRT BLOCK  *********************************

!  Store values of master variable
  phmax = 0.0
  phchg = maxval(DABS(sp(ikmast,1:nx,1:ny,1:nz)/clg-spno2(1:nx,1:ny,1:nz)/clg))
  phloc = maxloc(DABS(sp(ikmast,1:nx,1:ny,1:nz)/clg-spno2(1:nx,1:ny,1:nz)/clg))
  jxph = phloc(1)
  jyph = phloc(2)
  jzph = phloc(3)

  spold = sp
  spexold = spex
  spsurfold = spsurf


!**********************************
  IF (iprint3 == 1) THEN
    
    IF (OS3D) THEN
      IF (MOD(nn,ScreenInterval) == 0) THEN
        IF (dtmaxcour == 0.00) THEN       !  Diffusion only problem
          WRITE(*,*) 'Time step # ',nn
          IF (OutputTimeUnits == 'years') THEN
            WRITE(*,225) time,delt
          ELSE IF (OutputTimeUnits == 'days') THEN
            WRITE(*,2251) time*OutputTimeScale,delt*OutputTimeScale
          ELSE IF (OutputTimeUnits == 'hours') THEN
            WRITE(*,2252) time*OutputTimeScale,delt*OutputTimeScale
          ELSE IF (OutputTimeUnits == 'minutes') THEN
            WRITE(*,2253) time*OutputTimeScale,delt*OutputTimeScale
          ELSE IF (OutputTimeUnits == 'seconds') THEN
            WRITE(*,2254) time*OutputTimeScale,delt*OutputTimeScale
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Time units not recognized'
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
          WRITE(*,217) iterat
          WRITE(*,227) phchg,phloc(1),phloc(2),phloc(3)
          WRITE(*,*)
        ELSE
          check =  courfactor*delt/dtmaxcour
          WRITE(*,*) 'Time step # ',nn
          IF (OutputTimeUnits == 'years') THEN
            WRITE(*,225) time,delt
          ELSE IF (OutputTimeUnits == 'days') THEN
            WRITE(*,2251) time*OutputTimeScale,delt*OutputTimeScale
          ELSE IF (OutputTimeUnits == 'hours') THEN
            WRITE(*,2252) time*OutputTimeScale,delt*OutputTimeScale
          ELSE IF (OutputTimeUnits == 'minutes') THEN
            WRITE(*,2253) time*OutputTimeScale,delt*OutputTimeScale
          ELSE IF (OutputTimeUnits == 'seconds') THEN
            WRITE(*,2254) time*OutputTimeScale,delt*OutputTimeScale
          ELSE
            WRITE(*,*)
            WRITE(*,*) ' Time units not recognized'
            WRITE(*,*)
            READ(*,*)
            STOP
          END IF
          WRITE(*,552) check
          WRITE(*,217) iterat
          WRITE(*,227) phchg,phloc(1),phloc(2),phloc(3)
          WRITE(*,*)
        END IF
      END IF
      
    ELSE
    
        
  177 format(1x,1PE12.4,1x, 1PE12.4,1x,1PE12.4,1x,1PE12.4)
              
      IF (MOD(nn,ScreenInterval) == 0) THEN
        
        WRITE(*,*) 'Time step # ',nn

        IF (SerpentineFracture) THEN

          GasFlux_FaceWest = 0.0
          GasFlux_FaceEast = 0.0
          GasFlux_FaceSouth = 0.0
          GasFlux_FaceNorth = 0.0

          jz = 1
          jx = 1
          DO jy = 1,ny
            DO kk = 1,ngas
              GasFlux_FaceWest(kk) = GasFlux_FaceWest(kk) + ag(jx,jy,1) * (spgas10(kk,jx-1,jy,1)-spgas10(kk,jx,jy,jz)) * delt
            END DO
          END DO    
      
          jz = 1
          jx = nx
          DO jy = 1,ny
            DO kk = 1,ngas
              GasFlux_FaceEast(kk) = GasFlux_FaceEast(kk) +  cg(jx,jy,1) * (spgas10(kk,jx+1,jy,1)-spgas10(kk,jx,jy,jz)) * delt
            END DO
          END DO
      
          jz = 1
          jy = 1
          DO jx = 1,nx
            DO kk = 1,ngas
              GasFlux_FaceSouth(kk) = GasFlux_FaceSouth(kk) + fg(jx,jy,1) * (spgas10(kk,jx,jy-1,1)-spgas10(kk,jx,jy,jz)) * delt
            END DO
          END DO 
      
          jz = 1
          jy = ny
          DO jx = 1,nx
            DO kk = 1,ngas
              GasFlux_FaceNorth(kk) = GasFlux_FaceNorth(kk) + dg(jx,jy,1) * (spgas10(kk,jx,jy+1,1)-spgas10(kk,jx,jy,jz)) * delt
            END DO
          END DO
      
        !!! Transport coefficient "dg" = m^3/yr so multiplying by mol/m^3 gives mol/yr. Multiplying by Delta t gives moles
      
          write(202,2255) time, GasFlux_FaceWest(1)
          write(203,2255) time, GasFlux_FaceEast(1)
          write(204,2255) time, GasFlux_FaceSouth(1)
          write(205,2255) time, GasFlux_FaceNorth(1)
          
          AqueousFlux_FaceWest = 0.0
          AqueousFlux_FaceEast = 0.0
          AqueousFlux_FaceSouth = 0.0
          AqueousFlux_FaceNorth = 0.0
      
          kk = 21
      
          jz = 1
          jx = 1
          DO jy = 1,ny
            AqueousFlux_FaceWest(kk) = AqueousFlux_FaceWest(kk) + a(jx,jy,1) * (sp10(kk,jx-1,jy,1)-sp10(kk,jx,jy,jz)) * delt
          END DO  
      
          jz = 1
          jx = nx
          DO jy = 1,ny
            AqueousFlux_FaceEast(kk) = AqueousFlux_FaceEast(kk) + c(jx,jy,1) * (sp10(kk,jx+1,jy,1)-sp10(kk,jx,jy,jz)) * delt
          END DO 
      
          jz = 1
          jy = 1
          DO jx = 1,nx
            AqueousFlux_FaceSouth(kk) = AqueousFlux_FaceSouth(kk) + f(jx,jy,1) * (sp10(kk,jx,jy-1,1)-sp10(kk,jx,jy,jz)) * delt
          END DO 
      
          jz = 1
          jy = ny
          DO jx = 1,nx
            AqueousFlux_FaceNorth(kk) = AqueousFlux_FaceNorth(kk) + d(jx,jy,1) * (sp10(kk,jx,jy+1,1)-sp10(kk,jx,jy,jz)) * delt
          END DO
      
          write(212,2255) time, AqueousFlux_FaceWest(1)
          write(213,2255) time, AqueousFlux_FaceEast(1)
          write(214,2255) time, AqueousFlux_FaceSouth(1)
          write(215,2255) time, AqueousFlux_FaceNorth(1)
          
        END IF
          
        IF (OutputTimeUnits == 'years') THEN
          WRITE(*,225) time,delt
        ELSE IF (OutputTimeUnits == 'days') THEN
          WRITE(*,2251) time*OutputTimeScale,delt*OutputTimeScale
        ELSE IF (OutputTimeUnits == 'hours') THEN
          WRITE(*,2252) time*OutputTimeScale,delt*OutputTimeScale
        ELSE IF (OutputTimeUnits == 'minutes') THEN
          WRITE(*,2253) time*OutputTimeScale,delt*OutputTimeScale
        ELSE IF (OutputTimeUnits == 'seconds') THEN
          WRITE(*,2254) time*OutputTimeScale,delt*OutputTimeScale
        ELSE
          WRITE(*,*)
          WRITE(*,*) ' Time units not recognized'
          WRITE(*,*)
          READ(*,*)
          STOP
        END IF
          
        WRITE(*,217) iterat
        WRITE(*,227) phchg,phloc(1),phloc(2),phloc(3)
        WRITE(*,*)
        
      END IF
        
    ENDIF
    iteration_tot = iterat + iteration_tot

  END IF
!**************************

!  **********************  START OS3D BLOCK  *********************************

!  **********************  END OS3D BLOCK  *********************************

  IF (RunToSteady) THEN            !  Check for steady state
    CALL SteadyState(ncomp,nx,ny,nz,delt,steady)
  END IF


!  Calculate time step to be used based on various criteria
  IF (nn > 4) THEN
        dtmax = tstep
        IF (dtmaxcour < deltmin .AND. dtmaxcour /= 0.0) THEN   !  Reset the minimum DELT if the Courant-dictated time step is even smaller
          deltmin = dtmaxcour
        END IF
        IF (dtmaxcour /= 0.0) THEN
          dtmax = MIN(dtmax,dtmaxcour)
        END IF
        dtmax = MAX(dtmax,deltmin)
        IF (ReadNuft) THEN
          dtNuft = timeNuft - time
          IF (dtNuft < eps) THEN
            dtmax = delt
          ELSE
            dtmax = MIN(dtNuft,dtmax)
          END IF
        END IF
        IF (modflow) THEN
          dtModFlow = timeModFlow - time
          IF (dtModFlow(NumModFlowSteps) < eps) THEN
            dtmax = delt
          ELSE
            dtmax = MIN(dtModFlow(NumModFlowSteps),dtmax)
          END IF
        END IF
        CALL timestep(nx,ny,nz,delt,dtold,ttol,tstep,dtmax,ikmast)
        dtold = delt
    END IF


  !IF (time+delt > prtint(nint) .AND. prtint(nint) /= time) THEN
  IF (time+delt > prtint(nint) .AND. ABS(prtint(nint) - time) > 1.0d-10) THEN ! 1.0d-14 is a small number
    delt = prtint(nint) - time
    WRITE(*,*) ' Adjusting time step to match output file'
    WRITE(*,5085) delt*OutputTimeScale
    WRITE(*,*)
  END IF

  spnno2 = spno2
  DO jz = 1,nz
    DO jy = 1,ny
      DO jx = 1,nx
        spno2(jx,jy,jz) = sp(ikmast,jx,jy,jz)
      END DO
    END DO
  END DO

!! Write breakthrough files here
  CALL BreakthroughWrite(ncomp,nspec,nkin,nrct,ngas,npot,nx,ny,nz,nseries,  &
                        nexchange,nexch_sec,nsurf,nsurf_sec,nn,ikpH,time,nplotsurface,nplotexchange )

  CALL MineralBreakthroughWrite(ncomp,nkin,nrct,nx,ny,nz,minseries,nn,time)

  CALL IsotopeBreakthroughWrite(ncomp,nspec,nkin,nrct,ngas,npot, &
    nx,ny,nz,nseries,nexchange,nexch_sec,nsurf,nsurf_sec,nn,ikpH,time,nplotsurface,nplotexchange )

  call AqueousFluxWrite(ncomp,nspec,nx,ny,nz,nn,time,delt )
  call FluxWeightedConcentrationWrite(ncomp,nspec,nx,ny,nz,nn,time,delt )

  IF (MOD(nn,ScreenInterval) == 0 .AND. MakeMovie) THEN
    write(*,*) ' ----> Writing to movie file'
    CALL movie(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,  &
          ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,delt,jpor,FirstCall)
    FirstCall = .FALSE.
  END IF

  IF (OutputTimeCounter <= NumOutputTimes .AND. OutputTime(OutputTimeCounter) /= 0.0) THEN
    IF (time+delt > OutputTime(OutputTimeCounter)) THEN
      delt = OutputTime(OutputTimeCounter) - time
      WRITE(*,*) ' Adjusting time step to match time series output'
      WRITE(*,5085) delt*OutputTimeScale
      WRITE(*,*)
    END IF
  END IF

185 FORMAT(1PE12.5,12x,100(1X,1PE16.8))

!        if (ulab(iplot(1)).eq.'h+' .or. ulab(iplot(1).eq.'H+' .and. nplot.eq.1) then
!          phwrite = -(sp(iplot(1),j)+lngamma(iplot(1),j) )/clg

  !IF (time >= prtint(nint) .OR. steady) THEN
  IF (time >= prtint(nint) .OR. steady .OR. ABS(prtint(nint) - time) <= 1.0d-10) THEN ! 1.0d-14 is a small number to avoid numerical issues in the Richards solver

    iprnt = 1
    WRITE(*,*)
    WRITE(*,*) '  WRITING OUTPUT FILES'
    IF (OutputTimeUnits == 'years') THEN
      WRITE(*,2260) time
    ELSE IF (OutputTimeUnits == 'days') THEN
      WRITE(*,2261) time*OutputTimeScale
    ELSE IF (OutputTimeUnits == 'hours') THEN
      WRITE(*,2262) time*OutputTimeScale
    ELSE IF (OutputTimeUnits == 'minutes') THEN
      WRITE(*,2263) time*OutputTimeScale
    ELSE IF (OutputTimeUnits == 'seconds') THEN
      WRITE(*,2264) time*OutputTimeScale
    ELSE
      WRITE(*,*)
      WRITE(*,*) ' Time units not recognized'
      WRITE(*,*)
      READ(*,*)
      STOP
    END IF
    WRITE(*,*) '  File number  = ', nint
    WRITE(*,*)

  END IF

  IF (FractureNetwork .and. CubicLaw) THEN
    call rmesh51(nx,ny)
  END IF

  IF (iprnt == 1) THEN
    !!Stolze Lucien: to create and overwrite the restart file at each printout
    iures = 131

    INQUIRE(FILE=RestartOutputFile,OPENED=truefalse)
    IF (truefalse) then
      CLOSE(UNIT=iures)
      OPEN(UNIT=iures,FILE=RestartOutputFile,FORM='unformatted',STATUS='unknown')
    ELSE
      OPEN(UNIT=iures,FILE=RestartOutputFile,FORM='unformatted',STATUS='unknown')
    END IF

    WRITE(iures) time
    WRITE(iures) nn
    WRITE(iures) nint
    WRITE(iures) delt,dtold,tstep,deltmin,dtmaxcour,dtmax
    WRITE(iures) keqaq
    WRITE(iures) keqgas
    WRITE(iures) keqsurf
    WRITE(iures) xgram
    WRITE(iures) spnO2
    WRITE(iures) spnnO2
    WRITE(iures) sp
    WRITE(iures) s
    WRITE(iures) sn
    WRITE(iures) sp10
    WRITE(iures) spold
    WRITE(iures) spex
    WRITE(iures) spex10
    WRITE(iures) lngamma
    WRITE(iures) exchangesites
    WRITE(iures) spexold
    WRITE(iures) spgas
    WRITE(iures) spgasold
    WRITE(iures) spgas10
    IF (isaturate==1) then
      WRITE(iures) sgas
      WRITE(iures) sgasn
    ENDIF
    if (ierode==1) then
      WRITE(iures) ssurf
      WRITE(iures) ssurfn
    endif
    WRITE(iures) sexold
    WRITE(iures) ssurfold
    WRITE(iures) spsurf
    WRITE(iures) spsurf10
    WRITE(iures) spsurfold
    WRITE(iures) raq_tot
    WRITE(iures) sion
    WRITE(iures) jinit   !! Read IntegerDummyArray(nxyz)
    WRITE(iures) keqmin
    WRITE(iures) volfx
    WRITE(iures) dppt
    WRITE(iures) area
    WRITE(iures) areainByGrid
    WRITE(iures) volinByGrid
    WRITE(iures) specificByGrid
    WRITE(iures) LogPotential
    WRITE(iures) t
    WRITE(iures) told
    WRITE(iures) ro
    WRITE(iures) por
    WRITE(iures) satliq
    WRITE(iures) qxgas
    WRITE(iures) qygas
    WRITE(iures) qzgas
    WRITE(iures) pres
    
    WRITE(iures) ActiveCell
    WRITE(iures) VolSaveByTimeStep
    WRITE(iures) Volsave
    WRITE(iures) ncounter
    
    !********************************************
    ! Edit by Toshiyuki Bandai 2024 Oct.
    IF (Richards) THEN
        WRITE(iures) Richards_State%psi
        WRITE(iures) Richards_State%theta
    END IF
    ! End of Edit by Toshiyuki Bandai 2024 Oct.
    !*********************************************

    CLOSE(UNIT=iures,STATUS='keep')

    IF (ny == 1 .AND. nz == 1) THEN

      CALL speciation(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,npot,  &
         ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,delt)
      IF (kaleidagraph) THEN
        CALL GraphicsKaleidagraph(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,  &
           ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,delt,jpor)
      END IF
      IF (xmgr) THEN
!!        CALL GraphicsXmgr(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,nsurf_sec,  &
!!           ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,delt)
          write(*,*) ' XMGR graphics option no longer supported'
          write(*,*)
          read(*,*)
          stop
      END IF
      IF (xtool) THEN
        CALL xtoolOutput(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,  &
          nsurf_sec,ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,iko2,master,delt)
      END IF
      IF (nview) THEN
        CALL NviewOutput(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,  &
          nsurf_sec,ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,ikO2,master,delt)
      END IF
      IF (tecplot) THEN

        CALL GraphicsVisit(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,  &
          ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,delt,jpor,FirstCall)
        
      END IF

    ELSE

      IF (nz == 1 .OR. ny == 1) THEN
        CALL speciation(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,npot,  &
           ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,delt)
      END IF
      IF (xtool) THEN
        CALL xtoolOutput(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,  &
          nsurf_sec,ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,iko2,master,delt)
      ELSE IF (nview) THEN
        CALL NviewOutput(ncomp,nrct,nkin,nspec,nexchange,nexch_sec,nsurf,  &
          nsurf_sec,ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,ikO2,master,delt)
      ELSE IF (visit) THEN

        CALL GraphicsVisit(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,  &
          ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,delt,jpor,FirstCall)
      ELSE

        CALL GraphicsVisit(ncomp,nrct,nkin,nspec,ngas,nexchange,nexch_sec,nsurf,nsurf_sec,  &
          ndecay,ikin,nx,ny,nz,time,nn,nint,ikmast,ikph,delt,jpor,FirstCall)

      END IF

      TotalMass = 0.0d0
      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            TotalMass = TotalMass + sp10(1,jx,jy,jz)*dxx(jx)*dyy(jy)*dzz(jx,jy,jz)
          END DO
        END DO
      END DO

      WRITE(*,*)
      WRITE(*,*)   ' Initial mass in system           = ',TotalMass
      IF (InitialTotalMass /= 0.0d0) THEN
        WRITE(*,*) ' Ratio of final to initial mass   =', TotalMass/InitialTotalMass
      END IF

      jxmax = 0
      jymax = 0
      MaxDivergence = 0.00
      DO jz = 1,nz
        DO jy = 1,ny
          DO jx = 1,nx
            Coordinate = 'X'
            call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            checkw = RoAveLeft*qx(jx-1,jy,jz)*dyy(jy)*dzz(jx,jy,jz)
            checke = RoAveRight*qx(jx,jy,jz)*dyy(jy)*dzz(jx,jy,jz)
            Coordinate = 'Y'
            call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            checkn = RoAveRight*qy(jx,jy,jz)*dxx(jx)*dzz(jx,jy,jz)
            checks = RoAveLeft*qy(jx,jy-1,jz)*dxx(jx)*dzz(jx,jy,jz)
            Coordinate = 'Z'
            call AverageRo(Coordinate,jx,jy,jz,RoAveRight,RoAveLeft)
            checkPlus = RoAveRight*qz(jx,jy,jz)*dxx(jx)*dyy(jy)
            checkMinus = RoAveLeft*qz(jx,jy,jz-1)*dxx(jx)*dyy(jy)

            !! Allocate pump
            pumpterm = 0.0d0

            !! 1) pump time series
            IF (pumptimeseries) THEN
              IF (npump(jx,jy,jz)>0) THEN
                qg(1,jx,jy,jz)=qgdum
              ELSE
             qg(1,jx,jy,jz)=0
              END IF
            pumpterm = pumpterm + qg(1,jx,jy,jz)

            !! 2) normal pump
            ELSEIF (wells) THEN

              DO npz = 1,npump(jx,jy,jz)
                pumpterm = pumpterm + qg(npz,jx,jy,jz)
              END DO

            END IF

            RealSum = ro(jx,jy,jz)* pumpterm + checkw+checks+checkMinus-checkn-checke-CheckPlus

            IF (DABS(RealSum) > MaxDivergence) THEN
                jxmax = jx
                jymax = jy
            END IF
            MaxDivergence = DMAX1(MaxDivergence,DABS(RealSum))
          END DO
        END DO
      END DO

      WRITE(*,*) ' Maximum divergence in flow field = ', MaxDivergence
      write(*,*) ' At grid cells: ',jxmax,jymax
    !!  write(*,*) qx(jxmax,jymax,1), qx(jxmax-1,jymax,1)
    !!  write(*,*) qy(jxmax,jymax,1), qy(jxmax,jymax-1,1)
    !!  read(*,*)
      WRITE(*,*)
    !!  READ(*,*)


  END IF

    IF (nint >= nstop .OR. steady) THEN
      IF (steady) THEN
        WRITE(*,*)
        WRITE(*,*) '  **** STEADY STATE ACHIEVED **** '
        WRITE(*,*)
        WRITE(*,*) ' Total Newton iterations = ',iteration_tot
        WRITE(*,*)
      ELSE
        WRITE(*,*)
        WRITE(*,*) ' Total Newton iterations = ',iteration_tot
        WRITE(*,*)
        WRITE(*,*) '  *** RUN SUCCESSFULLY COMPLETED *** '
        WRITE(*,*)
        
      END IF

      CALL date_and_time(dumm1,dumm2,dumm3,curr_time)
      end_mon = curr_time(2)
      end_day = curr_time(3)
      end_hr  = curr_time(5)
      end_min = curr_time(6)
      end_sec = curr_time(7)
      end_millisec = curr_time(8)

      MilliSeconds = DFLOAT(end_millisec)
      EndSeconds = DFLOAT(end_sec) + MilliSeconds/1000.0d0

      WRITE(*,*)
      WRITE(*,*)
      WRITE(*,*) '====== Running time on this computer ======'
      WRITE(*,*)

      IF (end_mon /= str_mon) THEN
        end_day = str_day + 1
      END IF

      call CPU_TIME(PrintSeconds)
      isimu_hr = DINT(PrintSeconds/3600.0d0)
      isimu_min = DINT(PrintSeconds/60.0d0) - isimu_hr*60.0d0
      isimu_sec = DINT(PrintSeconds) - isimu_hr - isimu_min
      ExtraSeconds = PrintSeconds - DFLOAT(isimu_hr*3600)  - DFLOAT(isimu_min*60)

!!!      WRITE(*,771) isimu_hr, isimu_min, isimu_sec
      WRITE(*,771) isimu_hr, isimu_min, ExtraSeconds
      nint = nint + 1

!   Write a restart file

    iures = 131

    INQUIRE(FILE=RestartOutputFile,OPENED=truefalse)
    IF (truefalse) then
      CLOSE(UNIT=iures)
      OPEN(UNIT=iures,FILE=RestartOutputFile,FORM='unformatted',STATUS='unknown')
    ELSE
      OPEN(UNIT=iures,FILE=RestartOutputFile,FORM='unformatted',STATUS='unknown')
    END IF

    WRITE(iures) time
    WRITE(iures) nn
    WRITE(iures) nint
    WRITE(iures) delt,dtold,tstep,deltmin,dtmaxcour,dtmax

    WRITE(iures) keqaq
    WRITE(iures) keqgas
    WRITE(iures) keqsurf
    WRITE(iures) xgram
    WRITE(iures) spnO2
    WRITE(iures) spnnO2
    WRITE(iures) sp
    WRITE(iures) s
    WRITE(iures) sn
    WRITE(iures) sp10
    WRITE(iures) spold
    WRITE(iures) spex
    WRITE(iures) spex10
    WRITE(iures) lngamma
    WRITE(iures) exchangesites
    WRITE(iures) spexold
!!    WRITE(iures) sch
    WRITE(iures) spgas
    WRITE(iures) spgasold
    WRITE(iures) spgas10
    IF (isaturate==1) then
      WRITE(iures) sgas
      WRITE(iures) sgasn
    ENDIF
    if (ierode==1) then
      WRITE(iures) ssurf
      WRITE(iures) ssurfn
    endif
    WRITE(iures) sexold
    WRITE(iures) ssurfold
    WRITE(iures) spsurf
    WRITE(iures) spsurf10
    WRITE(iures) spsurfold
    WRITE(iures) raq_tot
    WRITE(iures) sion
    WRITE(iures) jinit

!!    WRITE(iures) mumin_decay
    WRITE(iures) keqmin
    WRITE(iures) volfx
    WRITE(iures) dppt
    WRITE(iures) area
    WRITE(iures) areainByGrid
    WRITE(iures) volinByGrid
    WRITE(iures) specificByGrid
    WRITE(iures) LogPotential

    WRITE(iures) t
    WRITE(iures) told
    WRITE(iures) ro

    WRITE(iures) por

    WRITE(iures) satliq
    WRITE(iures) qxgas
    WRITE(iures) qygas
    WRITE(iures) qzgas
    WRITE(iures) pres
  !!!  WRITE(iures) dspy
  !!!  WRITE(iures) dspz
  !!!  WRITE(iures) qg
    WRITE(iures) ActiveCell
    WRITE(iures) VolSaveByTimeStep
    WRITE(iures) Volsave
    WRITE(iures) ncounter
    
    !********************************************
    ! Edit by Toshiyuki Bandai 2023 May
    IF (Richards) THEN
        WRITE(iures) Richards_State%psi
        WRITE(iures) Richards_State%head
        WRITE(iures) Richards_State%theta
        WRITE(iures) Richards_State%theta_prev
    END IF
    ! End of Edit by Toshiyuki Bandai 2023 May
    !*********************************************

!!    WRITE(iures) tauZero

!    CLOSE(UNIT=intfile,STATUS='keep')

    CLOSE(UNIT=iures,STATUS='keep')

!         ***** PETSc closeout**************
         IF (InputFileCounter == NumInputFiles) THEN
           IF (petscon) THEN
             call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
           END IF
           IF (CalculateFlow) THEN
             call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
           END IF
           IF (OS3Dpetsc) THEN
             call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
           END IF
           call PetscFinalize(ierr)
         ELSE
           IF (petscon) THEN
             call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
           END IF
           IF (CalculateFlow) THEN
             call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
           END IF
           IF (OS3Dpetsc) THEN
             call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
           END IF
         END IF

!       ****** PETSc closeout finished *********

       IF (InputFileCounter < NumInputFiles ) THEN
         InputFileCounter = InputFileCounter + 1
         dumstring = InputFile(InputFileCounter)
         CALL stringlen(dumstring,ls)
         WRITE(*,*)
         WRITE(*,*) ' Restarting CRUNCH with input file: ', dumstring(1:ls)
         WRITE(*,*)
         NewInput = .TRUE.
       ELSE
         NewInput = .FALSE.
       END IF

       CLOSE(iunit2,status='keep')

       IF (nplot > 0) THEN
         DO ll = 1,nseries
           intfile = 100+ll
           CLOSE(intfile)
         END DO
       END IF
       
!!       READ(*,*)
       RETURN
!!    ********  NORMAL STOP HERE  **************
    END IF
    
    nint = nint + 1
    IF (time+delt > prtint(nint) .AND. ABS(prtint(nint) - time) > 1.0d-14) THEN
      delt = prtint(nint) - time
      WRITE(*,*) ' Adjusting time step to match output file'
      WRITE(*,5085) delt*OutputTimeScale
      WRITE(*,*)
    END IF
    iprnt = 0
  END IF

END DO    ! End of time loop

!  ***************  END OF TIME LOOP  **********************

201 WRITE(*,*) '   Maximum time steps completed--contact Carl Steefel at CISteefel@lbl.gov'
    WRITE(*,*)
    GOTO 202

555 WRITE(*,*) ' Error opening saturation file'
556 WRITE(*,*) ' Error in reading breakthrough file'
!         ***** PETSc closeout**************
           IF (petscon) THEN
             call CrunchPETScFinalizeSolver(xvec,bvec,amatpetsc,userC,ierr)
           END IF
           IF (CalculateFlow) THEN
             call CrunchPETScFinalizeSolver(xvecP,bvecP,amatP,userP,ierr)
           END IF
           IF (OS3Dpetsc) THEN
             call CrunchPETScFinalizeSolver(xvecD,bvecD,amatD,userD,ierr)
           END IF
           call PetscFinalize(ierr)
!       ****** PETSc closeout finished *********
202 READ(*,*)
STOP


!********************************
771 FORMAT(5X, 'hr:', i4, '    min:', i4, '    sec: ', f7.2 ,/)
1009 FORMAT(10(1X,1PE12.4))
211 FORMAT(2X,f10.6,2X,f9.3,2X,f10.6)
217 FORMAT(2X,'Number of Newton iterations = ',i2)
218 FORMAT(2X,'# of SOR iterations = ',i3)
224 FORMAT(2X,a18,2X,'Max residual = ',1PE12.4,2X,'Grid pt =',i3)
225 FORMAT(2X,'Time (yrs) = ',1PE12.5,2X,'Delt (yrs) =',1PE10.3)
2251 FORMAT(2X,'Time (days) = ',1PE12.5,2X,'Delt (days) =',1PE10.3)
2252 FORMAT(2X,'Time (hrs) = ',1PE12.5,2X,'Delt (hrs) =',1PE10.3)
2253 FORMAT(2X,'Time (mins) = ',1PE12.5,2X,'Delt (mins) =',1PE10.3)
2254 FORMAT(2X,'Time (secs) = ',1PE12.5,2X,'Delt (secs) =',1PE10.3)
2255 FORMAT( 2X,1PE12.5,2X,4(1X,1PE10.3) )
2260 FORMAT(2X,'Time (yrs) = ',1PE10.3)
2261 FORMAT(2X,'Time (days) = ',1PE10.3)
2262 FORMAT(2X,'Time (hrs) = ',1PE10.3)
2263 FORMAT(2X,'Time (mins) = ',1PE10.3)
2264 FORMAT(2X,'Time (secs) = ',1PE10.3)
226 FORMAT(2X,'ERST = ',1PE10.3)
227 FORMAT(2X,'Maximum change in master variable = ',1PE10.2,2X,' at grid pts ',i4,1x,i4,1x,i4)
228 FORMAT(2X,'Maximum number of Newton iterations exceeded')
229 FORMAT(2X,'New number of grid points ',i3)
230 FORMAT(2X,'Maximum change in master variable   =   ',1PE10.2)
231 FORMAT(2X,'Stationary state convergence criteria = ',1PE10.2)
2005 FORMAT(2X,1PE12.4,2X,10(1PE12.4))
505 FORMAT(2X,'Initial # of pts. in coarse grid =',i4)
506 FORMAT(2X,'# of pts. in fine grid  =',i5)
507 FORMAT(2X,'Initial grid spacing (m) =',1PE12.4)
2021 FORMAT(2X,'------ Segment',i2,'------')
508 FORMAT(2X,'Fine (fixed) grid spacing (m) =',1PE12.4)
801 FORMAT(2X,'Master variable: ',a18)
802 FORMAT(2X,1PE12.4,2X,1PE12.3)
551 FORMAT(2X,'Maximum grid Peclet number = ',1PE13.3)
552 FORMAT(2X,'Maximum Courant number = ',1PE13.3)
855 FORMAT(1X,10(1PE12.4))
5022 FORMAT(1X,'Initial porosity (%) = ',f8.4)
5023 FORMAT(1X,'Time for basis switching ',1PE12.2)
5024 FORMAT(1X,'Maximum conc. correction ',1PE12.4)
2003 FORMAT(2X,'Diffusion coeff (m**2/sec) =  ',1PE11.3)
2004 FORMAT(2X,'Formation factor =            ',f11.2)
2031 FORMAT(2X,'Porosity dependence in D* =   ',f11.2)
997 FORMAT(1X,f12.4,4(1X,1PE12.4))
901 FORMAT(1X,'# reaction steps = ',i4,1X,i3,1X, 1PE12.4,1X,1PE12.4)
501 FORMAT(13X,a80)
325 FORMAT(2X,'Maximum error = ',e10.3)
326 FORMAT(2X,'Errmax = ',e10.2)
100  FORMAT(1X,2(1PE12.4))
705 FORMAT(1X,1PE12.5,150(1X,1PE22.14))
707 FORMAT(1X,1PE12.5,11x,150(1X,1PE22.14))
706 FORMAT(1X,'#   Time      ',100(1X,a12))
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
5031 FORMAT(30(1X,1PE10.2))
5032 FORMAT(30(1X,i3))
5085 FORMAT(1X,' ---> New time step = ',1PE11.4)
5086 FORMAT(1X,' ---> Old time step = ',1PE11.4)

!         ***** PETSc closeout**************
       if( petscon) then
        call VecDestroy(bvec,ierr)
        call VecDestroy(xvec,ierr)
        call MatDestroy(amatpetsc,ierr)
        call KSPDestroy(ksp,ierr)
        call PetscFinalize(ierr)
       endif
       if( OS3Dpetsc) then
        call VecDestroy(bvec,ierr)
        call VecDestroy(xvec,ierr)
        call MatDestroy(amatD,ierr)
        call KSPDestroy(ksp,ierr)
        call PetscFinalize(ierr)
       endif
!       ****** PETSc closeout finished *********
READ(*,*)
STOP
END SUBROUTINE CrunchTope
!****************END DRIVER PROGRAM**********************
