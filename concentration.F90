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
    
MODULE  concentration

  USE crunchtype
  USE params

  LOGICAL(LGT)                                      :: O2Found  
  LOGICAL(LGT)                                      :: IncludeBdot
  
  LOGICAL(LGT)                                      :: ReadInitialConditions
  
  CHARACTER (LEN=mls)                               :: InitialConditionsFile
  

  INTEGER(I4B)                                      :: ikh2o  !  integer pointer for species number for H2O
  INTEGER(I4B)                                      :: nplot  !  number of time series plot files
  INTEGER(I4B)                                      :: interval  !  interval to write to screen
  INTEGER(I4B)                                      :: jxplot
  INTEGER(I4B)                                      :: jyplot
  INTEGER(I4B)                                      :: jzplot
  INTEGER(I4B)                                      :: ilabel
  INTEGER(I4B)                                      :: iexc
  INTEGER(I4B)                                      :: iplotall
  INTEGER(I4B)                                      :: iplotpH
  INTEGER(I4B)                                      :: nreactkinmax
  INTEGER(I4B)                                      :: nmonodmax
  INTEGER(I4B)                                      :: nmonodaqmax
  INTEGER(I4B)                                      :: ninhibitmax
  INTEGER(I4B)                                      :: ninhibitaqmax
  INTEGER(I4B)                                      :: npointO2gas
  INTEGER(I4B)                                      :: npointH2gas
  INTEGER(I4B)                                      :: npointO2aq
  INTEGER(I4B)                                      :: NPestExchange
  INTEGER(I4B)                                      :: UnitPestExchange
  INTEGER(I4B)                                      :: NPestSurface
  INTEGER(I4B)                                      :: UnitPestSurface

  INTEGER(I4B)                                      :: ikFe2
  INTEGER(I4B)                                      :: ikFe3
  INTEGER(I4B)                                      :: ikNa
  INTEGER(I4B)                                      :: ikCl

!  Temporary arrays used in the initialization (can be deallocated)

  REAL(DP), DIMENSION(:), ALLOCATABLE               :: keqaq_tmp
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: keqgas_tmp
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: keqsurf_tmp
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: sptmp  !  temporary array for ln species concentration
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: sptmp10
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: gamtmp
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: stmp
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: dxxt
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: dyyt
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: dzzt
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: eqgas
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: eqhom
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: sexch
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: spextmp
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: spextmp10
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: totextmp
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: spgastmp
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: spgastmp10
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: sgastmp
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: spsurftmp
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: spsurftmp10
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: ssurftmp
  REAL(DP), DIMENSION(:,:), ALLOCATABLE             :: dpsi
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: SurfaceCon
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: ExchangeCon 
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: TotChargeSave
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: PestScale

  REAL(DP), DIMENSION(:), ALLOCATABLE               :: DatabaseTemperature

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE           :: nbasin
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE           :: nbkin
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE           :: nvx
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE           :: nvy
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE           :: nvz
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE           :: SolidDensityFrom
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: SolidSolutionRatio
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: SolidDensity

  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE    :: stringarray
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE    :: PestExchangeList
  CHARACTER (LEN=mls)                               :: PestExchangeUnits
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE    :: PestSurfaceList
  CHARACTER (LEN=mls)                               :: PestSurfaceUnits

  LOGICAL(LGT), DIMENSION(:,:), ALLOCATABLE         :: equilibrate

!  Permanent arrays, but can be reallocated to smaller size

  REAL(DP), DIMENSION(5)                           :: adhcoeff 
  REAL(DP), DIMENSION(5)                           :: bdhcoeff 
  REAL(DP), DIMENSION(5)                           :: bdtcoeff 
  REAL(DP), DIMENSION(6)                           :: xgrambnd
  INTEGER(I4B), DIMENSION(6)                       :: jc 
  INTEGER(I4B), DIMENSION(2)                       :: MeanSalt

  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: ulab
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: namg
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: namkin
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: namexc
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: nam_exchsec
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: namsurf
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: namsurf_sec

  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: fweight
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: fanalyt
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: as1
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: as2
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spb
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spbgas
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: wtcomp
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: acmp
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: chg
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: adh
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: bdh
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: bdot
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: bdotParameter
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: alogkp
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: wtaq
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: zsurf
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: satkin
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spexb
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spsurfb
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: eqsurf
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: muaq
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: mugas
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: musurf
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: mukin
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: keqkin  !  Equilbrium constants for aqueous kinetic reactions
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: keqexc
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: bfit
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: muexc
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: ratek
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dependk
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: pre_raq
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: rdkin
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: raq
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: totex
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: sumactivity
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: tec
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: wt_aexch
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: aexch
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: gam_local
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: sNCexch_local
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: sNCsurf_local
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: ssurf_local
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: scond
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spcond
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spcond10
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spcondgas
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spcondgas10
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spcondex
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spcondex10
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spcondsurf
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: spcondsurf10
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: LogPotentialInit
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: ctot
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: guess
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: c_surf
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: guess_surf
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: gaspp
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: totexch
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: cec
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: halfsataq
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: termMonod
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: rinhibitaq
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: distrib
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: sbnd
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: s_local
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: half_life
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: temp_ratio  
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: ratio_isotope_init  
  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE    :: ratio_isotope
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: decay_correct
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: mole_rad
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: conversion
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: OneOverMassFraction
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: sumjackin

  REAL(DP), DIMENSION(:), ALLOCATABLE                :: sce
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: scw
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: scn
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: scs
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: sex_east
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: sex_west
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: sex_north
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: sex_south
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: schg
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: sdsp
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: sge
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: sgw
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: sgn
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: sgs
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: surf_east
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: surf_west
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: surf_north
  REAL(DP), DIMENSION(:), ALLOCATABLE                :: surf_south
  REAL(DP), DIMENSION(:),ALLOCATABLE                 :: sgaspump
  REAL(DP), DIMENSION(:),ALLOCATABLE                 :: sppTMP
  REAL(DP), DIMENSION(:),ALLOCATABLE                 :: sppTMP10
  REAL(DP), DIMENSION(:),ALLOCATABLE                 :: sppTMPperturb
  REAL(DP), DIMENSION(:),ALLOCATABLE                 :: sppTMP10perturb

  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: AqueousFluxSeriesFile
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: AqueousFluxSeriesSpecies

  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: FluxWeightedConcentrationFile
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE              :: FluxWeightedConcentrationSpecies

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: iplotAqueousFlux
!!  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jxAqueousFluxSeries
!!  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jyAqueousFluxSeries
!!  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jzAqueousFluxSeries
  INTEGER(I4B)                                     :: nAqueousFluxSeriesFile 
  INTEGER(I4B)                                     :: nplotAqueousFlux
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jxAqueousFluxSeries_lo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jxAqueousFluxSeries_hi
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jyAqueousFluxSeries_lo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jyAqueousFluxSeries_hi
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jzAqueousFluxSeries_lo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jzAqueousFluxSeries_hi

  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: CumulativeXflux
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: InstantaneousXflux
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: XfluxWeightedConcentration
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: YfluxWeightedConcentration
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: ZfluxWeightedConcentration

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: iplotFluxWeightedConcentration
  INTEGER(I4B)                                     :: nFluxWeightedConcentrationFile 
  INTEGER(I4B)                                     :: nplotFluxWeightedConcentration
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jxFluxWeightedConcentration_lo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jxFluxWeightedConcentration_hi
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jyFluxWeightedConcentration_lo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jyFluxWeightedConcentration_hi
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jzFluxWeightedConcentration_lo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jzFluxWeightedConcentration_hi

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: ksurf
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: iplot
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: nreactkin
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jxseries
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jyseries
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jzseries
  
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jxminseries
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jyminseries
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jzminseries
  
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: iexchange
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: kexch
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: ixlink
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: nclink
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: islink
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: nptlink
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: itot
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: itype
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: ndist
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jxxlo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jxxhi
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jyylo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jyyhi
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jzzlo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jzzhi
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: jjfix
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: iaqtype
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: nmonodaq
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: ninhibitaq
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: imonodaq
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: inhibitaq
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: itot_monodaq
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: itot_inhibitaq
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: idecay
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: nisotope 
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: kdecay
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: nmindecay
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: nrad_decay
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: kradpoint
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: npradpoint
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: iraddecay
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: icec

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: GasPressureTotal
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: GasPressureTotalInit

  CHARACTER (LEN=mls), DIMENSION(:,:), ALLOCATABLE :: decay_label
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: namtope
  CHARACTER (LEN=mls), DIMENSION(:,:), ALLOCATABLE :: ncon
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: condlabel
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: condtitle


! Allocatable arrays that are local (not dimensioned over spatial domain)

!  Allocatable arrays dimensioned over spatial domain

  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: keqaq
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: keqgas
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: keqsurf
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: xgram
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: xgramOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: H2Oreacted
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: spnO2
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: spnnO2
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: sp
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: s
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: sn
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: sp10
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: spold
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: spex
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: spex10
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: gam
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: exchangesites
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: spexold
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: sch
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: spgas
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: spgasold
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: spgas10
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: sgas
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: sgasn
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: ssurf
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: ssurfn
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: LogTotalSurface
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: sexold
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: ssurfold
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: skdold
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: spsurf
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: spsurf10
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: spsurfold 
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: raq_tot
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: sion

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: vrSave
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: vrInitial
  


  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: jinit

!!  Xtool stuff from Olivier

  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE          :: exflag
  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE          :: surflag


END MODULE concentration

