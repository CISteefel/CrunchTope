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
    
MODULE mineral

  USE crunchtype
  USE params

  LOGICAL(LGT)                                     :: JacobianNumerical
  INTEGER(I4B), PARAMETER                          :: mlen=mls
  INTEGER(I4B)                                     :: nreactmax
  INTEGER(I4B)                                     :: nradmax

  INTEGER(I4B)                                     :: ikCa
  INTEGER(I4B)                                     :: ikCO3
  INTEGER(I4B)                                     :: ik234U
  INTEGER(I4B)                                     :: ik238U
  INTEGER(I4B)                                     :: ikAl
  INTEGER(I4B)                                     :: kUCalcite
  INTEGER(I4B)                                     :: kMarineCalcite
  INTEGER(I4B)                                     :: kUPlag

  REAL(DP)                                         :: vdissmax
  REAL(DP)                                         :: vpptmax
  REAL(DP)                                         :: DistributionCalcite
  
  REAL(DP)                                         :: IntervalBelowEquilibrium
  REAL(DP)                                         :: IntervalAboveEquilibrium

! Temporary arrays used in the initialization (can be deallocated)

  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: keqmin_tmp
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: surfcharge_init
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: LogPotential_tmp

  CHARACTER (LEN=mlen), DIMENSION(:,:,:), ALLOCATABLE :: namdep_nyf

! Permanent arrays, but can be reallocated to smaller size

  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: umin
  CHARACTER (LEN=mls), DIMENSION(:,:), ALLOCATABLE :: rlabel
  CHARACTER (LEN=mls), DIMENSION(:,:), ALLOCATABLE :: crossaff
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE   :: namAssociate

  REAL(DP), DIMENSION(:,:), ALLOCATABLE              :: surf
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: volmol
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: wtmin
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: alnk
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: sumrd
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: volb
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: areab
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dependex
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dependsurf
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: mumin
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: jac_rmin
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: areain
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: volin
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: MineralMoles
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: si
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: sat1
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: sat2
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: Ea
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: rate0
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: silog
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: siln
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: snorm
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: ssa
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: actenergy
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: pre_rmin
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: rmin
  
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: Azero25C
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: Bnucleation
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: sigmaNucleation
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: SurfaceAreaNucleation
  

  LOGICAL(LGT), DIMENSION(:,:), ALLOCATABLE        :: HomogeneousNucleation
  LOGICAL(LGT), DIMENSION(:,:), ALLOCATABLE        :: SumMineralSurfaceArea
  
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: NucleationSurface
  
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: kNucleationPath
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: npNucleationPath
  
  
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: depend
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: halfsat
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: rinhibit
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: site_density
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: specific
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: surfcharge
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: AffinityDepend1
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: AffinityDepend2  
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: AffinityDepend3
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: voltemp

! biomass
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: BQ_min
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: BQ_kin
! biomass end

  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE          :: LocalEquilibrium
  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE          :: MineralAssociate

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: MineralID
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: ivolume
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: lenmin
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: ndependex
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: ndependsurf
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: ixdepend
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: isdepend
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: nreactmin
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: imintype
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: ndepend
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: nmonod
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: ninhibit
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: kcrossaff
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: idepend
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: imonod
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: kmonod
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: inhibit
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: itot_min
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: itot_monod
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: itot_inhibit
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: iarea
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: iedl
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: ispot

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: kpot
  LOGICAL(LGT), DIMENSION(:), ALLOCATABLE          :: kPotential
  
!! Hyperbolic Inhibition Term
!! Kformation(np,nkin),HyperbolicInhibitionName(np,nkin),HyperbolicInhibitionDepend(np,nkin)
  LOGICAL(LGT), DIMENSION(:,:), ALLOCATABLE          :: HyperbolicInhibition
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE          :: HyperbolicInhibitionPointer
  REAL(DP), DIMENSION(:,:), ALLOCATABLE              :: Kformation
  REAL(DP), DIMENSION(:,:), ALLOCATABLE              :: HyperbolicInhibitionDepend
  CHARACTER (LEN=mlen), DIMENSION(:,:), ALLOCATABLE  :: HyperbolicInhibitionName
  
! biomass
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: chi_min
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: direction_min
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: chi_kin
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: direction_kin
! biomass end

!  Allocatable arrays dimensioned over spatial domain

  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE    :: mumin_decay
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: keqmin
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: silogGlobal
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: volSaveByTimeStep
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: volSave
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: volfx
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: dppt
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: area
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: LogPotential
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: VolumeLastTimeStep
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: rminSaveForDePaolo
  
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: specificByGrid
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: areainByGrid
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: volinByGrid

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: muUranium234
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: muUranium238
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: muCalcium


  REAL(DP), DIMENSION(:), ALLOCATABLE              :: muUranium234Boundary
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: muUranium238Boundary
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: muCalciumBoundary

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: muUranium234Bulk
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: muUranium238Bulk
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: muCalciumBulk


!  Allocatable arrays dimensioned locally (not over spatial domain)

! biomass
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: muminTMP
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: keqminTMP
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: mukinTMP
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: keqkinTMP
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: LagTimeMineral
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: RampTimeMineral
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: ThresholdConcentrationMineral
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: LagTimeAqueous
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: RampTimeAqueous
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: ThresholdConcentrationAqueous

  integer(i4b),dimension(:,:),allocatable          :: SubstrateForLagMineral
  integer(i4b),dimension(:),allocatable            :: SubstrateForLagAqueous

  logical(lgt),dimension(:,:),allocatable          :: UseMetabolicLagMineral
  logical(lgt),dimension(:),allocatable            :: UseMetabolicLagAqueous

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE        :: tauZeroMineral
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE          :: tauZeroAqueous
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE        :: MetabolicLagMineral
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE          :: MetabolicLagAqueous
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE          :: SatLog

  INTEGER(I4B)                                     :: nMonodBiomassMineral
  INTEGER(I4B)                                     :: nMonodBiomassAqueous

  integer(i4b),dimension(:),allocatable            :: ibiomass_kin
  integer(i4b),dimension(:),allocatable            :: ibiomass_min
  integer(i4b),dimension(:),allocatable            :: p_cat_min
  integer(i4b),dimension(:),allocatable            :: p_cat_kin
  integer(i4b),dimension(:,:),allocatable          :: biomass_decay ! added for decay
  integer(i4b),dimension(:),allocatable            :: mintype       ! mintype=0 mineral, mintype=1 biomass
! biomass end

  REAL(DP), dimension(4,31)                        :: rate0Prime
  
  REAL(DP),DIMENSION(:,:,:), ALLOCATABLE    :: crankLogK
  
  LOGICAL(LGT)                              :: ContactPressureLogical
  LOGICAL(LGT)                              :: nmmLogical
  LOGICAL(LGT)                              :: SaltCreep
  LOGICAL(LGT)                              :: CalciteCreep
  LOGICAL(LGT)                              :: SerpentineFracture
  LOGICAL(LGT)                              :: CriticalZone
  
  CHARACTER (LEN=mlen)                      :: AqueousKineticFile
  
  CHARACTER (LEN=mlen)                      :: CatabolicKineticFile

END MODULE mineral
