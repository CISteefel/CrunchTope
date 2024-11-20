MODULE TYPES_AND_PARAMETERS

  IMPLICIT NONE

  ! ****************************** New KINDs ***********************************

  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: SIP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: LGT = KIND(.TRUE.)
  
  ! ****** Maximum string length to read an input keyword/parameter ************

  INTEGER(I4B), PARAMETER              :: max_len_str=500

  ! **************** Mathematic and physic constants ***************************

  REAL(DP), PARAMETER            :: rgasJ      = 8.3144621D0       !!  J/deg-mol
  REAL(DP), PARAMETER            :: rgas       = rgasJ / 1000D0    !!  kJ/deg-mol
  REAL(DP), PARAMETER            :: rgasKCAL   = rgas / 4.1855D0   !!  kcal/deg-mol
  REAL(DP), PARAMETER            :: tk25       = 1D0 / 298.15D0 
  REAL(DP), PARAMETER            :: clg        = LOG(10D0)
  REAL(DP), PARAMETER            :: secyr      = 365D0 * 24D0 * 60D0 * 60D0
  REAL(DP), PARAMETER            :: pi         = DACOS(-1.0D0)
  REAL(DP), PARAMETER            :: faraday    = 96485.3321233D0
  REAL(DP), PARAMETER            :: SQRT_debye = SQRT(rgasJ / (2D0 * faraday * faraday * 1000D0 )) !  Debye length=SQRT_debye*DSQRT(epsilonwater*T/I)

  ! ********************* Numerical constants ********************************

  REAL(DP), PARAMETER            :: atol   = 1D-08
  REAL(DP), PARAMETER            :: rtol   = 1D-06
  INTEGER(I4B), PARAMETER        :: nitmax = 100
  INTEGER(I4B), PARAMETER        :: north  = 1
  INTEGER(I4B), PARAMETER        :: nend   = HUGE(nitmax)
  INTEGER(I4B), PARAMETER        :: newton = 100
  INTEGER(I4B), PARAMETER        :: ione   = 1
  REAL(DP), PARAMETER            :: eps    = 1D-12 

  !************** Definition of the TYPEs ************************************


  TYPE TORTUOSITY_INFO
  
      INTEGER(I4B)                :: jxx_lo
      INTEGER(I4B)                :: jyy_lo
      INTEGER(I4B)                :: jzz_lo     
      INTEGER(I4B)                :: jxx_hi
      INTEGER(I4B)                :: jyy_hi
      INTEGER(I4B)                :: jzz_hi  
      REAL(DP)                    :: value

  END TYPE TORTUOSITY_INFO


  TYPE HETEROGENEITY_ARRAY
    ! nxyz, DEALLOCATED at the end of startclay
    INTEGER(I4B)                  :: jxxlo
    INTEGER(I4B)                  :: jxxhi
    INTEGER(I4B)                  :: jyylo
    INTEGER(I4B)                  :: jyyhi
    INTEGER(I4B)                  :: jzzlo
    INTEGER(I4B)                  :: jzzhi
    INTEGER(I4B)                  :: jjfix
    INTEGER(I4B)                  :: condition
    INTEGER(I4B)                  :: TypeOfBC

  END TYPE HETEROGENEITY_ARRAY


  TYPE OUTPUT_INFO_ARRAY  
      !scalar%nseries
      CHARACTER (LEN=max_len_str)      :: TimeSeriesFile
      ! nseries
      INTEGER(I4B)                     :: jxseries
      INTEGER(I4B)                     :: jyseries
      INTEGER(I4B)                     :: jzseries
      INTEGER(I4B)                     :: process_in_charge
      ! NumOUtputTimes+1
      REAL(DP)                         :: OutputTime
      ! scalar%nstop
      REAL(DP)                         :: prtint 
  END TYPE OUTPUT_INFO_ARRAY


TYPE SCALAR_V
   CHARACTER (LEN=2*max_len_str)  ::    ltitle
   CHARACTER (LEN=max_len_str)    ::    namefileout
   CHARACTER (LEN=max_len_str)    ::    data1
   CHARACTER (LEN=max_len_str)    ::    DensityModule
   CHARACTER (LEN=max_len_str)    ::    master
   CHARACTER (LEN=max_len_str)    ::    OutputTimeUnits
   CHARACTER (LEN=max_len_str)    ::    PermFileFormat
   CHARACTER (LEN=max_len_str)    ::    restartfile
   CHARACTER (LEN=max_len_str)    ::    RestartOutputFile
   CHARACTER (LEN=max_len_str)    ::    SaturationFile
   CHARACTER (LEN=max_len_str)    ::    SaturationFileFormat
   CHARACTER (LEN=max_len_str)    ::    TemperatureFileFormat
   CHARACTER (LEN=max_len_str)    ::    Tfile
   CHARACTER (LEN=max_len_str)    ::    TortuosityFileFormat
   CHARACTER (LEN=max_len_str)    ::    TortuosityOption
   CHARACTER (LEN=max_len_str)    ::    VelocityFileFormat
   INTEGER(I4B)    ::    icomplete
   INTEGER(I4B)    ::    icvg
   INTEGER(I4B)    ::    idiffus
   INTEGER(I4B)    ::    igamma
   INTEGER(I4B)    ::    ikh2o
   INTEGER(I4B)    ::    ikin
   INTEGER(I4B)    ::    ikmast
   INTEGER(I4B)    ::    ikO2
   INTEGER(I4B)    ::    ikH2
   INTEGER(I4B)    ::    ikph
   INTEGER(I4B)    ::    ikTracer
   INTEGER(I4B)    ::    interval
   INTEGER(I4B)    ::    intfile
   INTEGER(I4B)    ::    iprnt
   INTEGER(I4B)    ::    irestart
   INTEGER(I4B)    ::    isaturate
   INTEGER(I4B)    ::    ispeciate
   INTEGER(I4B)    ::    iterat
   INTEGER(I4B)    ::    iterat_previous
   INTEGER(I4B)    ::    itsiterate
   INTEGER(I4B)    ::    iunit1
   INTEGER(I4B)    ::    iunit2
   INTEGER(I4B)    ::    iunit3
   INTEGER(I4B)    ::    iunit5
   INTEGER(I4B)    ::    iures
   INTEGER(I4B)    ::    maxitsksp
   INTEGER(I4B)    ::    nb_electrodes
   INTEGER(I4B)    ::    nb_iterations_without_a_problem
   INTEGER(I4B)    ::    nb_potentials
   INTEGER(I4B)    ::    nchem
   INTEGER(I4B)    ::    ncomp
   INTEGER(I4B)    ::    ndl
   INTEGER(I4B)    ::    nexch_sec
   INTEGER(I4B)    ::    nexchange
   INTEGER(I4B)    ::    ngas
   INTEGER(I4B)    ::    nint
   INTEGER(I4B)    ::    nkin
   INTEGER(I4B)    ::    nkd
   INTEGER(I4B)    ::    nplot
   INTEGER(I4B)    ::    npointH2gas
   INTEGER(I4B)    ::    npointO2aq
   INTEGER(I4B)    ::    npointO2gas
   INTEGER(I4B)    ::    nprimary
   INTEGER(I4B)    ::    StartExch     
   INTEGER(I4B)    ::    StartSurf     
   INTEGER(I4B)    ::    StartPot      
   INTEGER(I4B)    ::    StartDL       
   INTEGER(I4B)    ::    posI          
   INTEGER(I4B)    ::    posGammaWater 
   INTEGER(I4B)    ::    posSelfPot    
   INTEGER(I4B)    ::    posCurrE      
   INTEGER(I4B)    ::    posCurrN      
   INTEGER(I4B)    ::    posCurrU      
   INTEGER(I4B)    ::    posLNval      
   INTEGER(I4B)    ::    posval
   INTEGER(I4B)    ::    nseries
   INTEGER(I4B)    ::    nspec
   INTEGER(I4B)    ::    nstop
   INTEGER(I4B)    ::    nsurf
   INTEGER(I4B)    ::    nsurf_sec
   INTEGER(I4B)    ::    ntemp
   INTEGER(I4B)    ::    nTortuosityMPZone
   INTEGER(I4B)    ::    nTortuosityZone
   INTEGER(I4B)    ::    number_of_surfaces
   INTEGER(I4B)    ::    NumOutputTimes
   INTEGER(I4B)    ::    numprocs
   INTEGER(I4B)    ::    nx
   INTEGER(I4B)    ::    nx_bc
   INTEGER(I4B)    ::    nxyz
   INTEGER(I4B)    ::    nxyz_bc
   INTEGER(I4B)    ::    ny
   INTEGER(I4B)    ::    ny_bc
   INTEGER(I4B)    ::    nz
   INTEGER(I4B)    ::    nz_bc
   INTEGER(I4B)    ::    O2H2production
   INTEGER(I4B)    ::    OMPthreads
   INTEGER(I4B)    ::    OutputTimeCounter
   INTEGER(I4B)    ::    OutputFileNumber
   INTEGER(I4B)    ::    ScreenInterval
   INTEGER(I4B)    ::    total_number_of_kinetics
   INTEGER(I4B)    ::    zero_electrode
   INTEGER(I4B)    ::    one_electrode
   INTEGER(I4B)    ::    more_electrodes
   INTEGER(I4B)    ::    oneD  
   INTEGER(I4B)    ::    twoD  
   INTEGER(I4B)    ::    threeD
   
   LOGICAL(LGT)    ::    Activity_gradient
   LOGICAL(LGT)    ::    AppendRestart
   LOGICAL(LGT)    ::    charge_balance_with_H
   LOGICAL(LGT)    ::    Constant_Tortuosity
   LOGICAL(LGT)    ::    Cylindrical
   LOGICAL(LGT)    ::    Fracture2D
   LOGICAL(LGT)    ::    IncludeBdot
   LOGICAL(LGT)    ::    mont_terri
   LOGICAL(LGT)    ::    O2Found
   LOGICAL(LGT)    ::    plot_porosities
   LOGICAL(LGT)    ::    plot_ionic_strength
   LOGICAL(LGT)    ::    plot_charge_balance
   LOGICAL(LGT)    ::    plot_gamma_water
   LOGICAL(LGT)    ::    plot_gas
   LOGICAL(LGT)    ::    plot_SIP
   LOGICAL(LGT)    ::    plot_DL_potentials
   LOGICAL(LGT)    ::    plot_temperature
   LOGICAL(LGT)    ::    plot_tot_aqueous
   LOGICAL(LGT)    ::    plot_conc_aqueous
   LOGICAL(LGT)    ::    plot_act_aqueous
   LOGICAL(LGT)    ::    plot_tot_DL
   LOGICAL(LGT)    ::    plot_tot_exch
   LOGICAL(LGT)    ::    plot_surf_area_ads
   LOGICAL(LGT)    ::    plot_surf_pot
   LOGICAL(LGT)    ::    plot_tot_surf
   LOGICAL(LGT)    ::    plot_minerals
   LOGICAL(LGT)    ::    plot_RSA_minerals
   LOGICAL(LGT)    ::    Rectangular
   LOGICAL(LGT)    ::    Rotational
   LOGICAL(LGT)    ::    Spherical
   LOGICAL(LGT)    ::    tecplot
   LOGICAL(LGT)    ::    TortuosityMapXEric
   LOGICAL(LGT)    ::    TortuosityMapYEric
   LOGICAL(LGT)    ::    Calculate_voltage
   
   REAL(DP)    ::    corrmax
   REAL(DP)    ::    dcoeffBulk
   REAL(DP)    ::    dcoeffDL
   REAL(DP)    ::    delt
   REAL(DP)    ::    deltmin
   REAL(DP)    ::    dtemp
   REAL(DP)    ::    dtmax
   REAL(DP)    ::    dtold
   REAL(DP)    ::    dzero
   REAL(DP)    ::    FixSaturation
   REAL(DP)    ::    fix_porosity
   REAL(DP)    ::    MinimumPorosity
   REAL(DP)    ::    OutputTimeScale
   REAL(DP)    ::    phwrite
   REAL(DP)    ::    PrintTime
   REAL(DP)    ::    rtolksp
   REAL(DP)    ::    time
   REAL(DP)    ::    tinit
   REAL(DP)    ::    time_scale_input
   REAL(DP)    ::    tstep
   REAL(DP)    ::    ttol
   REAL(DP)    ::    ttolmax
   !************* SWITCHES to accelerate functions evaluation with GPU
   REAL(DP)    ::    switch_helgeson
      
   !************* New counters to accelerate functions evaluation with GPU
   INTEGER(I4B)    ::    switch_temp_dependence
   INTEGER(I4B)    ::    switch_no_temp_dependence
   INTEGER(I4B)    ::    switch_igamma

   INTEGER(I4B)    ::    switch_Gaines_Thomas 
   INTEGER(I4B)    ::    switch_Vanselow      
   INTEGER(I4B)    ::    switch_Gapon         
   INTEGER(I4B)    ::    switch_DensityModule_temperature
   INTEGER(I4B)    ::    switch_DensityModule_sodium_chloride
   INTEGER(I4B)    ::    switch_DensityModule_sodium_nitrate
   INTEGER(I4B)    ::    switch_DensityModule_potassium_nitrate
   INTEGER(I4B)    ::    switch_DensityModule_calcium_nitrate
   INTEGER(I4B)    ::    switch_DensityModule_martin
   INTEGER(I4B)    ::    switch_activity_gradient
   
END TYPE SCALAR_V 



!***************************************************************************************************************

TYPE CONDITIONS_PARAMETERS

   CHARACTER (LEN=max_len_str)            :: condition_label
   CHARACTER (LEN=max_len_str)            :: condition_title
   LOGICAL(LGT)                           :: equilibrate_surface

END TYPE CONDITIONS_PARAMETERS

TYPE CELLS_PARAMETERS

   REAL(DP)                               :: porosity
   REAL(DP)                               :: porosity_init
   REAL(DP)                               :: volumetric_mass
   REAL(DP)                               :: tortuosity_x
   REAL(DP)                               :: tortuosity_y
   REAL(DP)                               :: tortuosity_z

END TYPE CELLS_PARAMETERS

TYPE POTENTIAL_PARAMETERS
    
   REAL(DP)                               :: Ie_East
   REAL(DP)                               :: Ie_North    
   REAL(DP)                               :: Ie_Up    
   REAL(DP)                               :: Id_East
   REAL(DP)                               :: Id_North    
   REAL(DP)                               :: Id_Up    
   REAL(DP)                               :: I_inject_pump    
    
END TYPE POTENTIAL_PARAMETERS   


TYPE SOLUTION_PARAMETERS

   INTEGER(I4B)                           :: Initial_condition_number
   INTEGER(I4B)                           :: boundary_inactive_active_cell
   INTEGER(I4B)                           :: SoliddensityFrom
   INTEGER(I4B)                           :: electrode
   REAL(DP)                               :: tempC
   REAL(DP)                               :: tempK
   REAL(DP)                               :: total_porosity
   REAL(DP)                               :: solid_volume_fraction  ! can be different from total_porosity depending on the way the solid density is defined
   REAL(DP)                               :: total_porosity_init
   REAL(DP)                               :: liquid_saturation
   REAL(DP)                               :: SolidSolutionRatio
   REAL(DP)                               :: Soliddensity
   REAL(DP)                               :: OneOverMassFraction
   REAL(DP)                               :: conversion 
   REAL(DP)                               :: average_volumetric_mass
   REAL(DP)                               :: forced_potential
   REAL(DP)                               :: forced_current

END TYPE SOLUTION_PARAMETERS




TYPE TEMPERATURE_DEPENDANCE_PARAMETERS

   REAL(DP)            :: as1

END TYPE TEMPERATURE_DEPENDANCE_PARAMETERS


TYPE AQUEOUS_SPECIES_GLOBAL_PARAMETERS

   CHARACTER (LEN=max_len_str)            :: name 
   REAL(DP)                               :: molar_mass
   REAL(DP)                               :: stoichio
   REAL(DP)                               :: charge
   REAL(DP)                               :: D0 !! 1 is for bulk, 2 is for DL1 etc. assume here that 1 + the number of DL is less than the number of primary species (should not be wrong most of time)
                                                !! Add an early check for that, OK done in startclay 
   REAL(DP)                               :: azero
   REAL(DP)                               :: bdotpar
   REAL(DP)                               :: lnkeqaq_standard
   
   INTEGER(I4B)                           :: Davies             ! 1/0
   INTEGER(I4B)                           :: DH_or_Helgeson     ! 1/0
   INTEGER(I4B)                           :: uncharged ! 1 if uncharged, 0 otherwise
   INTEGER(I4B)                           :: charged   ! 1 if charged  , 0 otherwise

END TYPE AQUEOUS_SPECIES_GLOBAL_PARAMETERS


TYPE AQUEOUS_SPECIES

   INTEGER(I4B)                           :: constraint  ! type of constraint: 1 = tot conc, 2 = charge etc.
   INTEGER(I4B)                           :: min_or_gas_num  ! mineral number or gas number depending on the constraint
   LOGICAL(LGT)                           :: equilibrate_sorption ! special type of constraint type 1
   REAL(DP)                               :: constraint_value ! tot conc, gas pressure, etc.
   REAL(DP)                               :: guess

END TYPE AQUEOUS_SPECIES


TYPE EXCHANGE_SPECIES_GLOBAL_PARAMETERS

   CHARACTER (LEN=max_len_str)            :: name 
   INTEGER(I4B)                           :: primlink  
   INTEGER(I4B)                           :: mineral_number
   INTEGER(I4B)                           :: convention ! to be implemented: having different convention for different exchange master species + add rothmund and Kornfeld??
   REAL(DP)                               :: stoichio
   REAL(DP)                               :: bexch
   REAL(DP)                               :: lnkeqexc_standard

END TYPE EXCHANGE_SPECIES_GLOBAL_PARAMETERS

TYPE EXCHANGE_SPECIES

   REAL(DP)                               :: CEC      ! in mol/m3
   REAL(DP)                               :: CECperkgw
   REAL(DP)                               :: CECperkgsolid
   REAL(DP)                               :: CECinput ! in mol per g (of mineral or of bulk) or in mol/kg water ! THIS COULD BE SIMPLIFIED BY ELIMINATING THIS VARIABLE IN READ_IONEXCHANGE...
   INTEGER(I4B)                           :: CEC_bulk_or_mineral
   INTEGER(I4B)                           :: nb_exchange_mineral 
   INTEGER(I4B)                           :: num_exchange_mineral
   INTEGER(I4B)                           :: nb_exchange_bulk 
   INTEGER(I4B)                           :: num_exchange_bulk

END TYPE EXCHANGE_SPECIES



TYPE SURFACE_FOR_SURFACE_COMPLEXATION_GLOBAL_PARAMETERS

   CHARACTER (LEN=max_len_str)            :: name 
   INTEGER(I4B)                           :: mineral_number 
   INTEGER(I4B)                           :: associated_DLnum  
   INTEGER(I4B)                           :: model_edl
   INTEGER(I4B)                           :: switch_DzombakMorel
   INTEGER(I4B)                           :: switch_no_edl
   INTEGER(I4B)                           :: switch_Kd
   
END TYPE SURFACE_FOR_SURFACE_COMPLEXATION_GLOBAL_PARAMETERS


TYPE SURFACE_FOR_SURFACE_COMPLEXATION

   REAL(DP)                               :: complexationRSA             ! m2/m3
   REAL(DP)                               :: ssa                         ! m2/g

END TYPE SURFACE_FOR_SURFACE_COMPLEXATION

TYPE SURFACE_SPECIES_GLOBAL_PARAMETERS

   CHARACTER (LEN=max_len_str)            :: name 
   INTEGER(I4B)                           :: surface_num
   REAL(DP)                               :: site_density  
   REAL(DP)                               :: lnkeqsurf_standard
   REAL(DP)                               :: stoichio ! relative to ncomp, nexch, nsurf
   REAL(DP)                               :: charge
   LOGICAL(LGT)                           :: KD
   REAL(DP)                               :: KDvalue
   INTEGER(I4B)                           :: icompKD
   INTEGER(I4B)                           :: switch_KD

END TYPE SURFACE_SPECIES_GLOBAL_PARAMETERS


TYPE MINERAL_SPECIES_GLOBAL_PARAMETERS

   CHARACTER (LEN=max_len_str)            :: name                 !
   CHARACTER (LEN=max_len_str)            :: rlabel               !
   CHARACTER (LEN=max_len_str)            :: crossaff             !
   CHARACTER (LEN=max_len_str)            :: link_surface         !
   INTEGER(I4B)                           :: number_of_kinetics   !
   INTEGER(I4B)                           :: type_of_kinetics     !
   INTEGER(I4B)                           :: kcrossaff            !
   INTEGER(I4B)                           :: klink_surface        !
   REAL(DP)                               :: klink_surface_factor ! g/mol
   REAL(DP)                               :: molar_mass           ! g/mol
   REAL(DP)                               :: molar_volume         ! cm3/mol in database file converted in m3/mol
   REAL(DP)                               :: stoichio             !
   REAL(DP)                               :: lnKeqmin_standard    !
   REAL(DP)                               :: rate0                !
   REAL(DP)                               :: thresh               !
   REAL(DP)                               :: activation_energy    !
   
   LOGICAL(LGT)                           :: LocalEquilibrium     !

   REAL(DP)                               :: AffinityDepend1      !
   REAL(DP)                               :: AffinityDepend2      !
   REAL(DP)                               :: AffinityDepend3      !
                                                                  !
   REAL(DP)                               :: halfsat              !
   INTEGER(I4B)                           :: imonod               !
   INTEGER(I4B)                           :: kmonod               !
   INTEGER(I4B)                           :: itot_monod           !
   INTEGER(I4B)                           :: nmonod               !
                                                                  !
   REAL(DP)                               :: rinhibit             !
   INTEGER(I4B)                           :: inhibit              !
   INTEGER(I4B)                           :: kinhibit             !
   INTEGER(I4B)                           :: itot_inhibit         !
   INTEGER(I4B)                           :: ninhibit             !

   ! dependence on primary component (total or not) concentration
   CHARACTER (LEN=max_len_str)            :: namdep_nyf
   INTEGER(I4B)                           :: ndepend
   INTEGER(I4B)                           :: ndependex
   INTEGER(I4B)                           :: ndependsurf
   INTEGER(I4B)                           :: itot_min
   REAL(DP)                               :: depend
   INTEGER(I4B)                           :: idepend
   REAL(DP)                               :: dependex
   INTEGER(I4B)                           :: ixdepend
   REAL(DP)                               :: dependsurf
   INTEGER(I4B)                           :: isdepend

   ! Switch to better parallelize the calculations 
   INTEGER(I4B)                           :: switch_monod         !
   INTEGER(I4B)                           :: switch_diss_only     !
   INTEGER(I4B)                           :: switch_precip_only   !
   INTEGER(I4B)                           :: switch_irreversible  !
   INTEGER(I4B)                           :: switch_cross_affinity!
   INTEGER(I4B)                           :: switch_tst           !
   INTEGER(I4B)                           :: switch_monod_self       !
   INTEGER(I4B)                           :: switch_monod_other_min  !
   INTEGER(I4B)                           :: switch_monod_other_surf !
   INTEGER(I4B)                           :: switch_monod_tot_aq     !
   INTEGER(I4B)                           :: switch_monod_conc_aq    !
   INTEGER(I4B)                           :: switch_inhib_other_min  !
   INTEGER(I4B)                           :: switch_inhib_other_surf !
   INTEGER(I4B)                           :: switch_inhib_tot_aq     !
   INTEGER(I4B)                           :: switch_inhib_conc_aq    !  
   
END TYPE MINERAL_SPECIES_GLOBAL_PARAMETERS


TYPE MINERAL_SPECIES

   INTEGER(I4B)                           :: surface_area_option
   REAL(DP)                               :: vol_fraction
   REAL(DP)                               :: vol_fraction_init
   REAL(DP)                               :: threshold_volume_fraction
   REAL(DP)                               :: RSA !  in m2/m3 to be updated in mineral_update = arr4Dm(k,jx,jy,jz)%area, only if not bulk surface area
   REAL(DP)                               :: RSA_init 
   REAL(DP)                               :: SSA   ! in m2/g
   REAL(DP)                               :: SSA_init 
   REAL(DP)                               :: dppt

END TYPE MINERAL_SPECIES


TYPE GAS_SPECIES_GLOBAL_PARAMETERS

   CHARACTER (LEN=max_len_str)            :: name 
   REAL(DP)                               :: stoichio
   REAL(DP)                               :: lnKeqgas_standard

END TYPE GAS_SPECIES_GLOBAL_PARAMETERS

TYPE GAS_SPECIES

   REAL(DP)                               :: pressure

END TYPE GAS_SPECIES

TYPE AQUEOUS_KINETICS_GLOBAL_PARAMETERS

   CHARACTER (LEN=max_len_str)            :: name 
   INTEGER(I4B)                           :: nreactkin
   INTEGER(I4B)                           :: reaction_type
   INTEGER(I4B)                           :: totconcflag
   INTEGER(I4B)                           :: itot_monodaq
   INTEGER(I4B)                           :: nmonodaq
   INTEGER(I4B)                           :: imonodaq
   INTEGER(I4B)                           :: itot_inhibitaq
   INTEGER(I4B)                           :: ninhibitaq
   INTEGER(I4B)                           :: inhibitaq
   REAL(DP)                               :: stoichio
   REAL(DP)                               :: ratek
   REAL(DP)                               :: lnKeqkin_standard
   REAL(DP)                               :: dependk
   REAL(DP)                               :: halfsataq
   REAL(DP)                               :: rinhibitaq

END TYPE AQUEOUS_KINETICS_GLOBAL_PARAMETERS

TYPE DISCRETIZATION_ARRAY

   REAL(DP)                               :: center
   REAL(DP)                               :: size
   REAL(DP)                               :: distance_next
   REAL(DP)                               :: distance_prev

END TYPE DISCRETIZATION_ARRAY

TYPE ELECT_ARRAY

   CHARACTER (LEN=max_len_str)            :: name 
   CHARACTER (LEN=max_len_str)            :: nameworking 
   CHARACTER (LEN=max_len_str)            :: namecounter 
   INTEGER(I4B)                           :: jx
   INTEGER(I4B)                           :: jy
   INTEGER(I4B)                           :: jz
   REAL(DP)                               :: fixed_potential
   REAL(DP)                               :: fixed_current
   LOGICAL(LGT)                           :: Found
   INTEGER(I4B)                           :: typeelectrode ! 1- Working electrode with fixed_current
                                                           ! 2- Counter electrode associated to working electrode with fixed current
                                                           ! 3- Ground electrode (only one electrode in the system)
                                                           ! Working electrode with fixed_potential
                                                           ! ...

END TYPE ELECT_ARRAY




TYPE RECEIVE_AND_SEND

   LOGICAL(LGT)                           :: Self
   LOGICAL(LGT)                           :: Sender
   LOGICAL(LGT)                           :: Receiver
   INTEGER(I4B)                           :: pos_vector
   INTEGER(I4B)                           :: pos_cell
   INTEGER(I4B)                           :: proc_sender
   INTEGER(I4B)                           :: proc_receiver

END TYPE RECEIVE_AND_SEND

TYPE RECEIVE_AND_SEND_GHOST

   LOGICAL(LGT)                           :: Sender
   LOGICAL(LGT)                           :: Receiver
   INTEGER(I4B)                           :: global_cell
   INTEGER(I4B)                           :: pos_cell_sender
   INTEGER(I4B)                           :: proc_sender
   INTEGER(I4B)                           :: rank_receiver       ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: pos_cell_receiver_1 ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: pos_cell_receiver_2 ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: pos_cell_receiver_3 ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: pos_cell_receiver_4 ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: pos_cell_receiver_5 ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: pos_cell_receiver_6 ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: pos_cell_receiver_7 ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: proc_receiver_1     ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: proc_receiver_2     ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: proc_receiver_3     ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: proc_receiver_4     ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: proc_receiver_5     ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: proc_receiver_6     ! For corners there can be seven proc that receive the information
   INTEGER(I4B)                           :: proc_receiver_7     ! For corners there can be seven proc that receive the information

END TYPE RECEIVE_AND_SEND_GHOST






TYPE PROCESSOR_OWNERSHIP 

   INTEGER(I4B)                           :: xbegin
   INTEGER(I4B)                           :: xend
   INTEGER(I4B)                           :: xghbegin
   INTEGER(I4B)                           :: xghend
   LOGICAL(LGT)                           :: xbeginghosted
   LOGICAL(LGT)                           :: xendghosted
   INTEGER(I4B)                           :: ybegin
   INTEGER(I4B)                           :: yend
   INTEGER(I4B)                           :: yghbegin
   INTEGER(I4B)                           :: yghend
   LOGICAL(LGT)                           :: ybeginghosted
   LOGICAL(LGT)                           :: yendghosted
   INTEGER(I4B)                           :: zbegin
   INTEGER(I4B)                           :: zend
   INTEGER(I4B)                           :: zghbegin
   INTEGER(I4B)                           :: zghend
   LOGICAL(LGT)                           :: zbeginghosted
   LOGICAL(LGT)                           :: zendghosted
   INTEGER(I4B)                           :: neigh_proc_E
   INTEGER(I4B)                           :: neigh_proc_W
   INTEGER(I4B)                           :: neigh_proc_N
   INTEGER(I4B)                           :: neigh_proc_S
   INTEGER(I4B)                           :: neigh_proc_U
   INTEGER(I4B)                           :: neigh_proc_D
   INTEGER(I4B)                           :: nbxMPI
   INTEGER(I4B)                           :: nbyMPI
   INTEGER(I4B)                           :: nbzMPI
   INTEGER(I4B)                           :: nbxMPI_noghost
   INTEGER(I4B)                           :: nbyMPI_noghost
   INTEGER(I4B)                           :: nbzMPI_noghost
   INTEGER(I4B)                           :: nbcellMPI
   INTEGER(I4B)                           :: nbcellMPI_noghost
   INTEGER(I4B)                           :: residual_ncell_begin
   INTEGER(I4B)                           :: residual_ncell_end
   INTEGER(I4B)                           :: residual_nbcell
   INTEGER(I4B)                           :: rowmin
   INTEGER(I4B)                           :: rowmax


END TYPE PROCESSOR_OWNERSHIP 


TYPE Map_Cells_ARRAY 

   INTEGER(I4B)                           :: jx
   INTEGER(I4B)                           :: jy
   INTEGER(I4B)                           :: jz
   INTEGER(I4B)                           :: jxref
   INTEGER(I4B)                           :: jyref
   INTEGER(I4B)                           :: jzref
   INTEGER(I4B)                           :: ncell_ref
   INTEGER(I4B)                           :: Ghostcell
   INTEGER(I4B)                           :: switch_N_current_no_charge_build_up
   INTEGER(I4B)                           :: switch_N_current_equate_potentials
   INTEGER(I4B)                           :: switch_U_current_no_charge_build_up
   INTEGER(I4B)                           :: switch_U_current_equate_potentials
   INTEGER(I4B)                           :: North_current_N_E_S_W
   INTEGER(I4B)                           :: Up_current_U_E_D_W
   INTEGER(I4B)                           :: Up_current_U_N_D_S
   INTEGER(I4B)                           :: E
   INTEGER(I4B)                           :: W
   INTEGER(I4B)                           :: N
   INTEGER(I4B)                           :: S
   INTEGER(I4B)                           :: U
   INTEGER(I4B)                           :: D
   INTEGER(I4B)                           :: NE
   INTEGER(I4B)                           :: UE
   INTEGER(I4B)                           :: UN
   INTEGER(I4B)                           :: next_to_E_bc
   INTEGER(I4B)                           :: next_to_W_bc
   INTEGER(I4B)                           :: next_to_N_bc
   INTEGER(I4B)                           :: next_to_S_bc
   INTEGER(I4B)                           :: next_to_U_bc
   INTEGER(I4B)                           :: next_to_D_bc
   INTEGER(I4B)                           :: E_bound
   INTEGER(I4B)                           :: W_bound
   INTEGER(I4B)                           :: N_bound
   INTEGER(I4B)                           :: S_bound
   INTEGER(I4B)                           :: U_bound
   INTEGER(I4B)                           :: D_bound
   INTEGER(I4B)                           :: Edge
   INTEGER(I4B)                           :: Corner
   INTEGER(I4B)                           :: Calc_East_Current
   INTEGER(I4B)                           :: Calc_North_Current
   INTEGER(I4B)                           :: Calc_Up_Current
   INTEGER(I4B)                           :: Calc_E_Diffusion 
   INTEGER(I4B)                           :: Calc_W_Diffusion 
   INTEGER(I4B)                           :: Calc_N_Diffusion 
   INTEGER(I4B)                           :: Calc_S_Diffusion 
   INTEGER(I4B)                           :: Calc_U_Diffusion 
   INTEGER(I4B)                           :: Calc_D_Diffusion 
   INTEGER(I4B)                           :: Constant_potential 
   INTEGER(I4B)                           :: Fixed_Current
   
   REAL(DP)                               :: x_center
   REAL(DP)                               :: y_center
   REAL(DP)                               :: z_center
   REAL(DP)                               :: dx
   REAL(DP)                               :: dy
   REAL(DP)                               :: dz
   REAL(DP)                               :: x_begin
   REAL(DP)                               :: y_begin
   REAL(DP)                               :: z_begin
   REAL(DP)                               :: x_end
   REAL(DP)                               :: y_end
   REAL(DP)                               :: z_end
   REAL(DP)                               :: surface_E
   REAL(DP)                               :: surface_W
   REAL(DP)                               :: surface_N
   REAL(DP)                               :: surface_S
   REAL(DP)                               :: surface_U
   REAL(DP)                               :: surface_D
   REAL(DP)                               :: volume   
   REAL(DP)                               :: Constant_Potential_Value
   REAL(DP)                               :: Fixed_Current_Value
   
   
END TYPE Map_Cells_ARRAY


TYPE DEBYE_HUCKEL_PARAMETERS

   REAL(DP)                               :: tempC
   REAL(DP)                               :: dh_a
   REAL(DP)                               :: dh_b
   REAL(DP)                               :: bdot
   REAL(DP)                               :: co2_coefs
   REAL(DP)                               :: dh_aPolynomial 
   REAL(DP)                               :: dh_bPolynomial 
   REAL(DP)                               :: bdotPolynomial 

END TYPE DEBYE_HUCKEL_PARAMETERS


!******************************************************************************************************************************************
!********* TYPES CREATED TO READ AND TRACK INFORMATION IN FILES ***************************************************************************

TYPE LIST_NAMES
   CHARACTER(LEN=max_len_str)           :: ListNames 
END TYPE LIST_NAMES


TYPE WORD_OR_NUMBER_IN_FILE

   INTEGER(I4B)                           :: StrLength
   CHARACTER(LEN=max_len_str)             :: WordOrNumber 
   LOGICAL(LGT)                           :: IsWord
   LOGICAL(LGT)                           :: IsNumber
   LOGICAL(LGT)                           :: IsInteger
   LOGICAL(LGT)                           :: IsComment
   LOGICAL(LGT)                           :: IsHyphRange
   INTEGER(I4B)                           :: IntValue
   REAL(DP)                               :: FloatValue
   
END TYPE WORD_OR_NUMBER_IN_FILE  


TYPE LINE_IN_FILE

   CHARACTER(LEN=max_len_str*10)                            :: FullLine 
   CHARACTER(LEN=max_len_str*10)                            :: CleanLine 
   INTEGER(I4B)                                             :: PosLineInFile
   CHARACTER(LEN=max_len_str)                               :: NameOfBlock 
   INTEGER(I4B)                                             :: PosLineInBlock
   INTEGER(I4B)                                             :: NumberOfItems
   TYPE(WORD_OR_NUMBER_IN_FILE), DIMENSION(:), ALLOCATABLE  :: LineWordsAndNumbers
   LOGICAL(LGT)                                             :: Empty
   LOGICAL(LGT)                                             :: Comment

   
END TYPE LINE_IN_FILE  

TYPE BLOCK_CONTENT

   LOGICAL(LGT)                                             :: Exists
   INTEGER(I4B)                                             :: PosBlockInFile
   INTEGER(I4B)                                             :: NbLinesInBlock
   CHARACTER(LEN=max_len_str)                               :: BlockDescription 
   TYPE(LINE_IN_FILE), DIMENSION(:), ALLOCATABLE            :: BlockLine
   
END TYPE BLOCK_CONTENT  

TYPE FILE_CONTENT

   INTEGER(I4B)                                             :: NbLinesInFile
   INTEGER(I4B)                                             :: NbLinesInFileContent
   TYPE(LINE_IN_FILE), DIMENSION(:), ALLOCATABLE            :: FileLine
   
END TYPE FILE_CONTENT  

TYPE SOLUTION_PRIMARY_SPECIES

   CHARACTER(LEN=max_len_str)           :: element 
   CHARACTER(LEN=max_len_str)           :: species 
   CHARACTER(LEN=max_len_str)           :: gfw_formula 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: ElementStoichio 
   REAL(DP)                             :: llnl_gamma 
   REAL(DP)                             :: alk 
   REAL(DP)                             :: element_gfw 
   REAL(DP)                             :: charge 
   REAL(DP)                             :: species_gfw 
   REAL(DP)                             :: RedoxState 
   INTEGER(I4B)                         :: SecondaryIsPrimary
   REAL(DP)                             :: Vm 

END TYPE SOLUTION_PRIMARY_SPECIES

TYPE SOLUTION_SECONDARY_SPECIES

   CHARACTER(LEN=max_len_str)           :: species 
   CHARACTER(LEN=max_len_str)           :: gfw_formula 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: PrimStoichio 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: ElementStoichio 
   REAL(DP)                             :: RedoxState 
   REAL(DP), DIMENSION(8)               :: log_kT 
   REAL(DP), DIMENSION(5)               :: analytical_logk 
   REAL(DP), DIMENSION(5)               :: LogKCrunchPolynomial 
   REAL(DP)                             :: alk 
   REAL(DP)                             :: charge 
   REAL(DP)                             :: species_gfw 
   REAL(DP)                             :: log_k25 
   REAL(DP)                             :: delta_h 
   REAL(DP)                             :: llnl_gamma 
   REAL(DP)                             :: electronstoichio 
   REAL(DP)                             :: Vm 
   LOGICAL(LGT)                         :: IsAnalyticallogk
   LOGICAL(LGT)                         :: IsLlnllogk
   LOGICAL(LGT)                         :: IsVantHoff
   INTEGER(I4B)                         :: PrimaryIsSecondary
   INTEGER(I4B)                         :: NbReactants


END TYPE SOLUTION_SECONDARY_SPECIES

TYPE SOLID_PHASES

   CHARACTER(LEN=max_len_str)           :: name 
   CHARACTER(LEN=max_len_str)           :: species 
   CHARACTER(LEN=max_len_str)           :: gfw_formula 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: PrimStoichio 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: ElementStoichio 
   REAL(DP), DIMENSION(8)               :: log_kT 
   REAL(DP), DIMENSION(5)               :: analytical_logk 
   REAL(DP), DIMENSION(5)               :: LogKCrunchPolynomial 
   REAL(DP)                             :: charge 
   REAL(DP)                             :: species_gfw 
   REAL(DP)                             :: log_k25 
   REAL(DP)                             :: delta_h 
   REAL(DP)                             :: electronstoichio 
   REAL(DP)                             :: Vm 
   LOGICAL(LGT)                         :: IsAnalyticallogk
   LOGICAL(LGT)                         :: IsLlnllogk
   LOGICAL(LGT)                         :: IsVantHoff

END TYPE SOLID_PHASES
 
TYPE SURFACE_PRIMARY_SPECIES

   CHARACTER(LEN=max_len_str)           :: name 
   CHARACTER(LEN=max_len_str)           :: species 
   REAL(DP)                             :: charge 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: ElementStoichio 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: SurfaceStoichio 
   CHARACTER(LEN=max_len_str)           :: NameOfSurface                    
   REAL(DP)                             :: SiteDensity            
   INTEGER(I4B)                         :: SurfaceId

END TYPE SURFACE_PRIMARY_SPECIES

TYPE SURFACE_SECONDARY_SPECIES

   CHARACTER(LEN=max_len_str)           :: species 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: PrimStoichio 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: SurfPrimStoichio 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: ElementStoichio 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: SurfaceStoichio 
   REAL(DP), DIMENSION(8)               :: log_kT 
   REAL(DP), DIMENSION(5)               :: analytical_logk 
   REAL(DP), DIMENSION(5)               :: LogKCrunchPolynomial 
   REAL(DP)                             :: charge 
   REAL(DP)                             :: log_k25 
   REAL(DP)                             :: delta_h 
   REAL(DP)                             :: electronstoichio 
   LOGICAL(LGT)                         :: IsAnalyticallogk
   LOGICAL(LGT)                         :: IsLlnllogk
   LOGICAL(LGT)                         :: IsVantHoff

END TYPE SURFACE_SECONDARY_SPECIES

TYPE EXCHANGE_PRIMARY_SPECIES

   CHARACTER(LEN=max_len_str)           :: Name 
   CHARACTER(LEN=max_len_str)           :: Species 
   REAL(DP)                             :: Charge 
   CHARACTER(LEN=max_len_str)           :: Convention                    
   CHARACTER(LEN=max_len_str)           :: Mineral                    
   INTEGER(I4B)                         :: MineralId

END TYPE EXCHANGE_PRIMARY_SPECIES

TYPE EXCHANGE_SECONDARY_SPECIES

   CHARACTER(LEN=max_len_str)           :: species 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: PrimStoichio 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: ExchPrimStoichio 
   REAL(DP), DIMENSION(:), ALLOCATABLE  :: ElementStoichio 
   REAL(DP), DIMENSION(8)               :: log_kT 
   REAL(DP), DIMENSION(5)               :: analytical_logk 
   REAL(DP), DIMENSION(5)               :: LogKCrunchPolynomial 
   REAL(DP)                             :: charge 
   REAL(DP)                             :: log_k25 
   REAL(DP)                             :: delta_h 
   REAL(DP)                             :: bexch 
   REAL(DP)                             :: electronstoichio 
   LOGICAL(LGT)                         :: IsAnalyticallogk
   LOGICAL(LGT)                         :: IsLlnllogk
   LOGICAL(LGT)                         :: IsVantHoff
   LOGICAL(LGT)                         :: Isbexch

END TYPE EXCHANGE_SECONDARY_SPECIES


TYPE ELEMENT

  CHARACTER(LEN=max_len_str)           :: name
  REAL(DP)                             :: MolarMass 
   
END TYPE ELEMENT   


TYPE MINERAL_KINETICS

   CHARACTER(LEN=max_len_str)                          :: MineralName 
   LOGICAL(LGT)                                        :: IsPresent
   INTEGER(I4B)                                        :: NbLabels
   TYPE(MINERAL_KINETICS_SUB),DIMENSION(:),ALLOCATABLE :: KinRateLaw
   
END TYPE MINERAL_KINETICS

TYPE MINERAL_KINETICS_SUB

   CHARACTER(LEN=max_len_str)                  :: Label 
   REAL(DP)                                    :: Rate 
   REAL(DP), DIMENSION(3)                      :: FreeEnergyDependence 
   REAL(DP)                                    :: ActEner 
   INTEGER(I4B)                                :: NbDependence
   TYPE(LIST_NAMES), DIMENSION(:), ALLOCATABLE :: DependenceName
   REAL(DP), DIMENSION(:), ALLOCATABLE         :: DependenceValue 
   INTEGER(I4B), DIMENSION(:), ALLOCATABLE     :: DependenceType ! 1: aq species, 2: aq tot, 3: surf, 4: min 
   INTEGER(I4B), DIMENSION(:), ALLOCATABLE     :: DependenceSpeciesId ! Id of the species 
   INTEGER(I4B)                                :: NbMonod
   TYPE(LIST_NAMES), DIMENSION(:), ALLOCATABLE :: MonodName
   REAL(DP), DIMENSION(:), ALLOCATABLE         :: MonodValue 
   INTEGER(I4B), DIMENSION(:), ALLOCATABLE     :: MonodType ! 1: aq species, 2: aq tot, 3: surf, 4: min 
   INTEGER(I4B), DIMENSION(:), ALLOCATABLE     :: MonodSpeciesId ! Id of the species 
   INTEGER(I4B)                                :: NbInhibition
   TYPE(LIST_NAMES), DIMENSION(:), ALLOCATABLE :: InhibitionName
   REAL(DP), DIMENSION(:), ALLOCATABLE         :: InhibitionValue 
   INTEGER(I4B), DIMENSION(:), ALLOCATABLE     :: InhibitionType ! 1: aq species, 2: aq tot, 3: surf, 4: min 
   INTEGER(I4B), DIMENSION(:), ALLOCATABLE     :: InhibitionSpeciesId ! Id of the species 
   LOGICAL(LGT)                                :: ConstantSurface
   LOGICAL(LGT)                                :: Threshold
   REAL(DP)                                    :: ThresholdValue 
   LOGICAL(LGT)                                :: CrossAffinity
   CHARACTER(LEN=max_len_str)                  :: CrossAffinityName
   INTEGER(I4B)                                :: CrossAffinitySpeciesId 
   LOGICAL(LGT)                                :: LinkSurface
   CHARACTER(LEN=max_len_str)                  :: LinkSurfaceName
   REAL(DP)                                    :: LinkSurfaceFraction 
   INTEGER(I4B)                                :: LinkSurfaceSpeciesId 
   LOGICAL(LGT)                                :: DissolutionOnly
   LOGICAL(LGT)                                :: PrecipitationOnly
   LOGICAL(LGT)                                :: Irreversible

   LOGICAL(LGT)                                :: NewRate 
   LOGICAL(LGT)                                :: NewFreeEnergyDependence 
   LOGICAL(LGT)                                :: NewActEner 
   LOGICAL(LGT)                                :: NewNbDependence
   LOGICAL(LGT)                                :: NewNbMonod
   LOGICAL(LGT)                                :: NewNbInhibition
   LOGICAL(LGT)                                :: NewConstantSurface
   LOGICAL(LGT)                                :: NewThreshold
   LOGICAL(LGT)                                :: NewCrossAffinity
   LOGICAL(LGT)                                :: NewLinkSurface
   LOGICAL(LGT)                                :: NewDissolutionOnly
   LOGICAL(LGT)                                :: NewPrecipitationOnly
   LOGICAL(LGT)                                :: NewIrreversible
   
END TYPE MINERAL_KINETICS_SUB

TYPE AVAILABLE_SURFACES

   CHARACTER(LEN=max_len_str)           :: Name 
   CHARACTER(LEN=max_len_str)           :: Mineral
   INTEGER(I4B)                         :: MineralId 
   INTEGER(I4B)                         :: DLId 
   INTEGER(I4B)                         :: EDLModel 
                    
END TYPE AVAILABLE_SURFACES

TYPE KD_SPECIES

   CHARACTER(LEN=max_len_str)           :: Species 
   CHARACTER(LEN=max_len_str)           :: SurfaceName 
   INTEGER(I4B)                         :: SurfaceId
   REAL(DP)                             :: KdValue            

END TYPE KD_SPECIES




TYPE GLOBAL_PARAMETERS_BUNDLE

  TYPE(TEMPERATURE_DEPENDANCE_PARAMETERS)                 ,DIMENSION(:,:)   ,ALLOCATABLE   :: TempDepend
  TYPE(DEBYE_HUCKEL_PARAMETERS)                           ,DIMENSION(:)     ,ALLOCATABLE   :: DHParam
  TYPE(AQUEOUS_SPECIES_GLOBAL_PARAMETERS)                 ,DIMENSION(:,:)   ,ALLOCATABLE   :: AqSpecies 
  TYPE(MINERAL_SPECIES_GLOBAL_PARAMETERS)                 ,DIMENSION(:,:,:) ,ALLOCATABLE   :: Minerals 
  TYPE(EXCHANGE_SPECIES_GLOBAL_PARAMETERS)                ,DIMENSION(:,:)   ,ALLOCATABLE   :: ExchSpecies 
  TYPE(GAS_SPECIES_GLOBAL_PARAMETERS)                     ,DIMENSION(:,:)   ,ALLOCATABLE   :: GasSpecies 
  TYPE(SURFACE_SPECIES_GLOBAL_PARAMETERS)                 ,DIMENSION(:,:)   ,ALLOCATABLE   :: SurfSpecies 
  TYPE(SURFACE_FOR_SURFACE_COMPLEXATION_GLOBAL_PARAMETERS),DIMENSION(:)     ,ALLOCATABLE   :: ComplexSurf 
  TYPE(AQUEOUS_KINETICS_GLOBAL_PARAMETERS)                ,DIMENSION(:,:,:) ,ALLOCATABLE   :: AqKinetics 
  INTEGER(I4B)                                            ,DIMENSION(:)     ,ALLOCATABLE   :: MeanSalt
  REAL(DP)                                                ,DIMENSION(:)     ,ALLOCATABLE   :: DBSTemp ! to be suppressed later 
  
END TYPE GLOBAL_PARAMETERS_BUNDLE





END MODULE TYPES_AND_PARAMETERS
