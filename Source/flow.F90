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
  
!! December 31, 2022:  Fixed flow.F90

!!!      ****************************************

MODULE flow

    USE crunchtype
    
    REAL(DP)                                     :: permmaxX
    REAL(DP)                                     :: permmaxY
    REAL(DP)                                     :: permmaxZ
    REAL(DP)                                     :: permminX
    REAL(DP)                                     :: permminY
    REAL(DP)                                     :: permminZ
    REAL(DP)                                     :: maxQx
    REAL(DP)                                     :: maxQy
    REAL(DP)                                     :: maxQz
    REAL(DP)                                     :: x_angle
    REAL(DP)                                     :: y_angle
    REAL(DP)                                     :: z_angle
    REAL(DP)                                     :: SignGravity
    
    LOGICAL(LGT)                                 :: CalculateFlow
    LOGICAL(LGT)                                 :: Neumann
    LOGICAL(LGT)                                 :: InitializeHydrostatic
    
    ! NavierStokes option added by Zhi Li
    LOGICAL(LGT)                                 :: NavierStokes
    
    INTEGER(I4B)                                 :: n_substep
    
    INTEGER(I4B)                                 :: infiltration
    
    
    INTEGER(I4B), DIMENSION(:,:,:,:), ALLOCATABLE  :: intbnd
    INTEGER(I4B), DIMENSION(:,:,:,:), ALLOCATABLE  :: intbndgas
    
    REAL(DP), DIMENSION(:,:), ALLOCATABLE          :: qrecharge
    
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE      :: qg
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE      :: gaspump
    
    LOGICAL(LGT)                                   :: wells
    LOGICAL(LGT)                                   :: ThereArePumpingWells
    
    REAL(DP), DIMENSION(:),ALLOCATABLE           :: permzonex
    REAL(DP), DIMENSION(:),ALLOCATABLE           :: permzoney
    REAL(DP), DIMENSION(:),ALLOCATABLE           :: permzonez
    
    INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE  :: activecellPressure
    
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermx_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermx_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermx_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermx_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermx_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermx_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermy_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermy_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermy_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermy_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermy_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermy_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermz_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxpermz_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermz_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyypermz_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermz_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzpermz_hi
    
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxvgn_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxvgn_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyyvgn_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyyvgn_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzvgn_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzvgn_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxvga_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxvga_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyyvga_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyyvga_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzvga_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzvga_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxwcr_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxwcr_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyywcr_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyywcr_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzwcr_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzwcr_hi
    
    REAL(DP), DIMENSION(:), ALLOCATABLE          :: PressureZone
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: PressureFix
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxPressure_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jxxPressure_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyyPressure_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jyyPressure_hi
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzPressure_lo
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE      :: jzzPressure_hi
    
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE       :: pres
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: harx
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: hary
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: harz
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: perminx
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: perminy
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: perminz
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permx
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permy
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permz
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permxOld
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permyOld
    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE         :: permzOld
    
    INTEGER(I4B), DIMENSION(:,:,:),ALLOCATABLE     :: npump
    
    ! *************************************************
    ! For Richards solver by Toshiyuki Bandai May, 2023
    LOGICAL(LGT)                                    :: Richards ! When solving the Richards equation with Toshi's code
    
    ! soil hydraulic parameters for van Genuchten model
    REAL(DP), DIMENSION(:),     ALLOCATABLE         :: VG_params_zone ! array to store VG parameters when reading input file
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE         :: jxx_VG_params_lo ! array to store the location of VG parameters when reading input file
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE         :: jxx_VG_params_hi ! array to store the location of VG parameters when reading input file
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE         :: jyy_VG_params_lo ! array to store the location of VG parameters when reading input file
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE         :: jyy_VG_params_hi ! array to store the location of VG parameters when reading input file
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE         :: jzz_VG_params_lo ! array to store the location of VG parameters when reading input file
    INTEGER(I4B), DIMENSION(:), ALLOCATABLE         :: jzz_VG_params_hi ! array to store the location of VG parameters when reading input file
        
    ! End of edits by Toshiyuki Bandai Oct, 2024
    ! *************************************************
    
    !!  PETSc arrays for solver
    
    REAL(DP), DIMENSION(:),ALLOCATABLE           :: XvecCrunchP
    REAL(DP), DIMENSION(:),ALLOCATABLE           :: BvecCrunchP
    
    !! Pump time series variable Lucien Stolze 20211022
    REAL(DP), DIMENSION(:), ALLOCATABLE          :: qgt
    REAL(DP), DIMENSION(:), ALLOCATABLE             :: tpump
    LOGICAL(LGT)                                   :: pumptimeseries
    
    !!!Evapotranspiration timeseries
    !REAL(DP), DIMENSION(:), ALLOCATABLE          :: qt_evapo
    !REAL(DP), DIMENSION(:), ALLOCATABLE             :: t_evapo
    !REAL(DP), DIMENSION(:), ALLOCATABLE          :: qt_transpi
    !REAL(DP), DIMENSION(:), ALLOCATABLE             :: t_transpi
    !! Lucien Stolze: evapotranspiration
    !LOGICAL(LGT)                                     :: evapofix !fixed evaporation
    !LOGICAL(LGT)                                     :: evapotimeseries !transient evaporation
    !LOGICAL(LGT)                                     :: transpifix !fixed transpiration
    !LOGICAL(LGT)                                     :: transpitimeseries !transient transpiration
    !REAL(DP)                                         :: evaporate !evaporation rate
    !REAL(DP)                                         :: transpirate !transpiration rate
    !INTEGER(I4B)                                     :: transpicells !nb of transpiration cells
    !! infiltration timeseries added by Toshiyuki Bandai June, 2023
    !LOGICAL(LGT)                                     :: infiltration_timeseries
    !LOGICAL(LGT)                                     :: infiltration_fix
    !REAL(DP)                                         :: infiltration_rate
    !REAL(DP), DIMENSION(:), ALLOCATABLE              :: qt_infiltration
    !REAL(DP), DIMENSION(:), ALLOCATABLE              :: t_infiltration
    !
    !REAL(DP), DIMENSION(:), ALLOCATABLE              :: transpirate_cell
    
    END MODULE flow
    