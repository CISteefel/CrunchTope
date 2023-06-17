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

MODULE temperature

  USE crunchtype
  USE params

  CHARACTER (LEN=mls)                              :: Tfile

  INTEGER(I4B)                                     :: jtemp
  INTEGER(I4B)                                     :: ntemp
  INTEGER(I4B)                                     :: TPointer
  
  LOGICAL(LGT)                                     :: RunIsothermal
    
  REAL(DP)                                         :: tinit
  REAL(DP)                                         :: tgrad

  REAL(DP), DIMENSION(:), ALLOCATABLE              :: tempcond
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: rocond

!  Allocatable arrays dimensioned over spatial domain

!fp! block(1);
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: t
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: told
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: ro
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: roOld

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: rogas
!fp! darray=0;

! Temperature time series (Lucien Stolze June 2023)

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE           :: temp_region ! region ID allocated to each cell
  REAL(DP), DIMENSION(:,:), ALLOCATABLE             :: temp_ts ! temperature of the time series
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: t_temp_ts ! time of the time series
  INTEGER(I4B)                                      :: nb_temp_ts ! nb of time series
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: reg_temp_ts ! region ID allocated to each time series
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: temp_fix ! temperature fix
  REAL(DP), DIMENSION(:), ALLOCATABLE               :: reg_temp_fix ! region ID allocated to each temp fix
  INTEGER(I4B)                                      :: nbreg ! region ID allocated to each temp fix
  INTEGER(I4B)                                      :: nb_temp_fix ! nb of temp fix
  LOGICAL(LGT)                                      :: RunTempts ! boolean to activate temperature from time series

END MODULE temperature
