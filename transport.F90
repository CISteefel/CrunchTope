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
    
MODULE transport
  
  USE crunchtype
  USE params

  LOGICAL(LGT)                                     :: species_diffusion
  LOGICAL(LGT)                                     :: xflow
  LOGICAL(LGT)                                     :: yflow
  LOGICAL(LGT)                                     :: zflow
  LOGICAL(LGT)                                     :: UseThresholdPorosity
  LOGICAL(LGT)                                     :: MillingtonQuirk
  LOGICAL(LGT)                                     :: Constant_Tortuosity
  CHARACTER (LEN=mls)                              :: TortuosityOption

  REAL(DP)                                         :: alfL
  REAL(DP)                                         :: alfT
  REAL(DP)                                         :: dcoeff
  REAL(DP)                                         :: dzero
  REAL(DP)                                         :: activation
  REAL(DP)                                         :: formation
  REAL(DP)                                         :: uli
  REAL(DP)                                         :: dgas
  REAL(DP)                                         :: d_25
  REAL(DP)                                         :: sumsigma_w
  REAL(DP)                                         :: sumsigma_e
  REAL(DP)                                         :: sumsigma_s
  REAL(DP)                                         :: sumsigma_n
  REAL(DP)                                         :: dgradw
  REAL(DP)                                         :: dgrade
  REAL(DP)                                         :: dgrads
  REAL(DP)                                         :: dgradn
  REAL(DP)                                         :: anisotropyY
  REAL(DP)                                         :: anisotropyZ
  REAL(DP)                                         :: ThresholdPorosity
  REAL(DP)                                         :: TortuosityBelowThreshold
  REAL(DP)                                         :: TortuosityAboveThreshold

  REAL(DP), DIMENSION(16)                                         :: TempFlux

  INTEGER(I4B)                                     :: idiffus
  INTEGER(I4B)                                     :: nTortuosityZone
  INTEGER(I4B)                                     :: MeanDiffusion

  REAL(DP), DIMENSION(:), ALLOCATABLE              :: d_sp
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: d_correct
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: sigma_w
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: sigma_e
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: sigma_s
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: sigma_n
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: cflowin
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: cflowout
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: FluidBuryX
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: FluidBuryY
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: SolidBuryX
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: SolidBuryY
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: aabur
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: bbbur
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: ccbur
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: rrbur
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: uubur
  REAL(DP), DIMENSION(:), ALLOCATABLE                           :: TortuosityZone

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: jxxTortuosity_lo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: jyyTortuosity_lo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: jzzTortuosity_lo
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: jxxTortuosity_hi
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: jyyTortuosity_hi
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                       :: jzzTortuosity_hi

! Allocatable arrays dimensioned over spatial domain

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: satliq
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: satliqold

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: a
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: b
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: c
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: d
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: e
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: f
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: ag
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: bg
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: cg
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dg
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: eg
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: fg
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: abu
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: bbu
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: cbu
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dbu
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: ebu
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: fbu
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: a_d
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: b_d
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: c_d
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: d_d
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: f_d
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: e_d
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: aDD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: bDD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: cDD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dDD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: fDD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: eDD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: gDD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: iDD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: hDD
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: ftvd
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: gtvd
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: htvd
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: sorp
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: ctvd
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: netflowx
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: netDiffuseX
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dstar
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: s_dsp
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: s_chg
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: sumwtchg
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: qx
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: qy 
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: qz

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: qxgas
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: qygas
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: qzgas

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dspx
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dspy
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dspz
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: qg
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: tortuosity

  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: fluxx
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE        :: fluxy
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: dfluxx
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: dfluxy
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: dflux_x
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: advflux_x
  INTEGER(I4B), DIMENSION(:,:,:), ALLOCATABLE      :: ActiveCell
  INTEGER(I4B), DIMENSION(:,:,:,:), ALLOCATABLE    :: jxe
  INTEGER(I4B), DIMENSION(:,:,:,:), ALLOCATABLE    :: jxw
  INTEGER(I4B), DIMENSION(:,:,:,:), ALLOCATABLE    :: jyn
  INTEGER(I4B), DIMENSION(:,:,:,:), ALLOCATABLE    :: jys

! Zhi Li
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: us
REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: vs

END module transport
