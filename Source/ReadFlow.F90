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
    
MODULE ReadFlow

  USE crunchtype
  USE params

  CHARACTER (LEN=mls)                                        :: NuftFile
  CHARACTER (LEN=mls), DIMENSION(:), ALLOCATABLE             :: NuftVariable

  LOGICAL(LGT)                                               :: NuftRead
  LOGICAL(LGT)                                               :: ReadNuft
  LOGICAL(LGT)                                               :: AlreadyUnsaturated

  REAL(DP)                                                   :: timeNuft
  REAL(DP)                                                   :: timeNuftOld

  INTEGER(I4B)                                               :: nxNuft
  INTEGER(I4B)                                               :: nyNuft
  INTEGER(I4B)                                               :: nzNuft
  INTEGER(I4B)                                               :: iunitNuft
  INTEGER(I4B)                                               :: NumNuftVariables
  INTEGER(I4B)                                               :: NumSourceTerms
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                    :: jxNuftSource
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                    :: jyNuftSource
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE                    :: jzNuftSource
  INTEGER(I4B)                                               :: NumNuftSteps
  INTEGER(I4B)                                               :: TotNuftSteps
  REAL(DP)                                                   :: NuftWeightNew
  REAL(DP)                                                   :: NuftWeightOld

  LOGICAL(LGT)                                               :: ReadNuftLiqFlux
  LOGICAL(LGT)                                               :: ReadNuftSaturation
  LOGICAL(LGT)                                               :: ReadNuftPorosity
  LOGICAL(LGT)                                               :: ReadNuftTemperature
  LOGICAL(LGT)                                               :: ReadNuftLiqDensity
  LOGICAL(LGT)                                               :: ReadNuftGasFlux
  LOGICAL(LGT)                                               :: ReadNuftGasDensity

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qxNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qyNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qzNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qxNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qyNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qzNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qxgasNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qygasNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qzgasNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qxgasNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qygasNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qzgasNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: satliqNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: satliqNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qgNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: qgNuftOld
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: roNuft
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE                    :: roNuftOld


END MODULE ReadFlow