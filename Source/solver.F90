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

    
MODULE solver

  USE crunchtype
  USE params

  INTEGER(I4B)                                     :: level
  INTEGER(I4B)                                     :: Gimrtlevel
  REAL(DP)                                         :: GimrtRTOLKSP
  
  CHARACTER (LEN=mls)                              :: GIMRT_SolverMethod
  CHARACTER (LEN=mls)                              :: SolverMethod
  CHARACTER (LEN=mls)                              :: PCMethod
  CHARACTER (LEN=mls)                              :: GIMRT_PCMethod

  REAL(DP), DIMENSION(:,:) , ALLOCATABLE           ::blockfl
  REAL(DP), DIMENSION(:,:) , ALLOCATABLE           ::blockl
  REAL(DP), DIMENSION(:,:) , ALLOCATABLE           ::blockm
  REAL(DP), DIMENSION(:,:) , ALLOCATABLE           ::blockr
  REAL(DP), DIMENSION(:,:) , ALLOCATABLE           ::blockfr
  
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: fsurf_local
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: fch_local
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: fjac_loc
  REAL(DP), DIMENSION(:,:,:),   ALLOCATABLE        :: alf
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: aaa
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: fexch
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: fsurftmp
  REAL(DP), DIMENSION(:),   ALLOCATABLE            :: jac_sat
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: jac_pre
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: jac_preKin
  REAL(DP), DIMENSION(:),   ALLOCATABLE            :: jac_rateFactor

  INTEGER(I4B), DIMENSION(:),    ALLOCATABLE       :: indd

  REAL(DP), DIMENSION(:),     ALLOCATABLE          :: fxmax
  REAL(DP), DIMENSION(:),     ALLOCATABLE          :: bb
  REAL(DP), DIMENSION(:),     ALLOCATABLE          :: xn
  REAL(DP), DIMENSION(:),     ALLOCATABLE          :: fxx

!  Arrays for WATSOLV

  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: list
  INTEGER(I4B), DIMENSION(:), ALLOCATABLE          :: lrow
  INTEGER(I4B), DIMENSION(:,:), ALLOCATABLE        :: levptr
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: row

!  ALLOCATABLE arrays DIMENSIONed over spatial domain

!fp! block_gb(3,1);
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fjac
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fjac_D
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fjac_chg
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fgas
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fsurf
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fch
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fjpotncomp
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE      :: fjpotnsurf
!fp! darray=0;

! Allocatable arrays that are local

  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: fjpotncomp_local
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: fjpotnsurf_local
  REAL(DP), DIMENSION(:,:), ALLOCATABLE            :: fgas_local

! One DIMENSIONal arrays for block tridiagonal solver

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: aah
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: bbh
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE          :: cch
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE          :: yh

  INTEGER(I4B), DIMENSION(:,:),  ALLOCATABLE       :: indexx

! Arrays for PETSc

  REAL(DP), DIMENSION(:), ALLOCATABLE              :: XvecCrunchD
  REAL(DP), DIMENSION(:), ALLOCATABLE              :: BvecCrunchD

 
END MODULE solver
     
