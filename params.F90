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
    
MODULE params

  USE crunchtype

  INTEGER(I4B), PARAMETER :: mcomp=50
  INTEGER(I4B), PARAMETER :: mrct=500
  INTEGER(I4B), PARAMETER :: mspec=1000
  INTEGER(I4B), PARAMETER :: mgas=75
  INTEGER(I4B), PARAMETER :: msurf=25
  INTEGER(I4B), PARAMETER :: mtot=mcomp+mspec
  INTEGER(I4B), PARAMETER :: mexch=15
  INTEGER(I4B), PARAMETER :: meqn=mcomp+mexch+msurf
  INTEGER(I4B), PARAMETER :: mexch_sec=150
  INTEGER(I4B), PARAMETER :: msurf_sec=250
  INTEGER(I4B), PARAMETER :: msurftot=msurf+msurf_sec
  INTEGER(I4B), PARAMETER :: mx=1000
  INTEGER(I4B), PARAMETER :: mp=meqn*mx*meqn*3
  INTEGER(I4B), PARAMETER :: my=1000
  INTEGER(I4B), PARAMETER :: mz=1
  INTEGER(I4B), PARAMETER :: mwidth=2
  INTEGER(I4B), PARAMETER :: mxyz=mx*my*mz
  INTEGER(I4B), PARAMETER :: nc=50
  INTEGER(I4B), PARAMETER :: meqneq=nc+mexch+msurf
  INTEGER(I4B), PARAMETER :: nm=mrct
  INTEGER(I4B), PARAMETER :: mcmplx=mspec
  INTEGER(I4B), PARAMETER :: ng=mgas
  INTEGER(I4B), PARAMETER :: msec=mcmplx+nm+ng+msurf_sec
  INTEGER(I4B), PARAMETER :: ndim=nc+mcmplx+nm+ng+msurf_sec
  INTEGER(I4B), PARAMETER :: nleq=4  
  INTEGER(I4B), PARAMETER :: ntmp=8 
  INTEGER(I4B), PARAMETER :: mbasis=1
  INTEGER(I4B), PARAMETER :: mls=132 
  INTEGER(I4B), PARAMETER :: mchem=10000
  INTEGER(I4B), PARAMETER :: mhet=10000
  INTEGER(I4B), PARAMETER :: mperm=10000
  INTEGER(I4B), PARAMETER :: mzone=60
  INTEGER(I4B), PARAMETER :: mreact=10
  INTEGER(I4B), PARAMETER :: mpre=35
  INTEGER(I4B), PARAMETER :: maqkin=25

  REAL(DP), PARAMETER :: rgasKCAL=0.001987d0      !!  kcal/deg-mol
  REAL(DP), PARAMETER :: rgas=.00831470       !!  kJ/deg-mol
  REAL(DP), PARAMETER :: tk25=1.0_dp/298.15_dp 
  REAL(DP), PARAMETER :: clg=2.30258509299405
  REAL(DP), PARAMETER :: secyr=365.0d0*24.0d0*60.0d0*60.0d0
  REAL(DP), PARAMETER :: ed=5.0
  REAL(DP), PARAMETER :: d0=1.e-09 
  REAL(DP), PARAMETER :: t_half = -0.6931471805599453

END MODULE params

