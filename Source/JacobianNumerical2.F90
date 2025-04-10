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

SUBROUTINE JacobianNumerical2(ncomp,nspec,nx,ny,nz)
USE crunchtype
USE params
USE concentration, only: s,sp10,sp,gam,keqaq,muaq,NumericalJac,sppTMP,sppTMP10,gamTMP,keqaqTMP
USE solver, only: fjac
USE temperature, only: T


IMPLICIT NONE

!  External variables

INTEGER(I4B), INTENT(IN)                             :: ncomp
INTEGER(I4B), INTENT(IN)                             :: nspec
INTEGER(I4B), INTENT(IN)                             :: nx
INTEGER(I4B), INTENT(IN)                             :: ny
INTEGER(I4B), INTENT(IN)                             :: nz

!  Internal variables

REAL(DP)                                             :: sum
REAL(DP)                                             :: spec_conc
REAL(DP)                                             :: mutemp

REAL(DP)                                             :: tempc

INTEGER(I4B)                                         :: i
INTEGER(I4B)                                         :: ik
INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: kk
INTEGER(I4B)                                         :: ksp
INTEGER(I4B)                                         :: jx
INTEGER(I4B)                                         :: jy
INTEGER(I4B)                                         :: jz


NumericalJac = 0.0d0

do jz = 1,nz
  do jy = 1,ny
    do jx = 1,nx   
      
      tempc = t(jx,jy,jz) 
        
        do i = 1,ncomp+nspec
          sppTMP10(i) = sp10(i,jx,jy,jz)
          sppTMP(i) = sp(i,jx,jy,jz)
          gamTMP(i) = gam(i,jx,jy,jz)
        end do
        do ksp = 1,nspec
          keqaqTMP(ksp) = keqaq(ksp,jx,jy,jz)
        end do        

      call GammaResidual(ncomp,nspec,tempc,jx,jy,jz)
                   
      do i = 1,ncomp+nspec
!!!         sp10(i,jx,jy,jz) = sppTMP10(i)
!!!         sp(i,jx,jy,jz) = sppTMP(i)
         gam(i,jx,jy,jz) = gamTMP(i)
      end do
      
      do ik = 1,ncomp
        do i = 1,ncomp
!!!          write(*,*) NumericalJac(ik,i),fjac(ik,i,jx,jy,jz)
          fjac(ik,i,jx,jy,jz) = NumericalJac(ik,i)
        end do
      end do

  
    END DO
  END DO
END DO


RETURN
END SUBROUTINE JacobianNumerical2
!***********************************************************
