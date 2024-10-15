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

SUBROUTINE jacobian(ncomp,nspec,nx,ny,nz)
USE crunchtype
USE params
USE concentration
USE solver
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

INTEGER(I4B)                                         :: i
INTEGER(I4B)                                         :: ik
INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: kk
INTEGER(I4B)                                         :: ksp
INTEGER(I4B)                                         :: jx
INTEGER(I4B)                                         :: jy
INTEGER(I4B)                                         :: jz

REAL(DP)                                             :: perturb=1.0E-09

REAL(DP)                                                    :: tempc
REAL(DP), ALLOCATABLE, DIMENSION(:)                         :: spppTMP
REAL(DP), ALLOCATABLE, DIMENSION(:)                         :: spppTMP10
REAL(DP), ALLOCATABLE, DIMENSION(:)                         :: keqaqTMP
REAL(DP), ALLOCATABLE, DIMENSION(:)                         :: gammaTMP
REAL(DP), ALLOCATABLE, DIMENSION(:)                         :: sTotTMP
REAL(DP), ALLOCATABLE, DIMENSION(:,:)                       :: NumericalJac 
;
ALLOCATE(spppTMP(ncomp+nspec))
ALLOCATE(spppTMP10(ncomp+nspec))
ALLOCATE(keqaqTMP(nspec))
ALLOCATE(gammaTMP(ncomp+nspec))
ALLOCATE(sTotTMP(ncomp))
ALLOCATE(NumericalJac(ncomp,ncomp))

fjac = 0.0d0

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      DO ksp = 1,nspec
        
        spec_conc = sp10(ksp+ncomp,jx,jy,jz)
        DO i = 1,ncomp
          IF (muaq(ksp,i) /= 0.0) THEN
            mutemp = muaq(ksp,i)
!            DO i2 = 1,ncomp

            DO i2 = 1,i-1                                         !!! Lower triangle of matrix

              fjac(i2,i,jx,jy,jz) = fjac(i2,i,jx,jy,jz) +           &
                 mutemp*muaq(ksp,i2)*spec_conc

            END DO

            fjac(i,i,jx,jy,jz) = fjac(i,i,jx,jy,jz) +               &
                 mutemp*mutemp*spec_conc 
          END IF
        END DO
        
      END DO

      DO i = 1,ncomp
        DO i2 = 1,i-1
          fjac(i,i2,jx,jy,jz) = fjac(i2,i,jx,jy,jz)
        END DO
        fjac(i,i,jx,jy,jz) = fjac(i,i,jx,jy,jz) + sp10(i,jx,jy,jz)
      END DO
      
      tempc = t(jx,jy,jz) 
      
      do ik = 1,ncomp+nspec
        
        spppTMP(ik)   = sp(ik,jx,jy,jz)
              spppTMP10(ik) = sp10(ik,jx,jy,jz)
              gammaTMP(ik)  = gam(ik,jx,jy,jz)
      end do
      do ksp = 1,nspec
        keqaqTMP(ksp)  = keqaq(ksp,jx,jy,jz)
      end do

      
      
      NumericalJac = 0.0d0
      do ik = 1,ncomp

         spppTMP(ik) = spppTMP(ik) + perturb
         spppTMP10(ik) = DEXP(spppTMP(ik))
        
!!!        call GammaResidual(ncomp,nspec,tempc,spppTMP,keqaqTMP,gammaTMP,sTotTMP)
        
        DO ksp = 1,nspec   
          sum = 0.0D0
          DO i = 1,ncomp
            sum = sum + muaq(ksp,i)*(spppTMP(i) + gammaTMP(i))
          END DO
          spppTMP(ksp+ncomp) = keqaqTMP(ksp) - gammaTMP(ksp+ncomp) + sum
          spppTMP10(ksp+ncomp) = DEXP(spppTMP(ksp+ncomp))           
        END DO
      

        DO i = 1,ncomp
          sum=0.0D0
          DO ksp = 1,nspec
            kk = ksp + ncomp
            sum = sum + muaq(ksp,i)*spppTMP10(kk)
          END DO
          stotTMP(i) = sum + spppTMP10(i)
        END DO
       
        do i = 1,ncomp
          NumericalJac(ik,i) = ( sTotTMP(i) - s(i,jx,jy,jz) )/perturb
        end do
        
         spppTMP(ik) = spppTMP(ik) - perturb
         spppTMP10(ik) = DEXP(spppTMP(ik))
         
         do i = 1,ncomp
           fjac(ik,i,jx,jy,jz) = NumericalJac(ik,i)
          end do
        
      END DO
      

      
    END DO
  END DO
END DO


RETURN
END SUBROUTINE jacobian
!***********************************************************
