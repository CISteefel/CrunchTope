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

SUBROUTINE jacobianCT(ncomp,nspec,nx,ny,nz)
USE crunchtype
USE params
USE concentration
USE solver

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
INTEGER(I4B)                                         :: i2
INTEGER(I4B)                                         :: kk
INTEGER(I4B)                                         :: ksp
INTEGER(I4B)                                         :: jx
INTEGER(I4B)                                         :: jy
INTEGER(I4B)                                         :: jz

INTEGER(I4B)                                   :: pos_der
INTEGER(I4B)                                   :: nk
INTEGER(I4B)                                   :: icomp

fjac = 0.0d0

DO jz = 1,nz
  DO jy = 1,ny
    DO jx = 1,nx

      DO icomp = 1,ncomp
        
 ! derivative of total concentration / icomp concentration
        ulabPrint = ulab(icomp)
        IF (ulabPrint(1:3) /= 'H2O' .and. ulabPrint(1:3) /= 'HHO') THEN
        DO ider = 1, ncomp
          ulabPrint = ulab(ider)
          IF (ulabPrint(1:3) /= 'H2O' .and. ulabPrint(1:3) /= 'HHO') THEN
            sum = 0D0
            DO ksp = 1,nspec
              nk = ksp + ncomp
              sum = sum + muaq(ksp,icomp) * deriv_conc(nk,ider,jx,jy,jz)
            END DO      
            fjac(ider,icomp,jx,jy,jz) = deriv_conc(icomp,ider,jx,jy,jz) + sum

          ELSE    !! For H2O
            
            fjac(ider,icomp,jx,jy,jz) = 0.d0
        ! deriv / gammawater
            pos_der = ncomp + nexchange + nsurf + npot + 1 + 1
            sum = 0D0
            DO ksp = 1,scalar%nspec
              nk = ksp + scalar%ncomp
              sum = sum + muaq(ksp,icomp) * deriv_conc(nk,pos_der,jx,jy,jz)
            END DO      
            
            aqspecies(ncell,icomp,pos_der,1)%deriv_totconc = sum
            
          END IF    
        END DO

    ! derivative of total concentration / ionic strength
        sum = 0D0
        pos_der = ncomp + nexchange + nsurf + npot + 1 
        DO ksp = 1,scalar%nspec
          nk = ksp + scalar%ncomp
          sum = sum + muaq(ksp,icomp) * deriv_conc(nk,pos_der,jx,jy,jz)
        END DO      
        
        aqspecies(ncell,icomp,pos_der,1)%deriv_totconc = sum
   
      END IF ! derivative is 0 for icomp = water

    END DO
  END DO
END DO

RETURN
END SUBROUTINE jacobianCT
!***********************************************************
