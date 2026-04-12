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


MODULE species_calculations
!==============================================================================
! MODULE: species_calculations
!
! DESCRIPTION:
!   This module provides pure functions for calculating aqueous species 
!   concentrations and their derivatives in geochemical systems.
!
! CONTAINS:
!   - species_result: Derived type for storing calculation results
!   - calculate_species: Pure function for species concentration calculations
!   - deallocate_species_result: Utility subroutine for memory cleanup
!
! AUTHOR: Converted from original SUBROUTINE species
! DATE: 2025
!==============================================================================

IMPLICIT NONE

! Precision parameters
INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: DP = SELECTED_REAL_KIND(15, 307)

PRIVATE
PUBLIC :: species_result, calculate_species, deallocate_species_result

!------------------------------------------------------------------------------
! Derived type to hold species calculation results
!------------------------------------------------------------------------------
TYPE :: species_result
  REAL(DP), ALLOCATABLE :: sp(:,:,:,:)           ! Species concentrations (log scale)
  REAL(DP), ALLOCATABLE :: sp10(:,:,:,:)         ! Species concentrations (linear scale)
  REAL(DP), ALLOCATABLE :: deriv_conc(:,:,:,:,:) ! Concentration derivatives
  LOGICAL :: allocated = .FALSE.                 ! Allocation status flag
END TYPE species_result

CONTAINS

!==============================================================================
PURE FUNCTION calculate_species(ncomp, nspec, nsurf, nexchange, npot, nx, ny, nz, &
                               sp_in, muaq, ulab, lngamma, lngammawater, keqaq, &
                               gamma, deriv_gamma) RESULT(result_data)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Calculates aqueous species concentrations and their derivatives for 
!   geochemical equilibrium calculations.
!
! INPUTS:
!   ncomp       - Number of components
!   nspec       - Number of secondary species  
!   nsurf       - Number of surface species
!   nexchange   - Number of exchange species
!   npot        - Number of potential species
!   nx,ny,nz    - Grid dimensions
!   sp_in       - Initial species concentrations
!   muaq        - Stoichiometric coefficients matrix
!   ulab        - Species labels
!   lngamma     - Natural log of activity coefficients
!   lngammawater- Natural log of water activity coefficient
!   keqaq       - Equilibrium constants
!   gamma       - Activity coefficients
!   deriv_gamma - Derivatives of activity coefficients
!
! OUTPUTS:
!   result_data - Structure containing calculated concentrations and derivatives
!------------------------------------------------------------------------------

IMPLICIT NONE

! Input parameters
INTEGER(I4B), INTENT(IN) :: ncomp, nspec, nsurf, nexchange, npot, nx, ny, nz

! Input arrays
REAL(DP), INTENT(IN) :: sp_in(:,:,:,:)              ! Initial species concentrations
REAL(DP), INTENT(IN) :: muaq(:,:)                   ! Stoichiometric coefficients
CHARACTER(LEN=*), INTENT(IN) :: ulab(:)             ! Species labels
REAL(DP), INTENT(IN) :: lngamma(:,:,:,:)            ! Natural log of activity coefficients
REAL(DP), INTENT(IN) :: lngammawater(:,:,:)         ! Natural log of water activity coefficient
REAL(DP), INTENT(IN) :: keqaq(:,:,:,:)              ! Equilibrium constants
REAL(DP), INTENT(IN) :: gamma(:,:,:,:)              ! Activity coefficients
REAL(DP), INTENT(IN) :: deriv_gamma(:,:,:,:,:)      ! Derivatives of activity coefficients

! Result
TYPE(species_result) :: result_data

! Local variables
INTEGER(I4B) :: ksp, icomp, jx, jy, jz, ik, pos_der
CHARACTER(LEN=3) :: ulabPrint
REAL(DP) :: sum, sumderiv_gamma

! Allocate result arrays
ALLOCATE(result_data%sp(SIZE(sp_in,1), nx, ny, nz))
ALLOCATE(result_data%sp10(SIZE(sp_in,1), nx, ny, nz))
ALLOCATE(result_data%deriv_conc(SIZE(sp_in,1), SIZE(deriv_gamma,2), nx, ny, nz))

! Set allocation flag
result_data%allocated = .TRUE.

! Initialize result arrays
result_data%sp = sp_in
result_data%sp10 = 0.0_DP
result_data%deriv_conc = 0.0_DP

! Main calculation loops
DO jz = 1, nz
  DO jy = 1, ny
    DO jx = 1, nx

      !------------------------------------------------------------------------
      ! Calculate secondary species concentrations
      !------------------------------------------------------------------------
      DO ksp = 1, nspec
        ik = ksp + ncomp
     
        sum = 0.0_DP
        DO icomp = 1, ncomp
          ulabPrint = ulab(icomp)
          IF (ulabPrint(1:3) /= 'H2O' .AND. ulabPrint(1:3) /= 'HHO') THEN
            sum = sum + muaq(ksp,icomp) * (result_data%sp(icomp,jx,jy,jz) + lngamma(icomp,jx,jy,jz))
          ELSE
            sum = sum + muaq(ksp,icomp) * lngammawater(jx,jy,jz)
          END IF
        END DO

        result_data%sp(ik,jx,jy,jz) = keqaq(ksp,jx,jy,jz) - lngamma(ik,jx,jy,jz) + sum
        result_data%sp10(ik,jx,jy,jz) = EXP(result_data%sp(ik,jx,jy,jz))
             
      END DO

      !------------------------------------------------------------------------
      ! Calculate derivatives for primary species
      !------------------------------------------------------------------------
      DO icomp = 1, ncomp
        ulabPrint = ulab(icomp)
        IF (ulabPrint(1:3) /= 'H2O' .AND. ulabPrint(1:3) /= 'HHO') THEN
          result_data%deriv_conc(icomp,icomp,jx,jy,jz) = result_data%sp10(icomp,jx,jy,jz)    
        ELSE 
          result_data%deriv_conc(icomp,icomp,jx,jy,jz) = 1.0_DP
        END IF
      END DO

      !------------------------------------------------------------------------
      ! Calculate derivatives for secondary species
      !------------------------------------------------------------------------
      DO ksp = 1, nspec
        ik = ksp + ncomp

        ! Derivative with respect to primary species concentrations
        DO icomp = 1, ncomp   
          ulabPrint = ulab(icomp)
          IF (ulabPrint(1:3) /= 'H2O' .AND. ulabPrint(1:3) /= 'HHO') THEN
            result_data%deriv_conc(ik,icomp,jx,jy,jz) = muaq(ksp,icomp) * result_data%sp10(ik,jx,jy,jz)
          ELSE
            ! Derivative with respect to water activity
            pos_der = ncomp + nexchange + nsurf + npot + 1 + 1
            result_data%deriv_conc(ik,pos_der,jx,jy,jz) = muaq(ksp,icomp) * result_data%sp10(ik,jx,jy,jz)
          END IF
        END DO
        
        ! Derivative with respect to ionic strength
        pos_der = ncomp + nexchange + nsurf + npot + 1
        sumderiv_gamma = 0.0_DP
        DO icomp = 1, ncomp
          ulabPrint = ulab(icomp)
          IF (ulabPrint(1:3) /= 'H2O' .AND. ulabPrint(1:3) /= 'HHO') THEN
            sumderiv_gamma = sumderiv_gamma + muaq(ksp,icomp) / gamma(icomp,jx,jy,jz) * &
                            deriv_gamma(icomp,pos_der,jx,jy,jz)   
          END IF
        END DO
        result_data%deriv_conc(ik,pos_der,jx,jy,jz) = result_data%sp10(ik,jx,jy,jz) * &
                                                      (sumderiv_gamma - deriv_gamma(ik,pos_der,jx,jy,jz) / &
                                                       gamma(ik,jx,jy,jz))

      END DO
      
    END DO 
  END DO
END DO

END FUNCTION calculate_species

!==============================================================================
SUBROUTINE deallocate_species_result(result_data)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   Utility subroutine to safely deallocate memory from species_result type
!
! INPUTS/OUTPUTS:
!   result_data - Structure to deallocate
!------------------------------------------------------------------------------

IMPLICIT NONE

TYPE(species_result), INTENT(INOUT) :: result_data

IF (result_data%allocated) THEN
  IF (ALLOCATED(result_data%sp)) DEALLOCATE(result_data%sp)
  IF (ALLOCATED(result_data%sp10)) DEALLOCATE(result_data%sp10)
  IF (ALLOCATED(result_data%deriv_conc)) DEALLOCATE(result_data%deriv_conc)
  result_data%allocated = .FALSE.
END IF

END SUBROUTINE deallocate_species_result

END MODULE species_calculations