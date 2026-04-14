module chemical_speciation
  implicit none
  
  ! Define precision
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: one = 1.0_dp
  
  private
  public :: calculate_secondary_species, calculate_species_derivatives, &
           calculate_primary_derivatives, is_water_species, &
           calculate_log_activity_sum, dp

contains

  !> Check if a species is water (H2O or HHO)
  pure function is_water_species(species_label) result(is_water)
    implicit none
    character(len=*), intent(in) :: species_label
    logical :: is_water
    
    character(len=20) :: label_trim
    
    label_trim = adjustl(species_label)
    is_water = (label_trim(1:3) == 'H2O' .or. label_trim(1:3) == 'HHO')
    
  end function is_water_species

  !> Calculate the sum term for equilibrium expression
  pure function calculate_log_activity_sum(ncomp, ksp, jx, jy, jz, &
                                         muaq, sp, lngamma, lngammawater, ulab) result(sum_val)
    implicit none
    
    ! Arguments
    integer, intent(in) :: ncomp, ksp, jx, jy, jz
    real(dp), intent(in) :: muaq(:,:)
    real(dp), intent(in) :: sp(:,:,:,:)
    real(dp), intent(in) :: lngamma(:,:,:,:)
    real(dp), intent(in) :: lngammawater(:,:,:)
    character(len=*), intent(in) :: ulab(:)
    
    ! Result
    real(dp) :: sum_val
    
    ! Local variables
    integer :: icomp
    
    sum_val = zero
    
    do icomp = 1, ncomp
      if (.not. is_water_species(ulab(icomp))) then
        sum_val = sum_val + muaq(ksp,icomp) * (sp(icomp,jx,jy,jz) + lngamma(icomp,jx,jy,jz))
      else
        sum_val = sum_val + muaq(ksp,icomp) * lngammawater(jx,jy,jz)
      end if
    end do
    
  end function calculate_log_activity_sum

  !> Calculate secondary species concentrations (log and linear scale)
  pure subroutine calculate_secondary_species(nspec, ncomp, jx, jy, jz, &
                                            muaq, keqaq, sp, lngamma, &
                                            lngammawater, ulab, sp_out, sp10_out)
    implicit none
    
    ! Arguments
    integer, intent(in) :: nspec, ncomp, jx, jy, jz
    real(dp), intent(in) :: muaq(:,:)
    real(dp), intent(in) :: keqaq(:,:,:,:)
    real(dp), intent(in) :: sp(:,:,:,:)
    real(dp), intent(in) :: lngamma(:,:,:,:)
    real(dp), intent(in) :: lngammawater(:,:,:)
    character(len=*), intent(in) :: ulab(:)
    real(dp), intent(out) :: sp_out(:)      ! Secondary species log concentrations
    real(dp), intent(out) :: sp10_out(:)    ! Secondary species concentrations
    
    ! Local variables
    integer :: ksp, ik
    real(dp) :: sum_val
    
    do ksp = 1, nspec
      ik = ksp + ncomp
      
      ! Calculate equilibrium log activity sum
      sum_val = calculate_log_activity_sum(ncomp, ksp, jx, jy, jz, &
                                          muaq, sp, lngamma, lngammawater, ulab)
      
      ! Calculate log concentration of secondary species
      sp_out(ksp) = keqaq(ksp,jx,jy,jz) - lngamma(ik,jx,jy,jz) + sum_val
      
      ! Calculate linear concentration
      sp10_out(ksp) = exp(sp_out(ksp))
    end do
    
  end subroutine calculate_secondary_species

  !> Calculate derivative of log activity sum with respect to primary species icomp
  pure function calculate_sum_derivative(ncomp, ksp, icomp_deriv, &
                                       muaq, ulab) result(deriv_sum)
    implicit none
    
    ! Arguments
    integer, intent(in) :: ncomp, ksp, icomp_deriv
    real(dp), intent(in) :: muaq(:,:)
    character(len=*), intent(in) :: ulab(:)
    
    ! Result
    real(dp) :: deriv_sum
    
    ! The derivative of sum with respect to ln(C_icomp_deriv) is simply muaq(ksp, icomp_deriv)
    ! if the species is not water, zero otherwise (water activity handled separately)
    if (.not. is_water_species(ulab(icomp_deriv))) then
      deriv_sum = muaq(ksp, icomp_deriv)
    else
      deriv_sum = zero
    end if
    
  end function calculate_sum_derivative

  !> Calculate derivatives of secondary species with respect to primary species
  pure subroutine calculate_species_derivatives(nspec, ncomp, jx, jy, jz, &
                                              muaq, sp10, ulab, deriv_matrix)
    implicit none
    
    ! Arguments
    integer, intent(in) :: nspec, ncomp, jx, jy, jz
    real(dp), intent(in) :: muaq(:,:)
    real(dp), intent(in) :: sp10(:,:,:,:)  ! Secondary species concentrations
    character(len=*), intent(in) :: ulab(:)
    real(dp), intent(out) :: deriv_matrix(:,:)  ! (nspec, ncomp) - d[S_k]/d[C_i]
    
    ! Local variables
    integer :: ksp, ik, icomp
    real(dp) :: deriv_sum
    
    do ksp = 1, nspec
      ik = ksp + ncomp
      
      do icomp = 1, ncomp
        if (.not. is_water_species(ulab(icomp))) then
          ! d[S_k]/d[C_i] = S_k * d(ln S_k)/d(ln C_i) = S_k * muaq(ksp, icomp)
          deriv_matrix(ksp, icomp) = muaq(ksp, icomp) * sp10(ik, jx, jy, jz)
        else
          ! Water species derivatives handled separately
          deriv_matrix(ksp, icomp) = zero
        end if
      end do
    end do
    
  end subroutine calculate_species_derivatives

  !> Calculate derivatives of primary species (identity matrix on diagonal)
  pure subroutine calculate_primary_derivatives(ncomp, jx, jy, jz, &
                                              sp10, ulab, deriv_matrix)
    implicit none
    
    ! Arguments
    integer, intent(in) :: ncomp, jx, jy, jz
    real(dp), intent(in) :: sp10(:,:,:,:)  ! Primary species concentrations
    character(len=*), intent(in) :: ulab(:)
    real(dp), intent(out) :: deriv_matrix(:,:)  ! (ncomp, ncomp) - d[C_i]/d[C_j]
    
    ! Local variables
    integer :: icomp, jcomp
    
    ! Initialize to zero
    deriv_matrix = zero
    
    do icomp = 1, ncomp
      do jcomp = 1, ncomp
        if (icomp == jcomp) then
          if (.not. is_water_species(ulab(icomp))) then
            ! d[C_i]/d[C_i] = C_i (for concentration-based derivatives)
            deriv_matrix(icomp, jcomp) = sp10(icomp, jx, jy, jz)
          else
            ! Water activity is usually fixed at 1.0
            deriv_matrix(icomp, jcomp) = one
          end if
        else
          ! d[C_i]/d[C_j] = 0 for i ? j
          deriv_matrix(icomp, jcomp) = zero
        end if
      end do
    end do
    
  end subroutine calculate_primary_derivatives

  !> Calculate derivative of secondary species with respect to water activity
  pure function calculate_water_derivative(ncomp, ksp, muaq, sp10_secondary, &
                                         jx, jy, jz, ulab) result(water_deriv)
    implicit none
    
    ! Arguments
    integer, intent(in) :: ncomp, ksp, jx, jy, jz
    real(dp), intent(in) :: muaq(:,:)
    real(dp), intent(in) :: sp10_secondary
    character(len=*), intent(in) :: ulab(:)
    
    ! Result
    real(dp) :: water_deriv
    
    ! Local variables
    integer :: icomp
    real(dp) :: water_coeff
    
    water_coeff = zero
    
    ! Sum stoichiometric coefficients for water species
    do icomp = 1, ncomp
      if (is_water_species(ulab(icomp))) then
        water_coeff = water_coeff + muaq(ksp, icomp)
      end if
    end do
    
    water_deriv = water_coeff * sp10_secondary
    
  end function calculate_water_derivative

  !> Master routine: Calculate all species concentrations and derivatives
  subroutine calculate_complete_speciation(nspec, ncomp, jx, jy, jz, &
                                         muaq, keqaq, sp, lngamma, &
                                         lngammawater, ulab, &
                                         sp_secondary, sp10_secondary, &
                                         deriv_primary, deriv_secondary)
    implicit none
    
    ! Arguments
    integer, intent(in) :: nspec, ncomp, jx, jy, jz
    real(dp), intent(in) :: muaq(:,:)
    real(dp), intent(in) :: keqaq(:,:,:,:)
    real(dp), intent(in) :: sp(:,:,:,:)
    real(dp), intent(in) :: lngamma(:,:,:,:)
    real(dp), intent(in) :: lngammawater(:,:,:)
    character(len=*), intent(in) :: ulab(:)
    real(dp), intent(out) :: sp_secondary(:)         ! Log concentrations
    real(dp), intent(out) :: sp10_secondary(:)       ! Concentrations
    real(dp), intent(out) :: deriv_primary(:,:)      ! Primary derivatives (ncomp, ncomp)
    real(dp), intent(out) :: deriv_secondary(:,:)    ! Secondary derivatives (nspec, ncomp)
    
    DO jz = 1,nz
      DO jy = 1,ny
        DO jx = 1,nx
    
    ! Calculate secondary species concentrations
    call calculate_secondary_species(nspec, ncomp, jx, jy, jz, &
                                   muaq, keqaq, sp, lngamma, &
                                   lngammawater, ulab, &
                                   sp_secondary, sp10_secondary)
    
    ! Calculate primary species derivatives
    call calculate_primary_derivatives(ncomp, jx, jy, jz, &
                                     sp, ulab, deriv_primary)
    
    ! Calculate secondary species derivatives
    call calculate_species_derivatives(nspec, ncomp, jx, jy, jz, &
                                     muaq, sp, ulab, deriv_secondary)
    
        END DO
      END DO
    END DO
    
  end subroutine calculate_complete_speciation

end module chemical_speciation