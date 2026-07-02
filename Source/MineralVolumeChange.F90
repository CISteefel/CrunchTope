
module mineral_volume_change
  implicit none
  
  ! Define precision
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: tiny_val = tiny(1.0_dp)
  
  private
  public :: calculate_normalized_mineral_change, calculate_mineral_mole_change, &
           calculate_mineral_volume_change, calculate_mineral_change_array, &
           calculate_mineral_change_with_timestep, dp

contains

  !> Calculate normalized mineral volume change: (vxx - vxxOld) / volmol
  !! This gives the change in moles from current and previous mineral volumes
  pure function calculate_normalized_mineral_change(vxx, vxxOld, volmol) result(mole_change)
    implicit none
    
    ! Arguments
    real(dp), intent(in) :: vxx     ! Current mineral volume [L³]
    real(dp), intent(in) :: vxxOld  ! Previous mineral volume [L³]
    real(dp), intent(in) :: volmol  ! Molar volume of mineral [L³/mol]
    
    ! Result
    real(dp) :: mole_change  ! Change in moles [mol]
    
    ! Local variables
    real(dp) :: vol  ! Volume change
    
    ! Calculate volume change
    vol = vxx - vxxOld
    
    ! Check for very small molar volume to avoid division by near-zero
    if (abs(volmol) < tiny_val) then
      mole_change = zero
    else
      mole_change = vol / volmol
    end if
    
  end function calculate_normalized_mineral_change
  
  
end module mineral_volume_change