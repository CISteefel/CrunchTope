SUBROUTINE interp(x_in, x, y, y_out, n)
USE crunchtype

implicit none
    
!  External variables and arrays
real(DP),intent(in) :: x_in
INTEGER(I4B),intent(in) :: n
real(DP),intent(in) :: y(n)
real(DP),intent(in) :: x(n)
real(DP),intent(out) :: y_out
integer(I4B) :: i

! Check if x_in is out of bounds
if (x_in <= x(1)) then
  y_out = y(1)
elseif (x_in >= x(n)) then
  y_out = y(n)
else
  ! Find the interval [x(i), x(i+1)] such that x(i) <= x_in < x(i+1)
  do i = 1, n-1
    if (x_in >= x(i) .and. x_in <= x(i+1)) then
      ! Perform linear interpolation
      y_out = y(i) + (y(i+1) - y(i)) * (x_in - x(i)) / (x(i+1) - x(i))
      exit
    end if
  end do
end if
END SUBROUTINE interp