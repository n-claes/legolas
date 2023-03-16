module mod_bg_profiles
  use mod_global_variables, only: dp
  implicit none

  private

  public :: zero
  public :: from_profile

contains

  real(dp) function zero()
    zero = 0.0_dp
  end function zero


  function from_profile(func, values) result(array)
    procedure(real(dp)), pointer :: func
    real(dp), intent(in) :: values(:)
    real(dp) :: array(size(values))
    integer :: i

    do i = 1, size(values)
      array(i) = func(values(i))
    end do
  end function from_profile

end module mod_bg_profiles
