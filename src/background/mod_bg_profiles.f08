module mod_bg_profiles
  use mod_global_variables, only: dp
  implicit none

  private

  public :: zero

contains

  real(dp) function zero()
    zero = 0.0_dp
  end function zero

end module mod_bg_profiles
