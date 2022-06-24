! =============================================================================
!> Module containing all dx profiles for the generation of non-uniform grids.
!! All profiles are imported in the <tt>get_accumulated_grid</tt> subroutine
!! using an _abstract interface_, so all new profiles should include 4 parameter
!! arguments p1, p2, p3, p4, even if they are not used in the prescription.
module mod_grid_profiles
  use mod_global_variables, only: dp
  use mod_physical_constants, only: dpi
  implicit none

  private

  public :: gaussian_dx
  public :: sine_dx

contains

  !> Gaussian dx function for use in the generation of a non-uniform grid
  function gaussian_dx(x, base, center, lowest, width) result(dx)
    real(dp) :: dx
    real(dp), intent(in) :: x, base, center, lowest, width

    dx = base - (base - lowest) * exp(-(x - center)**2 / (2.0d0 * width))
  end function gaussian_dx

  !> Sinusoidal dx function for use in the generation of a non-uniform grid
  function sine_dx(x, base, center, lowest, width) result(dx)
    real(dp) :: dx
    real(dp), intent(in) :: x, base, center, lowest, width

    dx = base + (base - lowest) * sin(dpi * (x - center) / width)
  end function sine_dx

end module mod_grid_profiles
