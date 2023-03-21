! =============================================================================
!> This submodule defines an equilibrium in Cartesian geometry with a
!! stratified equilibrium profile, representing a solar magnetic atmosphere.
!! The equilibrium is isothermal with a constant magnetic field.
!! The scale height is given by
!! $$ H = \frac{\text{cte_T0}}{g}  $$
!! The geometry is fixed to Cartesian, boundaries can be overridden using the parfile.
!!
!! This equilibrium is taken from
!! _Nye, A., & Thomas, J. (1976). Solar magneto-atmospheric waves. I. An exact solution
!! for a horizontal magnetic field. The Astrophysical Journal, 204._
!! [_573-581_](http://articles.adsabs.harvard.edu/pdf/1976ApJ...204..573N).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 2
!! - <tt>cte_rho0</tt> = 1 : used as a density prefactor.
!! - <tt>cte_B02</tt> = 0.25 : used to set the By-component
!! - <tt>cte_B03</tt> = 0.25 : used to set the Bz-component
!! - <tt>cte_T0</tt> = 1 : used to set the temperature (isothermal).
!! - <tt>g</tt> = 5 : used to set the gravity constant.
!!
!! and can all be changed in the parfile. @endnote
!! @warning This equilibrium has no regression test yet! @endwarning
submodule (mod_equilibrium) smod_equil_isothermal_atmosphere
  use mod_equilibrium_params, only: cte_rho0, cte_B02, cte_B03, cte_T0, g
  implicit none

  real(dp) :: scale_height

contains

  module procedure isothermal_atmosphere_eq
    real(dp)  :: x
    integer   :: i

    if (settings%equilibrium%use_defaults) then  ! LCOV_EXCL_START
      call settings%grid%set_geometry("Cartesian")
      call settings%grid%set_grid_boundaries(0.0_dp, 15.0_dp)
      call settings%physics%enable_gravity()
      cte_rho0 = 1.0_dp
      cte_B02 = 0.25_dp
      cte_B03 = 0.25_dp
      cte_T0 = 1.0_dp
      g = 5.0_dp

      k2 = 0.0_dp
      k3 = 2.0_dp
    end if  ! LCOV_EXCL_STOP
    call initialise_grid(settings)

    scale_height = cte_T0 / g
    grav_field % grav = g

    do i = 1, settings%grid%get_gauss_gridpts()
      x = grid_gauss(i)
      rho_field % rho0(i) = rho0(x)
      rho_field % d_rho0_dr(i) = drho0(x)
      T_field % T0 = T0()
      B_field % B02 = B02()
      B_field % B03 = B03()
      B_field % B0 = sqrt(B02()**2 + B03()**2)
    end do

  end procedure isothermal_atmosphere_eq


  real(dp) function rho0(x)
    real(dp), intent(in) :: x
    rho0 = cte_rho0 * exp(-x / scale_height)
  end function rho0

  real(dp) function drho0(x)
    real(dp), intent(in) :: x
    drho0 = -cte_rho0 * exp(-x / scale_height) / scale_height
  end function drho0

  real(dp) function T0()
    T0 = cte_T0
  end function T0

  real(dp) function B02()
    B02 = cte_B02
  end function B02

  real(dp) function B03()
    B03 = cte_B03
  end function B03

end submodule smod_equil_isothermal_atmosphere
