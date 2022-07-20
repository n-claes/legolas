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
  implicit none

contains

  module subroutine isothermal_atmosphere_eq()
    use mod_equilibrium_params, only: cte_rho0, cte_B02, cte_B03, cte_T0, g

    real(dp)  :: x, scale_height
    integer   :: i

    geometry = "Cartesian"
    call allow_geometry_override(default_x_start=0.0d0, default_x_end=15.0d0)
    call initialise_grid()

    if (use_defaults) then  ! LCOV_EXCL_START
      external_gravity = .true.
      cte_rho0 = 1.0d0
      cte_B02 = 0.25d0
      cte_B03 = 0.25d0
      cte_T0 = 1.0d0
      g = 5.0d0

      k2 = 0.0d0
      k3 = 2.0d0
    end if  ! LCOV_EXCL_STOP

    scale_height = cte_T0 / g

    T_field % T0 = cte_T0
    B_field % B02 = cte_B02
    B_field % B03 = cte_B03
    B_field % B0 = sqrt(cte_B02**2 + cte_B03**2)
    grav_field % grav = g

    do i = 1, gauss_gridpts
      x = grid_gauss(i)
      rho_field % rho0(i) = cte_rho0 * exp(-x / scale_height)
      rho_field % d_rho0_dr(i) = -cte_rho0 * exp(-x / scale_height) / scale_height
    end do

  end subroutine isothermal_atmosphere_eq
end submodule smod_equil_isothermal_atmosphere
