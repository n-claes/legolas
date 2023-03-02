! =============================================================================
!> This submodule defines internal kink modes in force-free magnetic fields.
!! The geometry is cylindrical with parabolic density and velocity profiles,
!! The geometry can be overridden in the parfile.
!!
!! This equilibrium is taken from section III.B in
!! _Goedbloed, J. P. "The Spectral Web of stationary plasma equilibria.
!! II. Internal modes." Physics of Plasmas 25.3 (2018): 032110_.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 1
!! - <tt>k3</tt> = \( 0.16\alpha \)
!! - <tt>cte_rho0</tt> = 1 : used as prefactor in setting the density.
!! - <tt>cte_v03</tt> = 1 : used as prefactor in setting the z-component of velocity.
!! - <tt>cte_p0</tt> = 3 : used to set the pressure.
!! - <tt>alpha</tt> = 5 / x_end : used in the Bessel functions.
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_internal_kink_instability
  implicit none

contains

  module procedure internal_kink_eq
    use mod_equilibrium_params, only: cte_rho0, cte_v03, cte_p0, alpha

    real(dp)      :: r, x, a0
    real(dp)      :: J0, J1, J2, DJ0, DJ1
    integer       :: i

    call settings%grid%set_geometry("cylindrical")
    call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
    call initialise_grid(settings)

    a0 = settings%grid%get_grid_end()
    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%physics%enable_flow()
      cte_rho0 = 1.0d0
      cte_v03  = 1.0d0
      cte_p0 = 9.0d0
      alpha = 5.0d0 / a0

      k2 = 1.0d0
      k3 = 0.16d0 * alpha
    end if ! LCOV_EXCL_STOP

    do i = 1, settings%grid%get_gauss_gridpts()
      r = grid_gauss(i)
      x = r / a0

      J0   = bessel_jn(0, alpha * x)
      J1   = bessel_jn(1, alpha * x)
      J2   = bessel_jn(2, alpha * x)
      DJ0  = -alpha * J1
      DJ1  = alpha * (0.5d0 * J0 - 0.5d0 * J2)

      rho_field % rho0(i) = cte_rho0 * (1.0d0 - x**2 / a0**2)
      v_field % v03(i)    = cte_v03 * (1.0d0 - x**2 / a0**2)
      B_field % B02(i)    = J1
      B_field % B03(i)    = J0
      B_field % B0(i)     = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i)     = cte_p0 / (rho_field % rho0(i))

      rho_field % d_rho0_dr(i) = -2.0d0 * cte_rho0 * x / a0
      T_field % d_T0_dr(i)     = 2.0d0 * x * cte_p0 / ( &
        a0**2 * cte_rho0 * (1.0d0 - x**2)**2 &
      )
      v_field % d_v03_dr(i)    = -2.0d0 * cte_v03 * x
      B_field % d_B02_dr(i)    = DJ1
      B_field % d_B03_dr(i)    = DJ0
    end do
  end procedure internal_kink_eq

end submodule smod_equil_internal_kink_instability
