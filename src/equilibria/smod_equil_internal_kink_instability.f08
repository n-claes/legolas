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
  use mod_equilibrium_params, only: cte_rho0, cte_v03, cte_p0, alpha
  implicit none

  real(dp) :: a0

contains

  module procedure internal_kink_eq
    real(dp)      :: r
    integer       :: i

    call settings%grid%set_geometry("cylindrical")
    call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
    call initialise_grid(settings)

    a0 = settings%grid%get_grid_end()
    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%physics%enable_flow()
      cte_rho0 = 1.0_dp
      cte_v03  = 1.0_dp
      cte_p0 = 9.0_dp
      alpha = 5.0_dp / a0

      k2 = 1.0_dp
      k3 = 0.16_dp * alpha
    end if ! LCOV_EXCL_STOP

    do i = 1, settings%grid%get_gauss_gridpts()
      r = grid_gauss(i)

      rho_field % rho0(i) = rho0(r)
      v_field % v03(i) = v03(r)
      B_field % B02(i) = B02(r)
      B_field % B03(i) = B03(r)
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i) = T0(r)

      rho_field % d_rho0_dr(i) = drho0(r)
      T_field % d_T0_dr(i) = dT0(r)
      v_field % d_v03_dr(i) = dv03(r)
      B_field % d_B02_dr(i) = dB02(r)
      B_field % d_B03_dr(i) = dB03(r)
    end do
  end procedure internal_kink_eq


  real(dp) function rho0(r)
    real(dp), intent(in) :: r
    real(dp) :: x
    x = r / a0
    rho0 = cte_rho0 * (1.0_dp - x**2 / a0**2)
  end function rho0

  real(dp) function drho0(r)
    real(dp), intent(in) :: r
    real(dp) :: x
    x = r / a0
    drho0 = -2.0_dp * cte_rho0 * x / a0
  end function drho0

  real(dp) function T0(r)
    real(dp), intent(in) :: r
    T0 = cte_p0 / rho0(r)
  end function T0

  real(dp) function dT0(r)
    real(dp), intent(in) :: r
    real(dp) :: x
    x = r / a0
    dT0 = 2.0_dp * x * cte_p0 / (a0**2 * cte_rho0 * (1.0_dp - x**2)**2)
  end function dT0

  real(dp) function v03(r)
    real(dp), intent(in) :: r
    real(dp) :: x
    x = r / a0
    v03 = cte_v03 * (1.0_dp - x**2 / a0**2)
  end function v03

  real(dp) function dv03(r)
    real(dp), intent(in) :: r
    real(dp) :: x
    x = r / a0
    dv03 = -2.0_dp * cte_v03 * x / a0
  end function dv03

  real(dp) function B02(r)
    real(dp), intent(in) :: r
    B02 = J1(r)
  end function B02

  real(dp) function dB02(r)
    real(dp), intent(in) :: r
    dB02 = dJ1(r)
  end function dB02

  real(dp) function B03(r)
    real(dp), intent(in) :: r
    B03 = J0(r)
  end function B03

  real(dp) function dB03(r)
    real(dp), intent(in) :: r
    dB03 = dJ0(r)
  end function dB03

  real(dp) function J0(r)
    real(dp), intent(in) :: r
    real(dp) :: x
    x = r / a0
    J0 = bessel_jn(0, alpha * x)
  end function J0

  real(dp) function J1(r)
    real(dp), intent(in) :: r
    real(dp) :: x
    x = r / a0
    J1 = bessel_jn(1, alpha * x)
  end function J1

  real(dp) function J2(r)
    real(dp), intent(in) :: r
    real(dp) :: x
    x = r / a0
    J2 = bessel_jn(2, alpha * x)
  end function J2

  real(dp) function dJ0(r)
    real(dp), intent(in) :: r
    real(dp) :: x
    x = r / a0
    dJ0 = -alpha * J1(r)
  end function dJ0

  real(dp) function dJ1(r)
    real(dp), intent(in) :: r
    real(dp) :: x
    x = r / a0
    dJ1 = alpha * (0.5_dp * J0(r) - 0.5_dp * J2(r))
  end function dJ1

end submodule smod_equil_internal_kink_instability
