! =============================================================================
!> This submodule defines an unperturbed magnetised jet model in cylindrical
!! geometry, giving rise to Kelvin-Helmholtz and current-driven instabilities.
!! The geometry is fixed for this problem; the cylinder wall is dependent
!! on the equilibrium parameters and is given by <tt>2rj</tt>.
!!
!! This equilibrium is taken from
!! _Baty, H., & Keppens, R. (2002). Interplay between Kelvin-Helmholtz and
!!  current-driven instabilities in jets. The Astrophysical Journal, 580(2), 800_.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = -1
!! - <tt>V</tt> = 1.63 : amplitude of the velocity shear
!! - <tt>cte_p0</tt> = 1 : used to set the pressure.
!! - <tt>cte_rho0</tt> = 1 : used to set the density.
!! - <tt>Bz0</tt> = 0.25 : used to set Bz.
!! - <tt>rc</tt> = 0.5 : length for radial variation.
!! - <tt>rj</tt> = 1 : jet radius
!!
!! and can all be changed in the parfile. @endnote
!! @note The default setup is _HEL2_ in the original paper.
!!       For _HEL1_ you can set <tt>rc = 2</tt>. @endnote
submodule (mod_equilibrium) smod_equil_kelvin_helmholtz_cd
  use mod_equilibrium_params, only: V, cte_rho0, cte_p0, Bz0, rc, Bth0, rj
  implicit none

  real(dp) :: a

contains

  !> Sets the equilibrium.
  module procedure kh_cd_instability_eq
    real(dp)    :: r
    integer     :: i

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%physics%enable_flow()

      V = 1.63_dp
      cte_rho0 = 1.0_dp
      cte_p0 = 1.0_dp
      Bz0 = 0.25_dp
      rj = 1.0_dp

      rc    = 0.5_dp
      k2  = -1.0_dp
    end if ! LCOV_EXCL_STOP

    Bth0  = 0.4_dp * (rc**2+rj**2) / (rj*rc)
    a = 0.1_dp * rj
    k3  = dpi / rj

    call settings%grid%set_geometry("cylindrical")
    call settings%grid%set_grid_boundaries(0.0_dp, 2.0_dp * rj)
    call initialise_grid(settings)

    do i = 1, settings%grid%get_gauss_gridpts()
      r = grid_gauss(i)

      rho_field % rho0(i) = rho0()
      v_field % v03(i) = v03(r)
      B_field % B02(i) = B02(r)
      B_field % B03(i) = B03()
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
      T_field % T0(i) = T0(r)
      B_field % d_B02_dr(i) = dB02(r)
      v_field % d_v03_dr(i) = dv03(r)
      T_field % d_T0_dr(i)  = dT0(r)
    end do
  end procedure kh_cd_instability_eq


  real(dp) function rho0()
    rho0 = cte_rho0
  end function rho0

  real(dp) function T0(r)
    real(dp), intent(in) :: r
    T0 = cte_p0 / rho0() &
      - (Bth0**2 / (2.0_dp * rho0())) * (1.0_dp - rc**4 / (rc**2 + r**2)**2)
  end function T0

  real(dp) function dT0(r)
    real(dp), intent(in) :: r
    dT0 = -(2.0_dp * Bth0**2 / rho0()) * rc**4 * r / (r**2 + rc**2)**3
  end function dT0

  real(dp) function v03(r)
    real(dp), intent(in) :: r
    v03 = (V / 2.0_dp) * tanh((rj - r) / a)
  end function v03

  real(dp) function dv03(r)
    real(dp), intent(in) :: r
    dv03 = -(V / (2.0_dp * a)) / cosh((rj - r) / a)**2
  end function dv03

  real(dp) function B02(r)
    real(dp), intent(in) :: r
    B02 = Bth0 * r * rc / (rc**2 + r**2 )
  end function B02

  real(dp) function dB02(r)
    real(dp), intent(in) :: r
    dB02 = Bth0 * rc * (rc**2 - r**2) / (r**2 + rc**2)**2
  end function dB02

  real(dp) function B03()
    B03 = Bz0
  end function B03

end submodule smod_equil_kelvin_helmholtz_cd
