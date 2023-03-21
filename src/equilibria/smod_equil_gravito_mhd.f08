! =============================================================================
!> This submodule defines an equilibrium in Cartesian geometry with a
!! stratified equilibrium profile, giving rise to gravito-MHD waves.
!! The geometry can be overridden using the parfile.
!!
!! This equilibrium is taken from section 7.3.3, p. 258 in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!!  and Astrophysical Plasmas. Cambridge University Press._ [DOI](http://doi.org/10.1017/9781316403679).
!! @note Default values are given by
!!
!! - <tt>k2</tt> = \(\pi\)
!! - <tt>k3</tt> = \(\pi\)
!! - <tt>cte_p0</tt> = 0.5 : used to set the pressure value.
!! - <tt>alpha</tt> = 20 : used to constrain the density.
!! - <tt>g</tt> = 0.5 : used to set the gravity constant.
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_gravito_mhd
  use mod_equilibrium_params, only: g, cte_rho0, cte_p0, alpha, beta
  implicit none

  real(dp) :: B0

contains

  module procedure gravito_mhd_eq
    real(dp)  :: x
    integer   :: i

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("Cartesian")
      call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
      call settings%physics%enable_gravity()

      k2 = dpi
      k3 = dpi
      cte_p0 = 0.5_dp
      g = 0.5_dp
      alpha = 20.0_dp
    end if ! LCOV_EXCL_STOP

    B0 = 1.0_dp
    beta  = 2.0_dp * cte_p0 / B0**2
    cte_rho0 = (alpha / g) * (cte_p0 + 0.5_dp * B0**2)

    call initialise_grid(settings)

    !! Equilibrium
    T_field % T0 = cte_p0 / cte_rho0
    grav_field % grav = g

    ! Full exponential prescription for the equilibrium configuration
    do i = 1, settings%grid%get_gauss_gridpts()
      x = grid_gauss(i)

      rho_field % rho0(i) = rho0(x)
      rho_field % d_rho0_dr(i) = drho0(x)
      B_field % B03(i) = B03(x)
      B_field % d_B03_dr(i) = dB03(x)
      B_field % B0(i) = sqrt((B_field % B02(i))**2 + (B_field % B03(i))**2)
    end do

  end procedure gravito_mhd_eq


  real(dp) function rho0(x)
    real(dp), intent(in) :: x
    rho0 = cte_rho0 * exp(-alpha * x)
  end function rho0

  real(dp) function drho0(x)
    real(dp), intent(in) :: x
    drho0 = -alpha * rho0(x)
  end function drho0

  real(dp) function T0()
    T0 = cte_p0 / cte_rho0
  end function T0

  real(dp) function B03(x)
    real(dp), intent(in) :: x
    B03 = B0 * exp(-0.5_dp * alpha * x)
  end function B03

  real(dp) function dB03(x)
    real(dp), intent(in) :: x
    dB03 = -0.5_dp * alpha * B03(x)
  end function dB03

end submodule smod_equil_gravito_mhd
