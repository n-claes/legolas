  ! =============================================================================
!> This submodule defines flow driven instabilities in a Cartesian geometry.
!! This equilibrium can not be called explicitly from the parfile, but rather acts
!! as a "parent setup" for the Rayleigh-Taylor and Kelvin-Helmholtz submodules which
!! use this specific kind of equilibrium but with different parameters.
!! This submodule is called within its implicit children.
!!
!! This equilibrium is taken from section 13.2, p. 486 in
!! _Goedbloed, H., Keppens, R., & Poedts, S. (2019). Magnetohydrodynamics of Laboratory
!!  and Astrophysical Plasmas. Cambridge University Press._ [DOI](http://doi.org/10.1017/9781316403679).
submodule (mod_equilibrium) smod_equil_flow_driven_instabilities
  use mod_equilibrium_params, only: g, delta, theta, p1, p2, p3, tau, p4, alpha, &
    cte_rho0, cte_p0
  implicit none

  real(dp) :: v0, v1, v2, phi0

contains

  module procedure flow_driven_instabilities_eq
    real(dp) :: x
    integer :: i

    call settings%grid%set_geometry("Cartesian")
    call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
    call initialise_grid(settings)
    call settings%physics%enable_flow()
    v0 = p1
    v1 = p2
    v2 = p3
    phi0 = p4

    do i = 1, settings%grid%get_gauss_gridpts()
      x = grid_gauss(i)
      rho_field%rho0(i) = rho0(x)
      rho_field%d_rho0_dr(i) = drho0()
      T_field%T0(i) = T0(x)
      T_field%d_T0_dr(i) = dT0(x)
      v_field%v02(i) = v02(x)
      v_field%v03(i) = v03(x)
      v_field%d_v02_dr(i) = dv02(x)
      v_field%d_v03_dr(i) = dv03(x)
      B_field%B02(i) = B02(x)
      B_field%B03(i) = B03(x)
      B_field%d_B02_dr(i) = dB02(x)
      B_field%d_B03_dr(i) = dB03(x)
      B_field%B0(i) = sqrt(B02(x)**2 + B03(x)**2)
      grav_field % grav(i) = g
    end do
  end procedure flow_driven_instabilities_eq


  real(dp) function rho0(x)
    real(dp), intent(in) :: x
    rho0 = cte_rho0 * (1.0_dp - delta*x)
  end function rho0

  real(dp) function drho0()
    drho0 = -cte_rho0 * delta
  end function drho0

  real(dp) function T0(x)
    real(dp), intent(in) :: x
    T0 = p_prof(x) / rho0(x)
  end function T0

  real(dp) function dT0(x)
    real(dp), intent(in) :: x
    dT0 = ( &
      -g * cte_rho0 * (1.0_dp - delta * x)**2 + cte_rho0 * delta * p_prof(x) &
    ) / rho0(x)**2
  end function dT0

  real(dp) function v02(x)
    real(dp), intent(in) :: x
    v02 = sin(theta) * v_prof(x)
  end function v02

  real(dp) function dv02(x)
    real(dp), intent(in) :: x
    dv02 = sin(theta) * (v1 + v2 * tau * cos(tau * (x - 0.5_dp)))
  end function dv02

  real(dp) function v03(x)
    real(dp), intent(in) :: x
    v03 = cos(theta) * v_prof(x)
  end function v03

  real(dp) function dv03(x)
    real(dp), intent(in) :: x
    dv03 = cos(theta) * (v1 + v2 * tau * cos(tau * (x - 0.5_dp)))
  end function dv03

  real(dp) function B02(x)
    real(dp), intent(in) :: x
    B02 = sin(phi_prof(x))
  end function B02

  real(dp) function dB02(x)
    real(dp), intent(in) :: x
    dB02 = cos(phi_prof(x)) * alpha
  end function dB02

  real(dp) function B03(x)
    real(dp), intent(in) :: x
    B03 = cos(phi_prof(x))
  end function B03

  real(dp) function dB03(x)
    real(dp), intent(in) :: x
    dB03 = -sin(phi_prof(x)) * alpha
  end function dB03

  real(dp) function p_prof(x)
    real(dp), intent(in) :: x
    p_prof = cte_p0 - (x - 0.5_dp * delta * x**2) * g
  end function p_prof

  real(dp) function phi_prof(x)
    real(dp), intent(in) :: x
    phi_prof = phi0 + alpha * (x - 0.5_dp)
  end function phi_prof

  real(dp) function v_prof(x)
    real(dp), intent(in) :: x
    v_prof = v0 + v1 * (x - 0.5_dp) + v2 * sin(tau * (x - 0.5_dp))
  end function v_prof

end submodule smod_equil_flow_driven_instabilities
