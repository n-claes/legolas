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
    call settings%grid%set_geometry("Cartesian")
    call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
    call initialise_grid(settings)
    call settings%physics%enable_flow()
    v0 = p1
    v1 = p2
    v2 = p3
    phi0 = p4

    grav_field%grav = g

    call background%set_density_funcs(rho0_func=rho0, drho0_func=drho0)
    call background%set_velocity_2_funcs(v02_func=v02, dv02_func=dv02)
    call background%set_velocity_3_funcs(v03_func=v03, dv03_func=dv03)
    call background%set_temperature_funcs(T0_func=T0, dT0_func=dT0)
    call background%set_magnetic_2_funcs(B02_func=B02, dB02_func=dB02)
    call background%set_magnetic_3_funcs(B03_func=B03, dB03_func=dB03)
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
