!> Submodule for user-defined equilibria.
!! Look at the examples in the equilibria subdirectory or consult the website for
!! more information.
submodule (mod_equilibrium) smod_user_defined
  implicit none

contains

  ! LCOV_EXCL_START <exclude this file from code coverage>
  module procedure user_defined_eq
    call settings%grid%set_geometry("Cartesian")
    call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)

    k2 = 0.0_dp
    k3 = 1.0_dp

    ! additional physics
    call settings%physics%enable_flow()
    call settings%physics%enable_gravity()

    ! Note: functions that are not set are automatically set to zero.
    call background%set_density_funcs(rho0_func=rho0, drho0_func=drho0)
    call background%set_velocity_2_funcs(v02_func=v02, dv02_func=dv02, ddv02_func=ddv02)
    call background%set_temperature_funcs(T0_func=T0)
    call background%set_magnetic_3_funcs(B03_func=B03)

    call physics%set_gravity_funcs(g0_func=g0)
  end procedure user_defined_eq

  real(dp) function rho0(x)
    real(dp), intent(in) :: x
    rho0 = x**2
  end function rho0

  real(dp) function drho0(x)
    real(dp), intent(in) :: x
    drho0 = 2.0_dp * x
  end function drho0

  real(dp) function v02(x)
    real(dp), intent(in) :: x
    v02 = 2.0_dp * x**2
  end function v02

  real(dp) function dv02(x)
    real(dp), intent(in) :: x
    dv02 = 4.0_dp * x
  end function dv02

  real(dp) function ddv02()
    ddv02 = 4.0_dp
  end function ddv02

  real(dp) function T0()
    T0 = 1.0_dp
  end function T0

  real(dp) function B03(x)
    real(dp), intent(in) :: x
    B03 = x
  end function B03

  real(dp) function g0(x, settings, background)
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    g0 = 1.0_dp / x**2
  end function g0

  ! LCOV_EXCL_STOP
end submodule smod_user_defined
