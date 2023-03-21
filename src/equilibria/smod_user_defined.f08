!> Submodule for user-defined equilibria.
!! Look at the examples in the equilibria subdirectory or consult the website for
!! more information.
submodule (mod_equilibrium) smod_user_defined
  implicit none

contains

  ! LCOV_EXCL_START <exclude this file from code coverage>
  module procedure user_defined_eq
    real(dp)    :: x
    integer     :: i

    call settings%grid%set_geometry("Cartesian")
    call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
    call initialise_grid(settings)

    k2 = 0.0_dp
    k3 = 1.0_dp

    ! additional physics
    call settings%physics%enable_flow()

    ! set up the grid
    do i = 1, settings%grid%get_gauss_gridpts()
      x = grid_gauss(i)

      ! Note: values that are not set/referenced are automatically set to zero.

      rho_field%rho0(i) = rho0(x)
      rho_field%d_rho0_dr(i) = drho0(x)
      T_field%T0(i) = T0()
      v_field%v02(i) = v02(x)
      v_field%d_v02_dr(i) = dv02(x)
      v_field%dd_v02_dr(i) = ddv02()
      B_field%B03(i) = B03(x)
      B_field%B0(i) = sqrt(B_field%B02(i)**2 + B_field%B03(i)**2)
    end do

    call background%set_density_funcs(rho0_func=rho0, drho0_func=drho0)
    call background%set_velocity_2_funcs(v02_func=v02, dv02_func=dv02, ddv02_func=ddv02)
    call background%set_temperature_funcs(T0_func=T0)
    call background%set_magnetic_3_funcs(B03_func=B03)
  end procedure user_defined_eq
  ! LCOV_EXCL_STOP

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

end submodule smod_user_defined
