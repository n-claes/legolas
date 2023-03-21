! =============================================================================
!> This submodule defines a steady Taylor-Couette flow in a cylindrical geometry
!! where a fluid is confined between two (rotating) coaxial cylinders
!! (without a magnetic field).
!!
!! This equilibrium is taken from
!! _Gebhardt, Thomas and Grossman, Siegfried.
!! "The Taylor-Couette eigenvalue problem with independently rotating cylinders.",
!! Z. Phys. B 90, 475--490 (1993)_.
!!
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 1
!! - <tt>cte_rho0</tt> = 1 : density (constant)
!! - <tt>alpha</tt> = 1 : rotational speed of the inner cylinder
!! - <tt>beta</tt> = 2 : rotational speed of the outer cylinder

!! and can all be changed in the parfile. @endnote
! SUBMODULE: smod_equil_taylor_couette
submodule (mod_equilibrium) smod_equil_taylor_couette
  use mod_equilibrium_params, only: cte_rho0, alpha, beta
implicit none

  real(dp) :: h, Rrat, A, B, Tstart
  real(dp) :: x_start, x_end

contains

  module procedure taylor_couette_eq

    real(dp) :: r, Ta, grid_middle
    real(dp) :: viscosity_value
    integer :: i, gauss_gridpts

    call settings%physics%enable_flow()
    settings%grid%coaxial = .true.

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("cylindrical")
      call settings%grid%set_grid_boundaries(1.0_dp, 2.0_dp)
      cte_rho0 = 1.0_dp
      alpha = 1.0_dp
      beta = 2.0_dp
      k2 = 0.0_dp
      k3 = 1.0_dp

      call settings%physics%enable_viscosity(viscosity_value=0.001_dp)
    end if ! LCOV_EXCL_STOP

    call initialise_grid(settings)
    x_start = settings%grid%get_grid_start()
    x_end = settings%grid%get_grid_end()
    gauss_gridpts = settings%grid%get_gauss_gridpts()

    viscosity_value = settings%physics%viscosity%get_viscosity_value()
    h = x_end - x_start
    Rrat = x_start / x_end
    A = (alpha * Rrat**2 - beta) / (Rrat**2 - 1.0_dp)
    B = x_start**2 * (alpha - beta) / (1.0_dp - Rrat**2)
    Tstart = 0.5_dp * ( &
      (A * x_start)**2 + 4.0_dp * A * B * log(x_start) - (B / x_start)**2 &
    )

    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      rho_field%rho0(i) = rho0()
      v_field % v02(i) = v02(r)
      v_field % d_v02_dr(i) = dv02(r)
      v_field % dd_v02_dr(i) = ddv02(r)
      T_field % T0(i) = T0(r)
      T_field % d_T0_dr(i) = dT0(r)
    end do

    grid_middle = grid_gauss(int(gauss_gridpts/2))
    Ta = ( &
      cte_rho0 * v02(grid_middle) * h / viscosity_value &
    )**2 * 2.0_dp * h / (x_start + x_end)
    call logger%info('Taylor number is ' // str(int(Ta)))
  end procedure taylor_couette_eq


  real(dp) function rho0()
    rho0 = cte_rho0
  end function rho0

  real(dp) function T0(r)
    real(dp), intent(in) :: r
    if (Tstart > 0.0_dp) then
      T0 = 0.5_dp * ((A * r)**2 + 4.0_dp * A * B * log(r) - (B / r)**2)
    else
      T0 = 2.0_dp * abs(Tstart) + 0.5_dp * ( &
        (A * r)**2 + 4.0_dp * A * B * log(r) - (B / r)**2 &
      )
    end if
  end function T0

  real(dp) function dT0(r)
    real(dp), intent(in) :: r
    dT0 = v02(r)**2 / r
  end function dT0

  real(dp) function v02(r)
    real(dp), intent(in) :: r
    v02 = A * r + B / r
  end function v02

  real(dp) function dv02(r)
    real(dp), intent(in) :: r
    dv02 = A - B / r**2
  end function dv02

  real(dp) function ddv02(r)
    real(dp), intent(in) :: r
    ddv02 = 2.0_dp * B / r**3
  end function ddv02

end submodule smod_equil_taylor_couette
