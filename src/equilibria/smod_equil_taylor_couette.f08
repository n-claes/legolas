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
  implicit none

contains

  module procedure taylor_couette_eq
    use mod_equilibrium_params, only: cte_rho0, alpha, beta

    real(dp) :: r, h, Rrat, A, B, Ta, Tstart
    real(dp) :: x_start, x_end
    real(dp) :: viscosity_value
    integer :: i, gauss_gridpts

    call settings%physics%enable_flow()
    settings%grid%coaxial = .true.

    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("cylindrical")
      call settings%grid%set_grid_boundaries(0.0_dp, 2.0_dp)
      cte_rho0 = 1.0d0
      alpha = 1.0d0
      beta = 2.0d0

      k2 = 0.0d0
      k3 = 1.0d0

      call settings%physics%enable_viscosity(viscosity_value=0.001_dp)
    end if ! LCOV_EXCL_STOP
    x_start = settings%grid%get_grid_start()
    x_end = settings%grid%get_grid_end()
    gauss_gridpts = settings%grid%get_gauss_gridpts()
    call initialise_grid(settings)

    viscosity_value = settings%physics%viscosity%get_viscosity_value()

    rho_field % rho0 = cte_rho0

    h = x_end - x_start
    Rrat = x_start / x_end
    A = (alpha * Rrat**2 - beta) / (Rrat**2 - 1.0d0)
    B = x_start**2 * (alpha - beta) / (1.0d0 - Rrat**2)

    Tstart = 0.5d0 * ( &
      (A * x_start)**2 + 4.0d0 * A * B * log(x_start) - (B / x_start)**2 &
    )
    do i = 1, gauss_gridpts
      r = grid_gauss(i)

      v_field % v02(i) = A * r + B / r
      v_field % d_v02_dr(i) = A - B / r**2
      v_field % dd_v02_dr(i) = 2.0d0 * B / r**3

      if (Tstart > 0) then
        T_field % T0(i) = 0.5d0 * ((A * r)**2 + 4.0d0 * A * B * log(r) - (B / r)**2)
      else
        T_field % T0(i) = 2.0d0 * abs(Tstart) + 0.5d0 * ( &
          (A * r)**2 + 4.0d0 * A * B * log(r) - (B / r)**2 &
        )
      end if

      T_field % d_T0_dr(i) = (v_field % v02(i))**2 / r
    end do

    Ta = (cte_rho0 * (v_field % v02(int(gauss_gridpts/2))) * h / viscosity_value)**2 &
          * 2.0d0 * h / (x_start + x_end)
    call logger%info('Taylor number is ' // str(int(Ta)))

  end procedure taylor_couette_eq

end submodule smod_equil_taylor_couette
