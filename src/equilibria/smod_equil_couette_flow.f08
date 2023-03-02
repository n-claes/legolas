! =============================================================================
!> This submodule defines a steady plane Couette flow in a Cartesian geometry
!! with flow and viscosity.
!! @note Default values are given by
!!
!! - <tt>k2</tt> = 0
!! - <tt>k3</tt> = 1
!! - <tt>cte_rho0</tt> = 1
!! - <tt>cte_T0</tt> = 1
!! - <tt>cte_v02</tt> = 0
!! - <tt>cte_v03</tt> = 1
!! - <tt>viscosity</tt> = True
!! - <tt>viscosity_value</tt> = 1e-3
!!
!! and can all be changed in the parfile. @endnote
submodule (mod_equilibrium) smod_equil_couette_flow
  implicit none

contains

  module procedure couette_flow_eq
    use mod_equilibrium_params, only: cte_rho0, cte_v02, cte_v03, cte_T0

    real(dp)    :: x, h
    integer     :: i

    call settings%physics%enable_flow()
    if (settings%equilibrium%use_defaults) then ! LCOV_EXCL_START
      call settings%grid%set_geometry("Cartesian")
      call settings%grid%set_grid_boundaries(0.0_dp, 1.0_dp)
      cte_v02  = 0.0d0
      cte_v03  = 1.0d0
      cte_T0   = 1.0d0
      cte_rho0 = 1.0d0

      k2 = 0.0d0
      k3 = 1.0d0

      call settings%physics%enable_viscosity(viscosity_value=0.001_dp)
    end if ! LCOV_EXCL_STOP

    call initialise_grid(settings)
    h = settings%grid%get_grid_end() - settings%grid%get_grid_start()

    do i = 1, settings%grid%get_gauss_gridpts()
      x = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0
      v_field % v02(i)    = cte_v02 * x / h
      v_field % v03(i)    = cte_v03 * x / h
      T_field % T0(i)     = cte_T0

      v_field % d_v02_dr(i)    = cte_v02 / h
      v_field % d_v03_dr(i)    = cte_v03 / h
    end do
  end procedure couette_flow_eq

end submodule smod_equil_couette_flow
