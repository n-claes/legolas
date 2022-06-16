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

  module subroutine couette_flow_eq()
    use mod_equilibrium_params, only: cte_rho0, cte_v02, cte_v03, cte_T0
    use mod_global_variables, only: viscosity_value

    real(dp)    :: x, h
    integer     :: i

    call allow_geometry_override( &
      default_geometry='Cartesian', default_x_start=0.0d0, default_x_end=1.0d0 &
    )
    call initialise_grid()

    flow = .true.
    if (use_defaults) then ! LCOV_EXCL_START
      cte_v02  = 0.0d0
      cte_v03  = 1.0d0
      cte_T0   = 1.0d0
      cte_rho0 = 1.0d0

      k2 = 0.0d0
      k3 = 1.0d0

      viscosity = .true.
      viscosity_value = 1.0d-3
    end if ! LCOV_EXCL_STOP

    h = x_end - x_start

    do i = 1, gauss_gridpts
      x = grid_gauss(i)

      rho_field % rho0(i) = cte_rho0
      v_field % v02(i)    = cte_v02 * x / h
      v_field % v03(i)    = cte_v03 * x / h
      T_field % T0(i)     = cte_T0

      v_field % d_v02_dr(i)    = cte_v02 / h
      v_field % d_v03_dr(i)    = cte_v03 / h
    end do
  end subroutine couette_flow_eq

end submodule smod_equil_couette_flow
