submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_flow
  implicit none

contains

  module procedure add_natural_flow_terms
    use mod_global_variables, only: flow
    use mod_equilibrium, only: v_field

    real(dp)  :: rho
    real(dp)  :: v01

    if (.not. flow) then
      return
    end if

    rho = rho_field % rho0(grid_idx)
    v01 = v_field % v01(grid_idx)

    ! ==================== Cubic * Cubic ====================
    call reset_factor_positions(new_size=1)
    ! Phi(2, 2)
    factors(1) = -ic * rho * v01
    positions(1, :) = [2, 2]
    call subblock(quadblock, factors, positions, weight, h_cubic, h_cubic)

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! Phi(4, 4)
    factors(1) = -ic * rho * v01
    positions(1, :) = [4, 4]
    ! Phi(5, 5)
    factors(2) = -ic * rho * v01
    positions(2, :) = [5, 5]
    call subblock(quadblock, factors, positions, weight, h_quad, h_quad)
  end procedure add_natural_flow_terms

end submodule smod_natural_bounds_flow
