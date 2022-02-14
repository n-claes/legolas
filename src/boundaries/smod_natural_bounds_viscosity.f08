submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_viscosity
  implicit none

contains

  module procedure add_natural_viscosity_terms
    use mod_global_variables, only: viscosity, viscous_heating, viscosity_value, &
                                    gamma_1, incompressible
    use mod_equilibrium, only: v_field

    real(dp)  :: eps, deps
    real(dp)  :: mu
    real(dp)  :: dv01, dv03

    if (.not. viscosity) then
      return
    end if

    eps = eps_grid(grid_idx)
    deps = d_eps_grid_dr(grid_idx)
    mu = viscosity_value
    dv01 = v_field % d_v01_dr(grid_idx)
    dv03 = v_field % d_v03_dr(grid_idx)

    ! ==================== Cubic * Cubic ====================
    call reset_factor_positions(new_size=1)
    ! Sigma(2, 2)
    factors(1) = -ic * mu * deps / eps
    positions(1, :) = [2, 2]
    call subblock(quadblock, factors, positions, weight, h_cubic, h_cubic)

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=1)
    ! Sigma(2, 2)
    factors(1) = 4.0d0 * ic * mu / 3.0d0
    positions(1, :) = [2, 2]
    call subblock(quadblock, factors, positions, weight, h_cubic, dh_cubic)

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! Sigma(2, 3)
    factors(1) = -ic * mu * k2 / 3.0d0
    positions(1, :) = [2, 3]
    ! Sigma(2, 4)
    factors(2) = -ic * mu * k3 / 3.0d0
    positions(2, :) = [2, 4]
    call subblock(quadblock, factors, positions, weight, h_cubic, h_quad)

    ! ==================== Quadratic * dQuadratic ====================
    call reset_factor_positions(new_size=2)
    ! Sigma(3, 3)
    factors(1) = ic * mu * eps
    positions(1, :) = [3, 3]
    ! Sigma(4, 4)
    factors(2) = ic * mu
    positions(2, :) = [4, 4]
    call subblock(quadblock, factors, positions, weight, h_quad, dh_quad)

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! Sigma(4, 4)
    factors(1) = -ic * mu * deps / eps
    positions(1, :) = [4, 4]
    ! Sigma(5, 4)
    factors(2) = (0.0d0, 0.0d0)
    if (viscous_heating .and. (.not. incompressible)) then
      factors(2) = 2.0d0 * ic * gamma_1 * mu * dv03
    end if
    positions(2, :) = [5, 4]
    call subblock(quadblock, factors, positions, weight, h_quad, h_quad)

    ! ==================== Quadratic * Cubic ====================
    call reset_factor_positions(new_size=1)
    ! Sigma(5, 2)
    factors(1) = (0.0d0, 0.0d0)
    if (viscous_heating .and. (.not. incompressible)) then
      factors(1) = 2.0d0 * gamma_1 * mu * dv01
    end if
    positions(1, :) = [5, 2]
    call subblock(quadblock, factors, positions, weight, h_quad, h_cubic)

  end procedure add_natural_viscosity_terms

end submodule
