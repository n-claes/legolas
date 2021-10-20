submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_hall
  implicit none

contains

  module procedure add_natural_hall_Bterms
    use mod_global_variables, only: hall_mhd, elec_inertia
    use mod_equilibrium, only: hall_field

    real(dp)  :: eps
    real(dp)  :: rho
    real(dp)  :: eta_e

    if (.not. (hall_mhd .and. elec_inertia)) then
      return
    end if

    eps = eps_grid(grid_idx)
    rho = rho_field % rho0(grid_idx)
    eta_e = hall_field % inertiafactor(grid_idx)

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! B_H(7, 6)
    factors(1) = eta_e * k2 / rho
    positions(1, :) = [7, 6]
    ! B_H(8, 6)
    factors(2) = eta_e * k3 * eps / rho
    positions(2, :) = [8, 6]
    call subblock(quadblock, factors, positions, weight, h_cubic, h_quad)

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! B_H(7, 7)
    factors(1) = -eta_e / rho
    positions(1, :) = [7, 7]
    ! B_H(8, 8)
    factors(2) = -eta_e * eps / rho
    positions(2, :) = [8, 8]
    call subblock(quadblock, factors, positions, weight, h_cubic, dh_cubic)

  end procedure add_natural_hall_Bterms

  module procedure add_natural_hall_terms
    use mod_global_variables, only: hall_mhd, viscosity, viscosity_value
    use mod_equilibrium, only: hall_field

    real(dp)  :: eps, deps
    real(dp)  :: rho
    real(dp)  :: eta_H, mu

    if (.not. (hall_mhd .and. viscosity)) then
      return
    end if

    eps = eps_grid(grid_idx)
    deps = d_eps_grid_dr(grid_idx)
    rho = rho_field % rho0(grid_idx)
    eta_H = hall_field % hallfactor(grid_idx)
    mu = viscosity_value

    ! ==================== Quadratic * Cubic ====================
    call reset_factor_positions(new_size=1)
    ! H(6, 2)
    factors(1) = -eta_H * ic * mu * deps / (eps * rho)
    positions(1, :) = [6, 2]
    call subblock(quadblock, factors, positions, weight, h_quad, h_cubic)

    ! ==================== Quadratic * dCubic ====================
    call reset_factor_positions(new_size=1)
    ! H(6, 2)
    factors(1) = 4.0d0 * eta_H * ic * mu / (3.0d0 * rho)
    positions(1, :) = [6, 2]
    call subblock(quadblock, factors, positions, weight, h_quad, dh_cubic)

    ! ==================== Cubic * dQuadratic ====================
    call reset_factor_positions(new_size=2)
    ! H(7, 3)
    factors(1) = eta_H * ic * mu * eps / rho
    positions(1, :) = [7, 3]
    ! H(8, 4)
    factors(2) = eta_H * ic * mu / rho
    positions(2, :) = [8, 4]
    call subblock(quadblock, factors, positions, weight, h_cubic, dh_quad)

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=1)
    ! H(8, 4)
    factors(1) = -eta_H * ic * mu * deps / (eps * rho)
    positions(1, :) = [8, 4]
    call subblock(quadblock, factors, positions, weight, h_cubic, h_quad)

  end procedure add_natural_hall_terms

end submodule
