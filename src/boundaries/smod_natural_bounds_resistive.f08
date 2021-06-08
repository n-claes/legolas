submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_resistive
  implicit none

contains

  module procedure add_natural_resistive_terms
    use mod_global_variables, only: gamma_1, resistivity
    use mod_equilibrium, only: eta_field

    real(dp)  :: eps, deps
    real(dp)  :: eta
    real(dp)  :: B02, dB02, drB02
    real(dp)  :: B03, dB03

    if (.not. resistivity) then
      return
    end if

    eps = eps_grid(grid_idx)
    deps = d_eps_grid_dr(grid_idx)
    eta = eta_field % eta(grid_idx)
    B02 = B_field % B02(grid_idx)
    dB02 = B_field % d_B02_dr(grid_idx)
    B03 = B_field % B03(grid_idx)
    dB03 = B_field % d_B03_dr(grid_idx)

    drB02 = deps * B02 + eps * dB02

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=1)
    ! R(5, 6)
    factors(1) = 2.0d0 * ic * gamma_1 * eta * (k3 * drB02 - k2 * dB03)
    positions(1, :) = [5, 6]
    call subblock(quadblock, factors, positions, weight, h_quad, h_quad)

    ! ==================== Quadratic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! R(5, 7)
    factors(1) = 2.0d0 * ic * gamma_1 * eta * dB03
    positions(1, :) = [5, 7]
    ! R(5, 8)
    factors(2) = -2.0d0 * ic * gamma_1 * eta * drB02
    positions(2, :) = [5, 8]
    call subblock(quadblock, factors, positions, weight, h_quad, dh_cubic)

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! R(7, 6)
    factors(1) = -ic * eta * k2
    positions(1, :) = [7, 6]
    ! R(8, 6)
    factors(2) = -ic * eta * eps * k3
    positions(2, :) = [8, 6]
    call subblock(quadblock, factors, positions, weight, h_cubic, h_quad)

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! R(7, 7)
    factors(1) = ic * eta
    positions(1, :) = [7, 7]
    ! R(8, 8)
    factors(2) = ic * eta * eps
    positions(2, :) = [8, 8]
    call subblock(quadblock, factors, positions, weight, h_cubic, dh_cubic)
  end procedure add_natural_resistive_terms

end submodule smod_natural_bounds_resistive
