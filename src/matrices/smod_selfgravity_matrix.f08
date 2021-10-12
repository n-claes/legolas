submodule (mod_matrix_manager) smod_selfgravity_matrix
  implicit none

contains

  module procedure add_selfgravity_bmatrix_terms
    real(dp)  :: eps

    eps = eps_grid(gauss_idx)

    ! Cubic * Cubic
    call reset_factor_positions(new_size=1)
    ! B(9, 9)
    factors(1) = -(k2**2 / eps + eps * k3**2)
    positions(1, :) = [9, 9]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_cubic)

    ! dCubic * dCubic
    call reset_factor_positions(new_size=1)
    ! B(9, 9)
    factors(1) = -eps
    positions(1, :) = [9, 9]
    call subblock(quadblock, factors, positions, current_weight, dh_cubic, dh_cubic)
  end procedure add_selfgravity_bmatrix_terms


  module procedure add_selfgravity_terms
    use mod_equilibrium, only: v_field

    real(dp)  :: rho, eps
    real(dp)  :: v01, v02, v03

    rho = rho_field % rho0(gauss_idx)
    eps = eps_grid(gauss_idx)

    v01 = v_field % v01(gauss_idx)
    v02 = v_field % v03(gauss_idx)
    v03 = v_field % v03(gauss_idx)

    ! Cubic * Quadratic
    call reset_factor_positions(new_size=3)
    ! G(9, 1)
    factors(1) = k2 * v02 / eps + k3 * v03
    positions(1, :) = [9, 1]
    ! G(9, 3)
    factors(2) = rho * k2
    positions(2, :) = [9, 3]
    ! G(9, 4)
    factors(3) = rho * k3
    positions(3, :) = [9, 4]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_quad)

    ! dCubic * Quadratic
    call reset_factor_positions(new_size=1)
    ! G(9, 1)
    factors(1) = ic * v01
    positions(1, :) = [9, 1]
    call subblock(quadblock, factors, positions, current_weight, dh_cubic, h_quad)

    ! dCubic * Cubic
    call reset_factor_positions(new_size=1)
    ! G(9, 2)
    factors(1) = rho
    positions(1, :) = [9, 2]
    call subblock(quadblock, factors, positions, current_weight, dh_cubic, h_cubic)

    ! Cubic * dCubic
    call reset_factor_positions(new_size=1)
    ! G(2, 9)
    factors(1) = eps * rho
    positions(1, :) = [9, 2]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, dh_cubic)

    ! Quadratic * Cubic
    call reset_factor_positions(new_size=2)
    ! G(3, 9)
    factors(1) = rho * k2
    positions(1, :) = [3, 9]
    ! G(4, 9)
    factors(2) = eps * rho * k3
    positions(2, :) = [4, 9]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_cubic)
  end procedure add_selfgravity_terms

end submodule smod_selfgravity_matrix
