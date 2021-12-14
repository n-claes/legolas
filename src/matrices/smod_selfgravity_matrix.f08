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
    use mod_selfgravity, only: get_gravity_prefactor

    real(dp)  :: rho, drho, eps
    real(dp)  :: v01, v02, v03, dv01
    real(dp)  :: gravity_prefactor

    rho = rho_field % rho0(gauss_idx)
    drho = rho_field % d_rho0_dr(gauss_idx)
    eps = eps_grid(gauss_idx)
    v01 = v_field % v01(gauss_idx)
    dv01 = v_field % d_v01_dr(gauss_idx)
    v02 = v_field % v03(gauss_idx)
    v03 = v_field % v03(gauss_idx)

    gravity_prefactor = get_gravity_prefactor()

    ! Cubic * Quadratic
    call reset_factor_positions(new_size=3)
    ! G(9, 1)
    factors(1) = k2 * v02 / eps + k3 * v03 - ic * dv01
    positions(1, :) = [9, 1]
    ! G(9, 3)
    factors(2) = rho * k2
    positions(2, :) = [9, 3]
    ! G(9, 4)
    factors(3) = rho * k3
    positions(3, :) = [9, 4]
    ! multiply with 4 pi G
    factors = factors * gravity_prefactor
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_quad)

    ! Cubic * dCubic
    call reset_factor_positions(new_size=2)
    ! G(2, 9)
    factors(1) = eps * rho
    positions(1, :) = [2, 9]
    ! G(9, 2)
    factors(2) = -gravity_prefactor * rho
    positions(2, :) = [9, 2]
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

    ! Cubic * Cubic
    call reset_factor_positions(new_size=1)
    ! G(9, 2)
    factors(1) = -gravity_prefactor * drho
    positions(1, :) = [9, 2]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_cubic)

    ! Cubic * dQuadratic
    call reset_factor_positions(new_size=1)
    ! G(9, 1)
    factors(1) = -gravity_prefactor * ic * v01
    positions(1, :) = [9, 1]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, dh_quad)
  end procedure add_selfgravity_terms

end submodule smod_selfgravity_matrix
