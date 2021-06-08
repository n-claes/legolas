submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_regular
  implicit none

contains

  module procedure add_natural_regular_terms
    use mod_matrix_shortcuts, only: get_G_operator

    real(dp)  :: eps
    real(dp)  :: rho, T0
    real(dp)  :: B01, B02, B03
    real(dp)  :: Gop_min

    eps = eps_grid(grid_idx)
    rho = rho_field % rho0(grid_idx)
    T0 = T_field % T0(grid_idx)
    B01 = B_field % B01
    B02 = B_field % B02(grid_idx)
    B03 = B_field % B03(grid_idx)
    Gop_min = get_G_operator(grid_idx, which="minus")

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=3)
    ! A(2, 1)
    factors(1) = T0
    positions(1, :) = [2, 1]
    ! A(2, 5)
    factors(2) = rho
    positions(2, :) = [2, 5]
    ! A(2, 6)
    factors(3) = eps * Gop_min
    positions(3, :) = [2, 6]
    call subblock(quadblock, factors, positions, weight, h_cubic, h_quad)

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! A(2, 7)
    factors(1) = B03
    positions(1, :) = [2, 7]
    ! A(2, 8)
    factors(2) = -eps * B02
    positions(2, :) = [2, 8]
    call subblock(quadblock, factors, positions, weight, h_cubic, dh_cubic)

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! A(3, 6)
    factors(1) = ic * eps * k3 * B01
    positions(1, :) = [3, 6]
    ! A(4, 6)
    factors(2) = -ic * k2 * B01
    positions(2, :) = [4, 6]
    call subblock(quadblock, factors, positions, weight, h_quad, h_quad)

    ! ==================== Quadratic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! A(3, 8)
    factors(1) = -ic * eps * B01
    positions(1, :) = [3, 8]
    ! A(4, 7)
    factors(2) = ic * B01
    positions(2, :) = [4, 7]
    call subblock(quadblock, factors, positions, weight, h_quad, dh_cubic)
  end procedure add_natural_regular_terms

end submodule smod_natural_bounds_regular
