submodule (mod_matrix_manager) smod_regular_matrix
  implicit none

contains

  module procedure add_bmatrix_terms
    real(dp)  :: rho, eps

    rho = rho_field % rho0(gauss_idx)
    eps = eps_grid(gauss_idx)

    ! Quadratic * Quadratic
    call reset_factor_positions(new_size=5)
    ! B(1, 1)
    factors(1) = 1.0d0
    positions(1, :) = [1, 1]
    ! B(3, 3)
    factors(2) = eps * rho
    positions(2, :) = [3, 3]
    ! B(4, 4)
    factors(3) = rho
    positions(3, :) = [4, 4]
    ! B(5, 5)
    factors(4) = rho
    positions(4, :) = [5, 5]
    ! B(6, 6)
    factors(5) = eps
    positions(5, :) = [6, 6]
    call subblock( &
      quadblock, factors, positions, current_weight, h_quad, h_quad, settings%dims &
    )

    ! Cubic * Cubic
    call reset_factor_positions(new_size=3)
    ! B(2, 2)
    factors(1) = rho
    positions(1, :) = [2, 2]
    ! B(7, 7)
    factors(2) = 1.0d0
    positions(2, :) = [7, 7]
    ! B(8, 8)
    factors(3) = eps
    positions(3, :) = [8, 8]
    call subblock( &
      quadblock, factors, positions, current_weight, h_cubic, h_cubic, settings%dims &
    )
  end procedure add_bmatrix_terms


  module procedure add_regular_matrix_terms
    use mod_equilibrium, only: grav_field

    real(dp)  :: eps, deps
    real(dp)  :: rho, drho
    real(dp)  :: T0, dT0
    real(dp)  :: B01, B02, dB02, drB02, B03, db03
    real(dp)  :: Fop_plus, Gop_plus, Gop_min, WVop
    real(dp) :: gamma_1

    gamma_1 = settings%physics%get_gamma_1()

    ! grid variables
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! density variables
    rho = rho_field % rho0(gauss_idx)
    drho = rho_field % d_rho0_dr(gauss_idx)
    ! temperature variables
    T0 = T_field % T0(gauss_idx)
    dT0 = T_field % d_T0_dr(gauss_idx)
    ! magnetic field variables
    B01 = B_field % B01
    B02 = B_field % B02(gauss_idx)
    dB02 = B_field % d_B02_dr(gauss_idx)
    drB02 = deps * B02 + eps * dB02
    B03 = B_field % B03(gauss_idx)
    dB03 = B_field % d_B03_dr(gauss_idx)
    ! operators
    Fop_plus = get_F_operator(gauss_idx, which="plus")
    Gop_plus = get_G_operator(gauss_idx, which="plus")
    Gop_min = get_G_operator(gauss_idx, which="minus")
    WVop = get_wv_operator(gauss_idx)

    ! ==================== Quadratic * Cubic ====================
    call reset_factor_positions(new_size=6)
    ! A(1, 2)
    factors(1) = -drho
    positions(1, :) = [1, 2]
    ! A(3, 7)
    factors(2) = k3 * (drB02 - ic * k2 * B01) / eps
    positions(2, :) = [3, 7]
    ! A(3, 8)
    factors(3) = k2 * (ic * k2 * B01 - drB02) / eps
    positions(3, :) = [3, 8]
    ! A(4, 7)
    factors(4) = k3 * (dB03 - ic * k3 * B01)
    positions(4, :) = [4, 7]
    ! A(4, 8)
    factors(5) = k2 * (ic * B01 * k3 - dB03)
    positions(5, :) = [4, 8]
    ! A(5, 2)
    factors(6) = 0.0d0
    if (.not. settings%physics%is_incompressible) then
      factors(6) = -dT0 * rho
    end if
    positions(6, :) = [5, 2]
    call subblock( &
      quadblock, factors, positions, current_weight, h_quad, h_cubic, settings%dims &
    )

    ! ==================== Quadratic * dCubic ====================
    call reset_factor_positions(new_size=6)
    ! A(1, 2)
    factors(1) = -rho
    positions(1, :) = [1, 2]
    ! A(3, 7)
    factors(2) = k2 * B03 / eps
    positions(2, :) = [3, 7]
    ! A(3, 8)
    factors(3) = eps * k3 * B03
    positions(3, :) = [3, 8]
    ! A(4, 7)
    factors(4) = -(k2 * B02 + ic * deps * B01) / eps
    positions(4, :) = [4, 7]
    ! A(4, 8)
    factors(5) = -eps * k3 * B02
    positions(5, :) = [4, 8]
    ! A(5, 2)
    factors(6) = -gamma_1 * T0 * rho
    positions(6, :) = [5, 2]
    call subblock( &
      quadblock, factors, positions, current_weight, h_quad, dh_cubic, settings%dims &
    )

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=12)
    ! A(1, 3)
    factors(1) = rho * k2
    positions(1, :) = [1, 3]
    ! A(1, 4)
    factors(2) = rho * k3
    positions(2, :) = [1, 4]
    ! A(3, 1)
    factors(3) = k2 * T0 / eps
    positions(3, :) = [3, 1]
    ! A(3, 5)
    factors(4) = k2 * rho / eps
    positions(4, :) = [3, 5]
    ! A(3, 6)
    factors(5) = -WVop * B03
    positions(5, :) = [3, 6]
    ! A(4, 1)
    factors(6) = k3 * T0
    positions(6, :) = [4, 1]
    ! A(4, 5)
    factors(7) = k3 * rho
    positions(7, :) = [4, 5]
    ! A(4, 6)
    factors(8) = ic * deps * k2 * B01 / eps + B02 * WVop
    positions(8, :) = [4, 6]
    ! A(5, 3)
    factors(9) = gamma_1 * k2 * rho * T0
    positions(9, :) = [5, 3]
    ! A(5, 4)
    factors(10) = gamma_1 * k3 * rho * T0
    positions(10, :) = [5, 4]
    ! A(6, 3)
    factors(11) = -eps * B03
    positions(11, :) = [6, 3]
    ! A(6, 4)
    factors(12) = B02
    positions(12, :) = [6, 4]
    call subblock( &
      quadblock, factors, positions, current_weight, h_quad, h_quad, settings%dims &
    )

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=5)
    ! A(2, 1)
    factors(1) = -deps * T0 / eps
    if (settings%physics%gravity%is_enabled()) then
      ! adds gravity term to A(2, 1) matrix element
      factors(1) = factors(1) + grav_field % grav(gauss_idx)
    end if
    positions(1, :) = [2, 1]
    ! A(2, 5)
    factors(2) = -deps * rho / eps
    positions(2, :) = [2, 5]
    ! A(2, 6)
    factors(3) = deps * Gop_plus
    positions(3, :) = [2, 6]
    ! A(7, 4)
    factors(4) = ic * B01
    positions(4, :) = [7, 4]
    ! A(8, 3)
    factors(5) = -ic * eps * B01
    positions(5, :) = [8, 3]
    call subblock( &
      quadblock, factors, positions, current_weight, h_cubic, h_quad, settings%dims &
    )

    ! ==================== dCubic * Quadratic ====================
    call reset_factor_positions(new_size=3)
    ! A(2, 1)
    factors(1) = -T0
    positions(1, :) = [2, 1]
    ! A(2, 5)
    factors(2) = -rho
    positions(2, :) = [2, 5]
    ! A(2, 6)
    factors(3) = -eps * Gop_min
    positions(3, :) = [2, 6]
    call subblock( &
      quadblock, factors, positions, current_weight, dh_cubic, h_quad, settings%dims &
    )

    ! ==================== Cubic * Cubic ====================
    call reset_factor_positions(new_size=4)
    ! A(2, 7)
    factors(1) = -k3 * Fop_plus
    positions(1, :) = [2, 7]
    ! A(2, 8)
    factors(2) = k2 * Fop_plus
    positions(2, :) = [2, 8]
    ! A(7, 2)
    factors(3) = -B03
    positions(3, :) = [7, 2]
    ! A(8, 2)
    factors(4) = B02
    positions(4, :) = [8, 2]
    call subblock( &
      quadblock, factors, positions, current_weight, h_cubic, h_cubic, settings%dims &
    )

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! A(2, 7)
    factors(1) = -deps * B03 / eps
    positions(1, :) = [2, 7]
    ! A(2, 8)
    factors(2) = -deps * B02
    positions(2, :) = [2, 8]
    call subblock( &
      quadblock, factors, positions, current_weight, h_cubic, dh_cubic, settings%dims &
    )

    ! ==================== dCubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! A(2, 7)
    factors(1) = -B03
    positions(1, :) = [2, 7]
    ! A(2, 8)
    factors(2) = eps * B02
    positions(2, :) = [2, 8]
    call subblock( &
      quadblock, factors, positions, current_weight, dh_cubic, dh_cubic, settings%dims &
    )

    ! ==================== dQuadratic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! A(3, 6)
    factors(1) = -ic * eps * k3 * B01
    positions(1, :) = [3, 6]
    ! A(4, 6)
    factors(2) = ic * k2 * B01
    positions(2, :) = [4, 6]
    call subblock( &
      quadblock, factors, positions, current_weight, dh_quad, h_quad, settings%dims &
    )

    ! ==================== dQuadratic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! A(3, 8)
    factors(1) = ic * eps * B01
    positions(1, :) = [3, 8]
    ! A(4, 7)
    factors(2) = -ic * B01
    positions(2, :) = [4, 7]
    call subblock( &
      quadblock, factors, positions, current_weight, dh_quad, dh_cubic, settings%dims &
    )

  end procedure add_regular_matrix_terms

end submodule smod_regular_matrix
