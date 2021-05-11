submodule (mod_matrix_manager) smod_hall_matrix
  use mod_equilibrium, only: hall_field
  implicit none

contains

  module procedure add_hall_bmatrix_terms
    use mod_global_variables, only: elec_inertia

    real(dp)  :: eps, deps
    real(dp)  :: rho, drho
    real(dp)  :: eta_e
    real(dp)  :: WVop

    if (elec_inertia) then
      eps = eps_grid(gauss_idx)
      deps = d_eps_grid_dr(gauss_idx)
      rho = rho_field % rho0(gauss_idx)
      drho = rho_field % d_rho0_dr(gauss_idx)
      eta_e = hall_field % inertiafactor(gauss_idx)
      WVop = get_wv_operator(gauss_idx)

      ! ==================== Quadratic * Quadratic ====================
      call reset_factor_positions(new_size=1)
      ! B_H(6, 6)
      factors(1) = eta_e * WVop / rho
      positions(1, :) = [6, 6]
      call subblock(quadblock, factors, positions, current_weight, h_quad, h_quad)

      ! ==================== Quadratic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(6, 7)
      factors(1) = -eta_e * k2 / (eps * rho)
      positions(1, :) = [6, 7]
      ! B_H(6, 8)
      factors(2) = -eta_e * eps * k3 / rho
      positions(2, :) = [6, 8]
      call subblock(quadblock, factors, positions, current_weight, h_quad, dh_cubic)

      ! ==================== Cubic * Quadratic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(7, 6)
      factors(1) = -eta_e * k2 * (deps / (eps * rho) - drho / rho**2)
      positions(1, :) = [7, 6]
      ! B_H(8, 6)
      factors(2) = eta_e * drho * eps * k3 / rho**2
      positions(2, :) = [8, 6]
      call subblock(quadblock, factors, positions, current_weight, h_cubic, h_quad)

      ! ==================== dCubic * Quadratic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(7, 6)
      factors(1) = -eta_e * k2 / rho
      positions(1, :) = [7, 6]
      ! B_H(8, 6)
      factors(2) = -eta_e * eps * k3 / rho
      positions(2, :) = [8, 6]
      call subblock(quadblock, factors, positions, current_weight, dh_cubic, h_quad)

      ! ==================== Cubic * Cubic ====================
      call reset_factor_positions(new_size=4)
      ! B_H(7, 7)
      factors(1) = eta_e * k3**2 / rho
      positions(1, :) = [7, 7]
      ! B_H(7, 8)
      factors(2) = -eta_e * k2 * k3 / rho
      positions(2, :) = [7, 8]
      ! B_H(8, 7)
      factors(3) = -eta_e * k2 * k3 / (eps * rho)
      positions(3, :) = [8, 7]
      ! B_H(8, 8)
      factors(4) = eta_e * k2**2 / (eps * rho)
      positions(4, :) = [8, 8]
      call subblock(quadblock, factors, positions, current_weight, h_cubic, h_cubic)

      ! ==================== Cubic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(7, 7)
      factors(1) = eta_e * (deps / (eps * rho) - drho / rho**2)
      positions(1, :) = [7, 7]
      ! B_H(8, 8)
      factors(2) = -eta_e * eps * drho / rho**2
      positions(2, :) = [8, 8]
      call subblock(quadblock, factors, positions, current_weight, h_cubic, dh_cubic)

      ! ==================== dCubic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(7, 7)
      factors(1) = eta_e / rho
      positions(1, :) = [7, 7]
      ! B_H(8, 8)
      factors(2) = eta_e * eps / rho
      positions(2, :) = [8, 8]
      call subblock(quadblock, factors, positions, current_weight, dh_cubic, dh_cubic)
    else
      return
    end if

  end procedure add_hall_bmatrix_terms


  module procedure add_hall_matrix_terms
    use mod_global_variables, only: elec_pressure

    real(dp)  :: eps, deps
    real(dp)  :: rho, drho
    real(dp)  :: T0, dT0
    real(dp)  :: B01, B02, B03, dB02, dB03
    real(dp)  :: drB02, dB02_r
    real(dp)  :: eta_H
    real(dp)  :: Fop_plus, Gop_plus, Gop_min, WVop

    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)

    rho = rho_field % rho0(gauss_idx)
    drho = rho_field % d_rho0_dr(gauss_idx)

    T0 = T_field % T0(gauss_idx)
    dT0 = T_field % d_T0_dr(gauss_idx)

    B01 = B_field % B01
    B02 = B_field % B02(gauss_idx)
    dB02 = B_field % d_B02_dr(gauss_idx)
    B03 = B_field % B03(gauss_idx)
    dB03 = B_field % d_B03_dr(gauss_idx)

    ! Calculate derivatives eps*B02, B02/eps
    drB02 = deps * B02 + eps * dB02
    dB02_r = dB02 / eps - B02 * deps / eps**2

    eta_H = hall_field % hallfactor(gauss_idx)

    Fop_plus = get_F_operator(gauss_idx, which="plus")
    Gop_plus = get_G_operator(gauss_idx, which="plus")
    Gop_min = get_G_operator(gauss_idx, which="minus")
    WVop = get_wv_operator(gauss_idx)

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=1)
    ! H(6, 6)
    factors(1) = -eta_H * (eps * drho * Gop_min / rho + deps * Gop_plus) / rho
    positions(1, :) = [6, 6]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_quad)
    if (elec_pressure) then
      call reset_factor_positions(new_size=2)
      ! H(6, 1)
      factors(1) = eta_H * T0 * (deps / eps - drho / rho) / rho
      positions(1, :) = [6, 1]
      ! H(6, 5)
      factors(2) = eta_H * (deps / eps - drho / rho)
      positions(2, :) = [6, 5]
      call subblock(quadblock, factors, positions, current_weight, h_quad, h_quad)
    end if

    ! ==================== dQuadratic * Quadratic ====================
    call reset_factor_positions(new_size=1)
    ! H(6, 6)
    factors(1) = eta_H * eps * Gop_min / rho
    positions(1, :) = [6, 6]
    call subblock(quadblock, factors, positions, current_weight, dh_quad, h_quad)
    if (elec_pressure) then
      call reset_factor_positions(new_size=2)
      ! H(6, 1)
      factors(1) = eta_H * T0 / rho
      positions(1, :) = [6, 1]
      ! H(6, 5)
      factors(2) = eta_H
      positions(2, :) = [6, 5]
      call subblock(quadblock, factors, positions, current_weight, dh_quad, h_quad)
    end if

    ! ==================== Quadratic * Cubic ====================
    call reset_factor_positions(new_size=2)
    ! H(6, 7)
    factors(1) = eta_H * k3 * Fop_plus / rho
    positions(1, :) = [6, 7]
    ! H(6, 8)
    factors(2) = -eta_H * k2 * Fop_plus / rho
    positions(2, :) = [6, 8]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_cubic)

    ! ==================== Quadratic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! H(6, 7)
    factors(1) = eta_H * B03 * (deps / eps - drho / rho) / rho
    positions(1, :) = [6, 7]
    ! H(6, 8)
    factors(2) = eta_H * eps * B02 * (deps / eps + drho / rho) / rho
    positions(2, :) = [6, 8]
    call subblock(quadblock, factors, positions, current_weight, h_quad, dh_cubic)

    ! ==================== dQuadratic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! H(6, 7)
    factors(1) = eta_H * B03 / rho
    positions(1, :) = [6, 7]
    ! H(6, 8)
    factors(2) = -eta_H * eps * B02 / rho
    positions(2, :) = [6, 8]
    call subblock(quadblock, factors, positions, current_weight, dh_quad, dh_cubic)

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! H(7, 6)
    factors(1) = eta_H * (WVop * B03 - ic * B01 * eps * k3 * drho / rho) / rho
    positions(1, :) = [7, 6]
    ! H(8, 6)
    factors(2) = -eta_H * (WVop * B02 + ic * B01 * k2 * (deps / eps - drho / rho)) / rho
    positions(2, :) = [8, 6]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_quad)
    if (elec_pressure) then
      call reset_factor_positions(new_size=4)
      ! H(7, 1)
      factors(1) = -eta_H * T0 * k2 / (eps * rho)
      positions(1, :) = [7, 1]
      ! H(7, 5)
      factors(2) = -eta_H * k2 / eps
      positions(2, :) = [7, 5]
      ! H(8, 1)
      factors(3) = -eta_H * T0 * k3 / rho
      positions(3, :) = [8, 1]
      ! H(8, 5)
      factors(4) = -eta_H * k3
      positions(4, :) = [8, 5]
      call subblock(quadblock, factors, positions, current_weight, h_cubic, h_quad)
    end if

    ! ==================== dCubic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! H(7, 6)
    factors(1) = eta_H * ic * B01 * eps * k3 / rho
    positions(1, :) = [7, 6]
    ! H(8, 6)
    factors(2) = -eta_H * ic * B01 * k2 / rho
    positions(2, :) = [8, 6]
    call subblock(quadblock, factors, positions, current_weight, dh_cubic, h_quad)

    ! ==================== Cubic * Cubic ====================
    call reset_factor_positions(new_size=4)
    ! H(7, 7)
    factors(1) = eta_H * k3 * (ic * B01 * k2 / eps - drB02 / eps) / rho
    positions(1, :) = [7, 7]
    ! H(7, 8)
    factors(2) = -eta_H * k2 * (ic * B01 * k2 / eps - drB02 / eps) / rho
    positions(2, :) = [7, 8]
    ! H(8, 7)
    factors(3) = eta_H * k3 * (ic * B01 * k3 - dB03) / rho
    positions(3, :) = [8, 7]
    ! H(8, 8)
    factors(4) = eta_H * k2 * (dB03 - ic * B01 * k3) / rho
    positions(4, :) = [8, 8]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_cubic)

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=4)
    ! H(7, 7)
    factors(1) = -eta_H * k2 * B03 / (eps * rho)
    positions(1, :) = [7, 7]
    ! H(7, 8)
    factors(2) = eta_H * eps * (ic * B01 * drho / rho - k3 * B03) / rho
    positions(2, :) = [7, 8]
    ! H(8, 7)
    factors(3) = eta_H * (ic * B01 * (deps / eps - drho / rho) + k2 * B02 / eps) / rho
    positions(3, :) = [8, 7]
    ! H(8, 8)
    factors(4) = eta_H * eps * k3 * B02 / rho
    positions(4, :) = [8, 8]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, dh_cubic)

    ! ==================== dCubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! H(7, 8)
    factors(1) = -eta_H * ic * B01 * eps / rho
    positions(1, :) = [7, 8]
    ! H(8, 7)
    factors(2) = eta_H * ic * B01 / rho
    positions(2, :) = [8, 7]
    call subblock(quadblock, factors, positions, current_weight, dh_cubic, dh_cubic)

  end procedure add_hall_matrix_terms

end submodule smod_hall_matrix
