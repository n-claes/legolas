submodule (mod_matrix_manager) smod_hall_matrix
  use mod_equilibrium, only: hall_field
  implicit none

contains

  module procedure add_hall_bmatrix_terms
    real(dp)  :: eps, deps
    real(dp)  :: rho, drho
    real(dp)  :: eta_H, eta_e
    real(dp)  :: WVop

    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    rho = rho_field % rho0(gauss_idx)
    drho = rho_field % d_rho0_dr(gauss_idx)
    eta_H = hall_field % hallfactor(gauss_idx)
    eta_e = hall_field % inertiafactor(gauss_idx)
    WVop = get_wv_operator(gauss_idx)

    if (settings%physics%hall%is_using_substitution()) then
      ! ==================== Quadratic * Cubic ====================
      call reset_factor_positions(new_size=1)
      ! B_H(6, 2)
      factors(1) = eta_H
      positions(1, :) = [6, 2]
      call subblock( &
        quadblock, factors, positions, current_weight, h_quad, h_cubic, settings%dims &
      )

      ! ==================== Cubic * Quadratic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(7, 3)
      factors(1) = eta_H * eps
      positions(1, :) = [7, 3]
      ! B_H(8, 4)
      factors(2) = eta_H
      positions(2, :) = [8, 4]
      call subblock( &
        quadblock, factors, positions, current_weight, h_cubic, h_quad, settings%dims &
      )
    end if

    if (settings%physics%hall%has_electron_inertia()) then
      ! ==================== Quadratic * Quadratic ====================
      call reset_factor_positions(new_size=1)
      ! B_H(6, 6)
      factors(1) = eta_e * WVop / rho
      positions(1, :) = [6, 6]
      call subblock( &
        quadblock, factors, positions, current_weight, h_quad, h_quad, settings%dims &
      )

      ! ==================== Quadratic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(6, 7)
      factors(1) = -eta_e * k2 / (eps * rho)
      positions(1, :) = [6, 7]
      ! B_H(6, 8)
      factors(2) = -eta_e * eps * k3 / rho
      positions(2, :) = [6, 8]
      call subblock( &
        quadblock, factors, positions, current_weight, h_quad, dh_cubic, settings%dims &
      )

      ! ==================== Cubic * Quadratic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(7, 6)
      factors(1) = -eta_e * k2 * (deps / (eps * rho) - drho / rho**2)
      positions(1, :) = [7, 6]
      ! B_H(8, 6)
      factors(2) = eta_e * drho * eps * k3 / rho**2
      positions(2, :) = [8, 6]
      call subblock( &
        quadblock, factors, positions, current_weight, h_cubic, h_quad, settings%dims &
      )

      ! ==================== dCubic * Quadratic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(7, 6)
      factors(1) = -eta_e * k2 / rho
      positions(1, :) = [7, 6]
      ! B_H(8, 6)
      factors(2) = -eta_e * eps * k3 / rho
      positions(2, :) = [8, 6]
      call subblock( &
        quadblock, factors, positions, current_weight, dh_cubic, h_quad, settings%dims &
      )

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
      call subblock( &
        quadblock, factors, positions, current_weight, h_cubic, h_cubic, settings%dims &
      )

      ! ==================== Cubic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(7, 7)
      factors(1) = eta_e * (deps / (eps * rho) - drho / rho**2)
      positions(1, :) = [7, 7]
      ! B_H(8, 8)
      factors(2) = -eta_e * eps * drho / rho**2
      positions(2, :) = [8, 8]
      call subblock( &
        quadblock, &
        factors, &
        positions, &
        current_weight, &
        h_cubic, &
        dh_cubic, &
        settings%dims &
      )

      ! ==================== dCubic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! B_H(7, 7)
      factors(1) = eta_e / rho
      positions(1, :) = [7, 7]
      ! B_H(8, 8)
      factors(2) = eta_e * eps / rho
      positions(2, :) = [8, 8]
      call subblock( &
        quadblock, &
        factors, &
        positions, &
        current_weight, &
        dh_cubic, &
        dh_cubic, &
        settings%dims &
      )
    end if

  end procedure add_hall_bmatrix_terms


  module procedure add_hall_matrix_terms
    use mod_equilibrium, only: v_field, rho_field, T_field, B_field

    real(dp)  :: eps, deps
    real(dp)  :: v01, v02, v03, dv01, dv02, dv03, ddv01, ddv02, ddv03
    real(dp)  :: rho, drho, T0, dT0, B01, B02, B03, dB02, dB03
    real(dp)  :: drB02, dB02_r, Fop_plus, Gop_plus, Gop_min, WVop
    real(dp)  :: eta_H, mu, efrac

    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)

    v01 = v_field % v01(gauss_idx)
    v02 = v_field % v02(gauss_idx)
    v03 = v_field % v03(gauss_idx)
    dv01 = v_field % d_v01_dr(gauss_idx)
    dv02 = v_field % d_v02_dr(gauss_idx)
    dv03 = v_field % d_v03_dr(gauss_idx)
    ddv01 = v_field % dd_v01_dr(gauss_idx)
    ddv02 = v_field % dd_v02_dr(gauss_idx)
    ddv03 = v_field % dd_v03_dr(gauss_idx)

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

    Fop_plus = get_F_operator(gauss_idx, which="plus")
    Gop_plus = get_G_operator(gauss_idx, which="plus")
    Gop_min = get_G_operator(gauss_idx, which="minus")
    WVop = get_wv_operator(gauss_idx)

    eta_H = hall_field % hallfactor(gauss_idx)
    mu = 0.0_dp
    efrac = settings%physics%hall%get_electron_fraction()

    ! Hall by substitution of the momentum equation
    if (settings%physics%hall%is_using_substitution()) then
      ! ==================== Quadratic * Cubic ====================
      call reset_factor_positions(new_size=1)
      ! H(6, 2)
      factors(1) = eta_H * (k2 * v02 / eps + k3 * v03)
      positions(1, :) = [6, 2]
      call subblock( &
        quadblock, factors, positions, current_weight, h_quad, h_cubic, settings%dims &
      )

      ! ==================== Quadratic * Quadratic ====================
      call reset_factor_positions(new_size=3)
      ! H(6, 1)
      factors(1) = -eta_H * (1.0d0-efrac) * dT0 / rho
      positions(1, :) = [6, 1]
      ! H(6, 3)
      factors(2) = -2.0d0 * eta_H * deps * v02
      positions(2, :) = [6, 3]
      ! H(6, 5)
      factors(3) = eta_H * (1.0d0-efrac) * drho / rho
      positions(3, :) = [6, 5]
      call subblock( &
        quadblock, factors, positions, current_weight, h_quad, h_quad, settings%dims &
      )

      ! ==================== Cubic * Cubic ====================
      call reset_factor_positions(new_size=2)
      ! H(7, 2)
      factors(1) = -eta_H * (dv02 - v02 * deps / eps)
      positions(1, :) = [7, 2]
      ! H(8, 2)
      factors(2) = -eta_H * dv03
      positions(2, :) = [8, 2]
      call subblock( &
        quadblock, factors, positions, current_weight, h_cubic, h_cubic, settings%dims &
      )

      ! ==================== Cubic * Quadratic ====================
      call reset_factor_positions(new_size=2)
      ! H(7, 3)
      factors(1) = eta_H * (k2 * v02 + eps * k3 * v03)
      positions(1, :) = [7, 3]
      ! H(8, 4)
      factors(2) = eta_H * (k2 * v02 / eps + k3 * v03)
      positions(2, :) = [8, 4]
      call subblock( &
        quadblock, factors, positions, current_weight, h_cubic, h_quad, settings%dims &
      )

      if (settings%physics%viscosity%is_enabled()) then
        mu = settings%physics%viscosity%get_viscosity_value()
        ! ==================== Quadratic * Cubic ====================
        call reset_factor_positions(new_size=1)
        ! H(6, 2)
        factors(1) = -eta_H * ic * mu * ((drho / rho + 1.0d0 / eps) * deps / eps &
                      + (k2 / eps)**2 + k3**2) / rho
        positions(1, :) = [6, 2]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          h_quad, &
          h_cubic, &
          settings%dims &
        )

        ! ==================== Quadratic * dCubic ====================
        call reset_factor_positions(new_size=1)
        ! H(6, 2)
        factors(1) = eta_H * ic * mu * (4.0d0 * drho / rho - deps / eps) / (3.0d0 * rho)
        positions(1, :) = [6, 2]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          h_quad, &
          dh_cubic, &
          settings%dims &
        )

        ! ==================== dQuadratic * Cubic ====================
        call reset_factor_positions(new_size=1)
        ! H(6, 2)
        factors(1) = eta_H * ic * mu * deps / (eps * rho)
        positions(1, :) = [6, 2]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          dh_quad, &
          h_cubic, &
          settings%dims &
        )

        ! ==================== dQuadratic * dCubic ====================
        call reset_factor_positions(new_size=1)
        ! H(6, 2)
        factors(1) = -4.0d0 * eta_H * ic * mu / (3.0d0 * rho)
        positions(1, :) = [6, 2]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          dh_quad, &
          dh_cubic, &
          settings%dims &
        )

        ! ==================== Quadratic * Quadratic ====================
        call reset_factor_positions(new_size=3)
        ! H(6, 1)
        factors(1) = eta_H * mu * 4.0d0 * (ddv01 + deps * (dv01 - v01 / eps) / eps) / (3.0d0 * rho**2)
        positions(1, :) = [6, 1]
        ! H(6, 3)
        factors(2) = 7.0d0 * eta_H * ic * mu * deps * k2 / (3.0d0 * eps * rho)
        positions(2, :) = [6, 3]
        ! H(6, 4)
        factors(3) = eta_H * ic * mu * k3 * deps / (3.0d0 * eps * rho)
        positions(3, :) = [6, 4]
        call subblock( &
          quadblock, factors, positions, current_weight, h_quad, h_quad, settings%dims &
        )

        ! ==================== Quadratic * dQuadratic ====================
        call reset_factor_positions(new_size=2)
        ! H(6, 3)
        factors(1) = -eta_H * ic * mu * k2 / (3.0d0 * rho)
        positions(1, :) = [6, 3]
        ! H(6, 4)
        factors(2) = -eta_H * ic * mu * k3 / (3.0d0 * rho)
        positions(2, :) = [6, 4]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          h_quad, &
          dh_quad, &
          settings%dims &
        )

        ! ==================== Cubic * Cubic ====================
        call reset_factor_positions(new_size=1)
        ! H(7, 2)
        factors(1) = 2.0d0 * eta_H * ic * mu * k2 * deps / (eps**2 * rho)
        positions(1, :) = [7, 2]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          h_cubic, &
          h_cubic, &
          settings%dims &
        )

        ! ==================== Cubic * dCubic ====================
        call reset_factor_positions(new_size=2)
        ! H(7, 2)
        factors(1) = eta_H * ic * mu * k2 / (3.0d0 * eps * rho)
        positions(1, :) = [7, 2]
        ! H(8, 2)
        factors(2) = eta_H * ic * mu * k3 / (3.0d0 * rho)
        positions(2, :) = [8, 2]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          h_cubic, &
          dh_cubic, &
          settings%dims &
        )

        ! ==================== Cubic * Quadratic ====================
        call reset_factor_positions(new_size=6)
        ! H(7, 1)
        factors(1) = -ic * eta_H * mu * (ddv02 + deps * (dv02 - v02 / eps) / eps) / rho**2
        positions(1, :) = [7, 1]
        ! H(7, 3)
        factors(2) = -eta_H * ic * mu * (4.0d0 * k2**2 / (3.0d0 * eps) &
                      + eps * k3**2 + deps / eps) / rho
        positions(2, :) = [7, 3]
        ! H(7, 4)
        factors(3) = -eta_H * ic * mu * k2 * k3 / (3.0d0 * eps * rho)
        positions(3, :) = [7, 4]
        ! H(8, 1)
        factors(4) = -ic * eta_H * mu * (ddv03 + deps * dv03 / eps) / rho**2
        positions(4, :) = [8, 1]
        ! H(8, 3)
        factors(5) = -eta_H * ic * mu * k2 * k3 / (3.0d0 * rho)
        positions(5, :) = [8, 3]
        ! H(8, 4)
        factors(6) = -eta_H * ic * mu * ((k2 / eps)**2 + 4.0d0 * k3**2 / 3.0d0 &
                      + drho * deps / (eps * rho)) / rho
        positions(6, :) = [8, 4]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          h_cubic, &
          h_quad, &
          settings%dims &
        )

        ! ==================== Cubic * dQuadratic ====================
        call reset_factor_positions(new_size=2)
        ! H(7, 3)
        factors(1) = eta_H * ic * mu * eps * drho / rho**2
        positions(1, :) = [7, 3]
        ! H(8, 4)
        factors(2) = eta_H * ic * mu * drho / rho**2
        positions(2, :) = [8, 4]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          h_cubic, &
          dh_quad, &
          settings%dims &
        )

        ! ==================== dCubic * dQuadratic ====================
        call reset_factor_positions(new_size=2)
        ! H(7, 3)
        factors(1) = -eta_H * ic * mu * eps / rho
        positions(1, :) = [7, 3]
        ! H(8, 4)
        factors(2) = -eta_H * ic * mu / rho
        positions(2, :) = [8, 4]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          dh_cubic, &
          dh_quad, &
          settings%dims &
        )

        ! ==================== dCubic * Quadratic ====================
        call reset_factor_positions(new_size=1)
        ! H(8, 4)
        factors(1) = eta_H * ic * mu * deps / (eps * rho)
        positions(1, :) = [8, 4]
        call subblock( &
          quadblock, &
          factors, &
          positions, &
          current_weight, &
          dh_cubic, &
          h_quad, &
          settings%dims &
        )
      end if

    ! Hall without substitution, only E redefinition up to a gradient
    else
      ! ==================== Quadratic * Quadratic ====================
      call reset_factor_positions(new_size=3)
      ! H(6, 1)
      factors(1) = eta_H * ((B02 * drB02 / eps + B03 * dB03) / rho + efrac * dT0) / rho
      positions(1, :) = [6, 1]
      ! H(6, 5)
      factors(2) = -eta_H * efrac * drho / rho
      positions(2, :) = [6, 5]
      ! H(6, 6)
      factors(3) = -eta_H * (eps * drho * Gop_min / rho + deps * Gop_plus) / rho
      positions(3, :) = [6, 6]
      call subblock( &
        quadblock, factors, positions, current_weight, h_quad, h_quad, settings%dims &
      )

      ! ==================== dQuadratic * Quadratic ====================
      call reset_factor_positions(new_size=1)
      ! H(6, 6)
      factors(1) = eta_H * eps * Gop_min / rho
      positions(1, :) = [6, 6]
      call subblock( &
        quadblock, factors, positions, current_weight, dh_quad, h_quad, settings%dims &
      )

      ! ==================== Quadratic * Cubic ====================
      call reset_factor_positions(new_size=2)
      ! H(6, 7)
      factors(1) = eta_H * k3 * Fop_plus / rho
      positions(1, :) = [6, 7]
      ! H(6, 8)
      factors(2) = -eta_H * k2 * Fop_plus / rho
      positions(2, :) = [6, 8]
      call subblock( &
        quadblock, factors, positions, current_weight, h_quad, h_cubic, settings%dims &
      )

      ! ==================== Quadratic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! H(6, 7)
      factors(1) = eta_H * B03 * (deps / eps - drho / rho) / rho
      positions(1, :) = [6, 7]
      ! H(6, 8)
      factors(2) = eta_H * B02 * (deps + eps * drho / rho) / rho
      positions(2, :) = [6, 8]
      call subblock( &
        quadblock, factors, positions, current_weight, h_quad, dh_cubic, settings%dims &
      )

      ! ==================== dQuadratic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! H(6, 7)
      factors(1) = eta_H * B03 / rho
      positions(1, :) = [6, 7]
      ! H(6, 8)
      factors(2) = -eta_H * eps * B02 / rho
      positions(2, :) = [6, 8]
      call subblock( &
        quadblock, &
        factors, &
        positions, &
        current_weight, &
        dh_quad, &
        dh_cubic, &
        settings%dims &
      )

      ! ==================== Cubic * Quadratic ====================
      call reset_factor_positions(new_size=4)
      ! H(7, 1)
      factors(1) = eta_H * ic * B01 * drB02 / rho**2
      positions(1, :) = [7, 1]
      ! H(7, 6)
      factors(2) = eta_H * (WVop * B03 - ic * B01 * eps * k3 * drho / rho) / rho
      positions(2, :) = [7, 6]
      ! H(8, 1)
      factors(3) = eta_H * ic * eps * B01 * dB03 / rho**2
      positions(3, :) = [8, 1]
      ! H(8, 6)
      factors(4) = -eta_H * (WVop * B02 + ic * B01 * k2 * (deps / eps - drho / rho)) / rho
      positions(4, :) = [8, 6]
      call subblock( &
        quadblock, factors, positions, current_weight, h_cubic, h_quad, settings%dims &
      )

      ! ==================== dCubic * Quadratic ====================
      call reset_factor_positions(new_size=2)
      ! H(7, 6)
      factors(1) = eta_H * ic * B01 * eps * k3 / rho
      positions(1, :) = [7, 6]
      ! H(8, 6)
      factors(2) = -eta_H * ic * B01 * k2 / rho
      positions(2, :) = [8, 6]
      call subblock( &
        quadblock, factors, positions, current_weight, dh_cubic, h_quad, settings%dims &
      )

      ! ==================== Cubic * Cubic ====================
      call reset_factor_positions(new_size=4)
      ! H(7, 7)
      factors(1) = eta_H * k3 * (ic * B01 * k2 - drB02) / (eps * rho)
      positions(1, :) = [7, 7]
      ! H(7, 8)
      factors(2) = -eta_H * k2 * (ic * B01 * k2 - drB02) / (eps * rho)
      positions(2, :) = [7, 8]
      ! H(8, 7)
      factors(3) = eta_H * k3 * (ic * B01 * k3 - dB03) / rho
      positions(3, :) = [8, 7]
      ! H(8, 8)
      factors(4) = -eta_H * k2 * (ic * B01 * k3 - dB03) / rho
      positions(4, :) = [8, 8]
      call subblock( &
        quadblock, factors, positions, current_weight, h_cubic, h_cubic, settings%dims &
      )

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
      call subblock( &
        quadblock, &
        factors, &
        positions, &
        current_weight, &
        h_cubic, &
        dh_cubic, &
        settings%dims &
      )

      ! ==================== dCubic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! H(7, 8)
      factors(1) = -eta_H * ic * B01 * eps / rho
      positions(1, :) = [7, 8]
      ! H(8, 7)
      factors(2) = eta_H * ic * B01 / rho
      positions(2, :) = [8, 7]
      call subblock( &
        quadblock, &
        factors, &
        positions, &
        current_weight, &
        dh_cubic, &
        dh_cubic, &
        settings%dims &
      )
    end if

  end procedure add_hall_matrix_terms

end submodule smod_hall_matrix
