submodule (mod_matrix_manager) smod_conduction_matrix
  use mod_equilibrium, only: kappa_field
  implicit none

contains

  module procedure add_conduction_matrix_terms
    use mod_matrix_shortcuts, only: get_Kp_operator, get_diffF_operator
    use mod_global_variables, only: incompressible

    real(dp)  :: eps, deps
    real(dp)  :: dT0, ddT0
    real(dp)  :: B0, B01, B02, B03
    real(dp)  :: diffKp, Kp, Kp_plus, Kp_plusplus
    real(dp)  :: WVop, Fop_plus, dFop_plus, Gop_min
    real(dp)  :: kappa_para, dkappa_para_dT
    real(dp)  :: kappa_perp, dkappa_perp_drho, dkappa_perp_dT, dkappa_perp_dB2
    complex(dp) :: Fop_B01

    if (incompressible) then
      return
    end if

    ! grid variables
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! temperature variables
    dT0 = T_field % d_T0_dr(gauss_idx)
    ddT0 = T_field % dd_T0_dr(gauss_idx)
    ! magnetic field variables
    B0 = B_field % B0(gauss_idx)
    B01 = B_field % B01
    B02 = B_field % B02(gauss_idx)
    B03 = B_field % B03(gauss_idx)
    ! prefactors
    Kp = kappa_field % prefactor(gauss_idx)
    diffKp = kappa_field % d_prefactor_dr(gauss_idx)
    Kp_plus = get_Kp_operator(gauss_idx, which="+")
    Kp_plusplus = get_Kp_operator(gauss_idx, which="++")
    ! parallel thermal conduction variables
    kappa_para = kappa_field % kappa_para(gauss_idx)
    dkappa_para_dT = kappa_field % d_kappa_para_dT(gauss_idx)
    ! perpendicular thermal conduction variables
    kappa_perp = kappa_field % kappa_perp(gauss_idx)
    dkappa_perp_drho = kappa_field % d_kappa_perp_drho(gauss_idx)
    dkappa_perp_dT = kappa_field % d_kappa_perp_dT(gauss_idx)
    dkappa_perp_dB2 = kappa_field % d_kappa_perp_dB2(gauss_idx)
    ! operators
    WVop = get_wv_operator(gauss_idx)
    Fop_plus = get_F_operator(gauss_idx, which="plus")
    dFop_plus = get_diffF_operator(gauss_idx, which="plus")
    Gop_min = get_G_operator(gauss_idx, which="minus")

    ! B01 modified F+ operator
    Fop_B01 = deps * ic * B01 / eps + Fop_plus

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=3)
    ! K(5, 1)
    factors(1) = gamma_1 * dT0 * (B01 / B0**2) * dkappa_perp_drho * Fop_B01
    positions(1, :) = [5, 1]
    ! K(5, 5)
    factors(2) = gamma_1 * ( &
      B01 * Fop_plus * diffKp &
      - ic * WVop * kappa_perp / eps &
      - B01 * dT0 * (dkappa_para_dT - dkappa_perp_dT) * Fop_B01 / B0**2 &
      + Kp * ( &
        B01 * deps * (2.0d0 * deps * ic * B01 / eps + 3.0d0 * Fop_plus) / eps &
        + B01 * dFop_plus &
        - ic * Fop_plus**2 &
      ) &
    )
    positions(2, :) = [5, 5]
    ! K(5, 6)
    factors(3) = 2.0d0 * gamma_1 * eps * dT0 * Gop_min * Kp_plus * B01 * Fop_B01 / B0**2
    positions(3, :) = [5, 6]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_quad)

    ! ==================== dQuadratic * Quadratic ====================
    call reset_factor_positions(new_size=3)
    ! K(5, 1)
    factors(1) = -ic * gamma_1 * dT0 * dkappa_perp_drho * (1.0d0 - B01**2 / B0**2)
    positions(1, :) = [5, 1]
    ! K(5, 5)
    factors(2) = gamma_1 * ( &
      B01 * Kp * (2.0d0 * deps * ic * B01 / eps + 3.0d0 * Fop_plus) &
      + ic * deps * kappa_perp / eps &
      - ic * dT0 * ( &
        B01**2 * dkappa_para_dT / B0**2 + dkappa_perp_dT * (1.0d0 - B01**2 / B0**2) &
      ) &
    )
    positions(2, :) = [5, 5]
    ! K(5, 6)
    factors(3) = -2.0d0 * ic * gamma_1 * eps * dT0 * Gop_min * Kp_plusplus
    positions(3, :) = [5, 6]
    call subblock(quadblock, factors, positions, current_weight, dh_quad, h_quad)

    ! ==================== Quadratic * dQuadratic ====================
    call reset_factor_positions(new_size=1)
    ! K(5, 5)
    factors(1) = -2.0d0 * ic * gamma_1 * deps * B01**2 * Kp / eps
    positions(1, :) = [5, 5]
    call subblock(quadblock, factors, positions, current_weight, h_quad, dh_quad)

    ! ==================== dQuadratic * dQuadratic ====================
    call reset_factor_positions(new_size=1)
    ! K(5, 5)
    factors(1) = -ic * gamma_1 * (2.0d0 * B01**2 * Kp + kappa_perp)
    positions(1, :) = [5, 5]
    call subblock(quadblock, factors, positions, current_weight, dh_quad, dh_quad)

    ! ==================== Quadratic * Cubic ====================
    call reset_factor_positions(new_size=2)
    ! K(5, 7)
    factors(1) = gamma_1 * k3 * ( &
      Kp * (B01 * ddT0 + ic * dT0 * Fop_plus) &
      + B01 * dT0 * (diffKp - 2.0d0 * ic * B01 * Kp_plus * Fop_B01 / B0**2) &
    )
    positions(1, :) = [5, 7]
    ! K(5, 8)
    factors(2) = -gamma_1 * k2 * ( &
      Kp * (B01 * ddT0 + ic * dT0 * Fop_plus) &
      + B01 * dT0 * (diffKp - 2.0d0 * ic * B01 * Kp_plus * Fop_B01 / B0**2) &
    )
    positions(2, :) = [5, 8]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_cubic)

    ! ==================== Quadratic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! K(5, 7)
    factors(1) = gamma_1 * dT0 * B01 * ( &
      2.0d0 * B03 * Kp_plus * Fop_B01 / B0**2 - k3 * Kp &
    )
    positions(1, :) = [5, 7]
    ! K(5, 8)
    factors(2) = gamma_1 * dT0 * B01 * ( &
      -2.0d0 * eps * B02 * Kp_plus * Fop_B01 / B0**2 + k2 * Kp &
    )
    positions(2, :) = [5, 8]
    call subblock(quadblock, factors, positions, current_weight, h_quad, dh_cubic)

    ! ==================== dQuadratic * Cubic ====================
    call reset_factor_positions(new_size=2)
    ! K(5, 7)
    factors(1) = -2.0d0 * gamma_1 * dT0 * B01 * k3 * Kp_plusplus
    positions(1, :) = [5, 7]
    ! K(5, 8)
    factors(2) = 2.0d0 * gamma_1 * dT0 * B01 * k2 * Kp_plusplus
    positions(2, :) = [5, 8]
    call subblock(quadblock, factors, positions, current_weight, dh_quad, h_cubic)

    ! ==================== dQuadratic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! K(5, 7)
    factors(1) = -2.0d0 * ic * gamma_1 * dT0 * B03 * Kp_plusplus
    positions(1, :) = [5, 7]
    ! K(5, 8)
    factors(2) = 2.0d0 * ic * gamma_1 * dT0 * eps * B02 * Kp_plusplus
    positions(2, :) = [5, 8]
    call subblock(quadblock, factors, positions, current_weight, dh_quad, dh_cubic)
  end procedure add_conduction_matrix_terms

end submodule smod_conduction_matrix
