submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_conduction
  implicit none

contains

  module procedure add_natural_conduction_terms
    use mod_global_variables, only: gamma_1, thermal_conduction
    use mod_equilibrium, only: kappa_field
    use mod_matrix_shortcuts, only: get_Kp_operator, get_F_operator, get_G_operator

    real(dp)  :: eps, deps
    real(dp)  :: dT0
    real(dp)  :: B0, B01, B02, B03
    real(dp)  :: dkappa_para_dT
    real(dp)  :: kappa_perp
    real(dp)  :: dkappa_perp_drho, dkappa_perp_dT
    real(dp)  :: Fop, Gop_min, Kp, Kp_plusplus

    if (.not. thermal_conduction) then
      return
    end if

    eps = eps_grid(grid_idx)
    deps = d_eps_grid_dr(grid_idx)
    dT0 = T_field % d_T0_dr(grid_idx)
    B0 = B_field % B0(grid_idx)
    B01 = B_field % B01
    B02 = B_field % B02(grid_idx)
    B03 = B_field % B03(grid_idx)
    dkappa_para_dT = kappa_field % d_kappa_para_dT(grid_idx)
    kappa_perp = kappa_field % kappa_perp(grid_idx)
    dkappa_perp_drho = kappa_field % d_kappa_perp_drho(grid_idx)
    dkappa_perp_dT = kappa_field % d_kappa_perp_dT(grid_idx)

    Gop_min = get_G_operator(grid_idx, which="minus")
    Fop = get_F_operator(grid_idx, which="plus")
    Kp = kappa_field % prefactor(grid_idx)
    Kp_plusplus = get_Kp_operator(grid_idx, which="++")

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=3)
    ! K(5, 1)
    factors(1) = ic * gamma_1 * dT0 * dkappa_perp_drho * (1.0d0 - B01**2 / B0**2)
    positions(1, :) = [5, 1]
    ! K(5, 5)
    factors(2) = gamma_1 * ( &
      -B01 * Kp * (2.0d0 * (deps / eps) * ic * B01 + 3.0d0 * Fop) &
      - deps * ic * kappa_perp / eps &
      + ic * dT0 * ( &
        B01**2 * dkappa_para_dT / B0**2 + dkappa_perp_dT * (1.0d0 - B01**2 / B0**2) &
      ) &
    )
    positions(2, :) = [5, 5]
    ! K(5, 6)
    factors(3) = 2.0d0 * ic * gamma_1 * eps * dT0 * Gop_min * Kp_plusplus
    positions(3, :) = [5, 6]
    call subblock(quadblock, factors, positions, weight, h_quad, h_quad)

    ! ==================== Quadratic * dQuadratic ====================
    call reset_factor_positions(new_size=1)
    ! K(5, 5)
    factors(1) = ic * gamma_1 * (2.0d0 * B01**2 * Kp + kappa_perp)
    positions(1, :) = [5, 5]
    call subblock(quadblock, factors, positions, weight, h_quad, dh_quad)

    ! ==================== Quadratic * Cubic ====================
    call reset_factor_positions(new_size=2)
    ! K(5, 7)
    factors(1) = 2.0d0 * gamma_1 * k3 * dT0 * B01 * Kp_plusplus
    positions(1, :) = [5, 7]
    ! K(5, 8)
    factors(2) = -2.0d0 * gamma_1 * k2 * dT0 * B01 * Kp_plusplus
    positions(2, :) = [5, 8]
    call subblock(quadblock, factors, positions, weight, h_quad, h_cubic)

    ! ==================== Quadratic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! K(5, 7)
    factors(1) = 2.0d0 * ic * gamma_1 * dT0 * B03 * Kp_plusplus
    positions(1, :) = [5, 7]
    ! K(5, 8)
    factors(2) = -2.0d0 * ic * gamma_1 * dT0 * eps * B02 * Kp_plusplus
    positions(2, :) = [5, 8]
    call subblock(quadblock, factors, positions, weight, h_quad, dh_cubic)
  end procedure add_natural_conduction_terms

end submodule smod_natural_bounds_conduction
