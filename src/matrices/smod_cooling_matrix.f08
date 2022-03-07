submodule (mod_matrix_manager) smod_cooling_matrix
  implicit none

contains

  module procedure add_cooling_matrix_terms
    use mod_equilibrium, only: rc_field
    use mod_global_variables, only: incompressible

    real(dp)  :: rho
    real(dp)  :: Lrho, LT, L0

    if (incompressible) then
      return
    end if

    rho = rho_field % rho0(gauss_idx)
    Lrho = rc_field % d_L_drho(gauss_idx)
    LT = rc_field % d_L_dT(gauss_idx)
    L0 = rc_field % heat_loss(gauss_idx)

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! Lambda(5, 1)
    factors(1) = -ic * gamma_1 * (L0 + rho * Lrho)
    positions(1, :) = [5, 1]
    ! Lambda(5, 5)
    factors(2) = -ic * gamma_1 * rho * LT
    positions(2, :) = [5, 5]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_quad)

  end procedure add_cooling_matrix_terms

end submodule smod_cooling_matrix
