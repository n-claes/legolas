submodule (mod_matrix_manager) smod_cooling_matrix
  implicit none

contains

  module procedure add_cooling_matrix_terms
    use mod_equilibrium, only: rc_field

    real(dp)  :: rho
    real(dp)  :: Lrho, LT, L0
    real(dp) :: gamma_1
    type(matrix_elements_t) :: elements

    if (settings%physics%is_incompressible) return

    gamma_1 = settings%physics%get_gamma_1()
    rho = rho_field % rho0(gauss_idx)
    Lrho = rc_field % d_L_drho(gauss_idx)
    LT = rc_field % d_L_dT(gauss_idx)
    L0 = rc_field % heat_loss(gauss_idx)

    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    call elements%add(-ic * gamma_1 * (L0 + rho * Lrho), "T", "rho", h_quad, h_quad)
    call elements%add(-ic * gamma_1 * rho * LT, "T", "T", h_quad, h_quad)

    call add_to_quadblock(quadblock, elements, current_weight, settings%dims)
    call elements%delete()
  end procedure add_cooling_matrix_terms

end submodule smod_cooling_matrix
