submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_resistive
  implicit none

contains

  module procedure add_natural_resistive_terms
    use mod_equilibrium, only: eta_field

    real(dp)  :: eps, deps
    real(dp)  :: eta
    real(dp)  :: B02, dB02, drB02
    real(dp)  :: B03, dB03
    real(dp) :: gamma_1
    type(matrix_elements_t) :: elements

    if (.not. settings%physics%resistivity%is_enabled()) return

    gamma_1 = settings%physics%get_gamma_1()

    eps = eps_grid(grid_idx)
    deps = d_eps_grid_dr(grid_idx)
    eta = eta_field % eta(grid_idx)
    B02 = B_field % B02(grid_idx)
    dB02 = B_field % d_B02_dr(grid_idx)
    B03 = B_field % B03(grid_idx)
    dB03 = B_field % d_B03_dr(grid_idx)

    drB02 = deps * B02 + eps * dB02
    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Quadratic * Quadratic ====================
    call elements%add( &
      2.0d0 * ic * gamma_1 * eta * (k3 * drB02 - k2 * dB03), &
      "T", &
      "a1", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    ! ==================== Quadratic * dCubic ====================
    call elements%add( &
      2.0d0 * ic * gamma_1 * eta * dB03, "T", "a2", spline1=h_quad, spline2=dh_cubic &
    )
    call elements%add( &
      -2.0d0 * ic * gamma_1 * eta * drB02, "T", "a3", spline1=h_quad, spline2=dh_cubic &
    )
    ! ==================== Cubic * Quadratic ====================
    call elements%add(-ic * eta * k2, "a2", "a1", spline1=h_cubic, spline2=h_quad)
    call elements%add(-ic * eta * eps * k3, "a3", "a1", spline1=h_cubic, spline2=h_quad)
    ! ==================== Cubic * dCubic ====================
    call elements%add(ic * eta, "a2", "a2", spline1=h_cubic, spline2=dh_cubic)
    call elements%add(ic * eta * eps, "a3", "a3", spline1=h_cubic, spline2=dh_cubic)

    call add_to_quadblock(quadblock, elements, weight, settings%dims)
    call elements%delete()
  end procedure add_natural_resistive_terms

end submodule smod_natural_bounds_resistive
