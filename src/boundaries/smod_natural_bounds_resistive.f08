submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_resistive
  implicit none

contains

  module procedure add_natural_resistive_terms
    real(dp)  :: eps, deps
    real(dp)  :: eta
    real(dp)  :: B02, dB02, drB02
    real(dp)  :: B03, dB03
    real(dp) :: gamma_1
    real(dp) :: x
    type(matrix_elements_t) :: elements

    if (.not. settings%physics%resistivity%is_enabled()) return

    gamma_1 = settings%physics%get_gamma_1()

    eps = eps_grid(grid_idx)
    deps = d_eps_grid_dr(grid_idx)
    x = grid_gauss(grid_idx)
    eta = physics%resistivity%eta(x)
    B02 = background%magnetic%B02(x)
    dB02 = background%magnetic%dB02(x)
    B03 = background%magnetic%B03(x)
    dB03 = background%magnetic%dB03(x)

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
