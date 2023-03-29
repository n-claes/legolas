submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_regular
  implicit none

contains

  module procedure add_natural_regular_terms
    real(dp)  :: eps
    real(dp)  :: rho, T0
    real(dp)  :: B01, B02, B03
    real(dp)  :: Gop_min
    type(matrix_elements_t) :: elements

    eps = eps_grid(grid_idx)
    rho = rho_field % rho0(grid_idx)
    T0 = T_field % T0(grid_idx)
    B01 = B_field % B01
    B02 = B_field % B02(grid_idx)
    B03 = B_field % B03(grid_idx)
    Gop_min = k3 * B02 - k2 * B03 / eps
    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Cubic * Quadratic ====================
    call elements%add(T0, "v1", "rho", spline1=h_cubic, spline2=h_quad)
    call elements%add(rho, "v1", "T", spline1=h_cubic, spline2=h_quad)
    call elements%add(eps * Gop_min, "v1", "a1", spline1=h_cubic, spline2=h_quad)
    ! ==================== Cubic * dCubic ====================
    call elements%add(B03, "v1", "a2", spline1=h_cubic, spline2=dh_cubic)
    call elements%add(-eps * B02, "v1", "a3", spline1=h_cubic, spline2=dh_cubic)
    ! ==================== Quadratic * Quadratic ====================
    call elements%add(ic * eps * k3 * B01, "v2", "a1", spline1=h_quad, spline2=h_quad)
    call elements%add(-ic * k2 * B01, "v3", "a1", spline1=h_quad, spline2=h_quad)
    ! ==================== Quadratic * dCubic ====================
    call elements%add(-ic * eps * B01, "v2", "a3", spline1=h_quad, spline2=dh_cubic)
    call elements%add(ic * B01, "v3", "a2", spline1=h_quad, spline2=dh_cubic)

    call add_to_quadblock(quadblock, elements, weight, settings%dims)
    call elements%delete()
  end procedure add_natural_regular_terms

end submodule smod_natural_bounds_regular
