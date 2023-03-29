submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_hall
  implicit none

contains

  module procedure add_natural_hall_Bterms
    real(dp)  :: eps
    real(dp)  :: rho
    real(dp)  :: eta_e
    real(dp) :: x
    type(matrix_elements_t) :: elements

    if (.not. settings%physics%hall%is_enabled()) return
    if (.not. settings%physics%hall%has_electron_inertia()) return

    x = grid_gauss(grid_idx)
    eps = eps_grid(grid_idx)
    rho = background%density%rho0(grid_gauss(grid_idx))
    eta_e = physics%hall%inertiafactor(x)
    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Cubic * Quadratic ====================
    call elements%add(eta_e * k2 / rho, "a2", "a1", spline1=h_cubic, spline2=h_quad)
    call elements%add( &
      eta_e * k3 * eps / rho, "a3", "a1", spline1=h_cubic, spline2=h_quad &
    )
    ! ==================== Cubic * dCubic ====================
    call elements%add(-eta_e / rho, "a2", "a2", spline1=h_cubic, spline2=dh_cubic)
    call elements%add( &
      -eta_e * eps / rho, "a3", "a3", spline1=h_cubic, spline2=dh_cubic &
    )

    call add_to_quadblock(quadblock, elements, weight, settings%dims)
    call elements%delete()
  end procedure add_natural_hall_Bterms

  module procedure add_natural_hall_terms
    real(dp)  :: eps, deps
    real(dp)  :: rho, T0, B01, B02, B03
    real(dp)  :: eta_H, mu, efrac
    real(dp) :: x
    type(matrix_elements_t) :: elements

    if (.not. settings%physics%hall%is_enabled()) return
    ! Hall by substitution of the momentum equation, only contains viscosity terms
    if (.not. settings%physics%viscosity%is_enabled()) return

    eps = eps_grid(grid_idx)
    deps = d_eps_grid_dr(grid_idx)
    x = grid_gauss(grid_idx)

    rho = background%density%rho0(x)
    T0 = background%temperature%T0(x)
    B01 = background%magnetic%B01(x)
    B02 = background%magnetic%B02(x)
    B03 = background%magnetic%B03(x)

    eta_H = physics%hall%hallfactor(x)
    mu = settings%physics%viscosity%get_viscosity_value()
    efrac = settings%physics%hall%get_electron_fraction()
    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Quadratic * Cubic ====================
    call elements%add( &
      -eta_H * ic * mu * deps / (eps * rho), &
      "a1", &
      "v1", &
      spline1=h_quad, &
      spline2=h_cubic &
    )
    ! ==================== Quadratic * dCubic ====================
    call elements%add( &
      4.0d0 * eta_H * ic * mu / (3.0d0 * rho), &
      "a1", &
      "v1", &
      spline1=h_quad, &
      spline2=dh_cubic &
    )
    ! ==================== Cubic * dQuadratic ====================
    call elements%add( &
      eta_H * ic * mu * eps / rho, "a2", "v2", spline1=h_cubic, spline2=dh_quad &
    )
    call elements%add( &
      eta_H * ic * mu / rho, "a3", "v3", spline1=h_cubic, spline2=dh_quad &
    )
    ! ==================== Cubic * Quadratic ====================
    call elements%add( &
      -eta_H * ic * mu * deps / (eps * rho), &
      "a3", &
      "v3", &
      spline1=h_cubic, &
      spline2=h_quad &
    )

    call add_to_quadblock(quadblock, elements, weight, settings%dims)
    call elements%delete()
  end procedure add_natural_hall_terms

end submodule
