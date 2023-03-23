submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_viscosity
  implicit none

contains

  module procedure add_natural_viscosity_terms
    real(dp)  :: eps, deps
    real(dp)  :: mu
    real(dp)  :: dv01, dv03
    real(dp) :: gamma_1
    logical :: viscous_heating, is_compressible
    type(matrix_elements_t) :: elements

    if (.not. settings%physics%viscosity%is_enabled()) return

    gamma_1 = settings%physics%get_gamma_1()
    viscous_heating = settings%physics%viscosity%has_viscous_heating()
    is_compressible = .not. settings%physics%is_incompressible
    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    eps = eps_grid(grid_idx)
    deps = d_eps_grid_dr(grid_idx)
    mu = settings%physics%viscosity%get_viscosity_value()
    dv01 = background%velocity%dv01(grid_gauss(grid_idx))
    dv03 = background%velocity%dv03(grid_gauss(grid_idx))

    ! ==================== Cubic * Cubic ====================
    call elements%add( &
      -ic * mu * deps / eps, "v1", "v1", spline1=h_cubic, spline2=h_cubic &
    )
    ! ==================== Cubic * dCubic ====================
    call elements%add( &
      4.0d0 * ic * mu / 3.0d0, "v1", "v1", spline1=h_cubic, spline2=dh_cubic &
    )
    ! ==================== Cubic * Quadratic ====================
    call elements%add( &
      -ic * mu * k2 / 3.0d0, "v1", "v2", spline1=h_cubic, spline2=h_quad &
    )
    call elements%add( &
      -ic * mu * k3 / 3.0d0, "v1", "v3", spline1=h_cubic, spline2=h_quad &
    )
    ! ==================== Quadratic * dQuadratic ====================
    call elements%add(ic * mu * eps, "v2", "v2", spline1=h_quad, spline2=dh_quad)
    call elements%add(ic * mu, "v3", "v3", spline1=h_quad, spline2=dh_quad)
    ! ==================== Quadratic * Quadratic ====================
    call elements%add(-ic * mu * deps / eps, "v3", "v3", spline1=h_quad, spline2=h_quad)

    if (viscous_heating .and. is_compressible) then
      ! ==================== Quadratic * Quadratic ====================
      call elements%add( &
        2.0d0 * ic * gamma_1 * mu * dv03, "T", "v3", spline1=h_quad, spline2=h_quad &
      )
      ! ==================== Quadratic * Cubic ====================
      call elements%add( &
        2.0d0 * gamma_1 * mu * dv01, "T", "v1", spline1=h_quad, spline2=h_cubic &
      )
    end if

    call add_to_quadblock(quadblock, elements, weight, settings%dims)
    call elements%delete()
  end procedure add_natural_viscosity_terms

end submodule
