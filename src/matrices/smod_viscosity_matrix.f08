submodule (mod_matrix_manager) smod_viscosity_matrix
  use mod_equilibrium, only: v_field
  implicit none

contains

  module procedure add_viscosity_matrix_terms
    real(dp)  :: eps, deps
    real(dp)  :: v01, dv01, ddv01
    real(dp)  :: v02, dv02
    real(dp)  :: v03, dv03, ddv03
    real(dp)  :: mu
    real(dp)  :: WVop
    real(dp) :: gamma_1
    logical :: viscous_heating, is_compressible
    type(matrix_elements_t) :: elements

    gamma_1 = settings%physics%get_gamma_1()
    viscous_heating = settings%physics%viscosity%has_viscous_heating()
    is_compressible = .not. settings%physics%is_incompressible
    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! grid variables
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! viscous heating variables
    v01 = v_field % v01(gauss_idx)
    dv01 = v_field % d_v01_dr(gauss_idx)
    ddv01 = v_field % dd_v01_dr(gauss_idx)
    v02 = v_field % v02(gauss_idx)
    dv02 = v_field % d_v02_dr(gauss_idx)
    v03 = v_field % v03(gauss_idx)
    dv03 = v_field % d_v03_dr(gauss_idx)
    ddv03 = v_field % dd_v03_dr(gauss_idx)
    ! viscosity value
    mu = settings%physics%viscosity%get_viscosity_value()
    ! operators
    WVop = get_wv_operator(gauss_idx)

    ! ==================== Cubic * Cubic ====================
    call elements%add( &
      -ic * mu * (deps / eps + WVop) / eps, &
      "v1", &
      "v1", &
      spline1=h_cubic, &
      spline2=h_cubic &
    )
    ! ==================== Cubic * dCubic ====================
    call elements%add( &
      -ic * mu * deps / (3.0d0 * eps), "v1", "v1", spline1=h_cubic, spline2=dh_cubic &
    )
    ! ==================== dCubic * Cubic ====================
    call elements%add( &
      ic * mu * deps / eps, "v1", "v1", spline1=dh_cubic, spline2=h_cubic &
    )
    ! ==================== dCubic * dCubic ====================
    call elements%add( &
      -4.0d0 * ic * mu / 3.0d0, "v1", "v1", spline1=dh_cubic, spline2=dh_cubic &
    )
    ! ==================== Cubic * Quadratic ====================
    call elements%add( &
      7.0d0 * deps * ic * mu * k2 / (3.0d0 * eps), &
      "v1", &
      "v2", &
      spline1=h_cubic, &
      spline2=h_quad &
    )
    call elements%add( &
      ic * mu * deps * k3 / (3.0d0 * eps), "v1", "v3", spline1=h_cubic, spline2=h_quad &
    )
    ! ==================== dCubic * Quadratic ====================
    call elements%add( &
      ic * mu * k2 / 3.0d0, "v1", "v2", spline1=dh_cubic, spline2=h_quad &
    )
    call elements%add( &
      ic * mu * k3 / 3.0d0, "v1", "v3", spline1=dh_cubic, spline2=h_quad &
    )
    ! ==================== Quadratic * Cubic ====================
    call elements%add( &
      ic * mu * deps * 2.0d0 * k2 / eps**2, &
      "v2", &
      "v1", &
      spline1=h_quad, &
      spline2=h_cubic &
    )
    ! ==================== Quadratic * dCubic ====================
    call elements%add( &
      ic * mu * k2 / (3.0d0 * eps), "v2", "v1", spline1=h_quad, spline2=dh_cubic &
    )
    call elements%add( &
      ic * mu * k3 / 3.0d0, "v3", "v1", spline1=h_quad, spline2=dh_cubic &
    )
    ! ==================== Quadratic * Quadratic ====================
    call elements%add( &
      -ic * mu * (deps / eps + 4.0d0 * k2**2 / (3.0d0 * eps) + eps * k3**2), &
      "v2", &
      "v2", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      -ic * mu * k2 * k3 / (3.0d0 * eps), &
      "v2", &
      "v3", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      -ic * mu * k2 * k3 / 3.0d0, &
      "v3", &
      "v2", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      -ic * mu * (k2**2 / eps**2 + 4.0d0 * k3**2 / 3.0d0), &
      "v3", &
      "v3", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    ! ==================== dQuadratic * dQuadratic ====================
    call elements%add(-ic * mu * eps, "v2", "v2", spline1=dh_quad, spline2=dh_quad)
    call elements%add(-ic * mu, "v3", "v3", spline1=dh_quad, spline2=dh_quad)
    ! ==================== dQuadratic * Quadratic ====================
    call elements%add(ic * mu * deps / eps, "v3", "v3", spline1=dh_quad, spline2=h_quad)

    if (viscous_heating .and. is_compressible) then
      ! ==================== Quadratic * Cubic ====================
      call elements%add( &
        2.0d0 * gamma_1 * mu * ( &
          (deps**2 * v01 - ic * deps * k2 * v02) / eps**2 - deps * dv01 / eps - ddv01 &
        ), &
        "T", &
        "v1", &
        spline1=h_quad, &
        spline2=h_cubic &
      )
      ! ==================== Quadratic * Quadratic ====================
      call elements%add( &
        2.0d0 * gamma_1 * mu * (deps**2 * ic * v02 - deps * k2 * v01) / eps, &
        "T", &
        "v2", &
        spline1=h_quad, &
        spline2=h_quad &
      )
      call elements%add( &
        -2.0d0 * ic * gamma_1 * mu * (deps * dv03 / eps + ddv03), &
        "T", &
        "v3", &
        spline1=h_quad, &
        spline2=h_quad &
      )
      ! ==================== dQuadratic * Quadratic ====================
      call elements%add( &
        -2.0d0 * ic * gamma_1 * mu * dv03, "T", "v3", spline1=dh_quad, spline2=h_quad &
      )
      ! ==================== dQuadratic * Cubic ====================
      call elements%add( &
        -2.0d0 * gamma_1 * mu * dv01, "T", "v1", spline1=dh_quad, spline2=h_cubic &
      )
      ! ==================== Quadratic * dQuadratic ====================
      call elements%add( &
        2.0d0 * ic * gamma_1 * mu * eps * dv02, &
        "T", &
        "v2", &
        spline1=h_quad, &
        spline2=dh_quad &
      )
    end if

    call add_to_quadblock(quadblock, elements, current_weight, settings%dims)
    call elements%delete()
  end procedure add_viscosity_matrix_terms

end submodule smod_viscosity_matrix
