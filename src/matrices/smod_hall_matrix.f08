submodule (mod_matrix_manager) smod_hall_matrix
  use mod_equilibrium, only: hall_field
  implicit none

contains

  module procedure add_hall_bmatrix_terms
    real(dp)  :: eps, deps
    real(dp)  :: rho, drho
    real(dp)  :: eta_H, eta_e
    real(dp)  :: WVop
    type(matrix_elements_t) :: elements

    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    rho = rho_field % rho0(gauss_idx)
    drho = rho_field % d_rho0_dr(gauss_idx)
    eta_H = hall_field % hallfactor(gauss_idx)
    eta_e = hall_field % inertiafactor(gauss_idx)
    WVop = get_wv_operator(gauss_idx)
    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Quadratic * Cubic ====================
    call elements%add(eta_H, "a1", "v1", spline1=h_quad, spline2=h_cubic)
    ! ==================== Cubic * Quadratic ====================
    call elements%add(eta_H * eps, "a2", "v2", spline1=h_cubic, spline2=h_quad)
    call elements%add(eta_H, "a3", "v3", spline1=h_cubic, spline2=h_quad)

    if (settings%physics%hall%has_electron_inertia()) then
      ! ==================== Quadratic * Quadratic ====================
      call elements%add(eta_e * WVop / rho, "a1", "a1", spline1=h_quad, spline2=h_quad)
      ! ==================== Quadratic * dCubic ====================
      call elements%add( &
        -eta_e * k2 / (eps * rho), "a1", "a2", spline1=h_quad, spline2=dh_cubic &
      )
      call elements%add( &
        -eta_e * eps * k3 / rho, "a1", "a3", spline1=h_quad, spline2=dh_cubic &
      )
      ! ==================== Cubic * Quadratic ====================
      call elements%add( &
        -eta_e * k2 * (deps / (eps * rho) - drho / rho**2), &
        "a2", &
        "a1", &
        spline1=h_cubic, &
        spline2=h_quad &
      )
      call elements%add( &
        eta_e * drho * eps * k3 / rho**2, "a3", "a1", spline1=h_cubic, spline2=h_quad &
      )
      ! ==================== dCubic * Quadratic ====================
      call elements%add(-eta_e * k2 / rho, "a2", "a1", spline1=dh_cubic, spline2=h_quad)
      call elements%add( &
        -eta_e * eps * k3 / rho, "a3", "a1", spline1=dh_cubic, spline2=h_quad &
      )
      ! ==================== Cubic * Cubic ====================
      call elements%add( &
        eta_e * k3**2 / rho, "a2", "a2", spline1=h_cubic, spline2=h_cubic &
      )
      call elements%add( &
        -eta_e * k2 * k3 / rho, "a2", "a3", spline1=h_cubic, spline2=h_cubic &
      )
      call elements%add( &
        -eta_e * k2 * k3 / (eps * rho), "a3", "a2", spline1=h_cubic, spline2=h_cubic &
      )
      call elements%add( &
        eta_e * k2**2 / (eps * rho), "a3", "a3", spline1=h_cubic, spline2=h_cubic &
      )
      ! ==================== Cubic * dCubic ====================
      call elements%add( &
        eta_e * (deps / (eps * rho) - drho / rho**2), &
        "a2", &
        "a2", &
        spline1=h_cubic, &
        spline2=dh_cubic &
      )
      call elements%add( &
        -eta_e * eps * drho / rho**2, "a3", "a3", spline1=h_cubic, spline2=dh_cubic &
      )
      ! ==================== dCubic * dCubic ====================
      call elements%add(eta_e / rho, "a2", "a2", spline1=dh_cubic, spline2=dh_cubic)
      call elements%add( &
        eta_e * eps / rho, "a3", "a3", spline1=dh_cubic, spline2=dh_cubic &
      )
    end if
    call add_to_quadblock(quadblock, elements, current_weight, settings%dims)
    call elements%delete()
  end procedure add_hall_bmatrix_terms


  module procedure add_hall_matrix_terms
    use mod_equilibrium, only: v_field, rho_field, T_field, B_field

    real(dp)  :: eps, deps
    real(dp)  :: v01, v02, v03, dv01, dv02, dv03, ddv01, ddv02, ddv03
    real(dp)  :: rho, drho, T0, dT0, B01, B02, B03, dB02, dB03
    real(dp)  :: drB02, dB02_r, Fop_plus, Gop_plus, Gop_min, WVop
    real(dp)  :: eta_H, mu, efrac
    type(matrix_elements_t) :: elements
    logical :: has_viscosity

    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)

    v01 = v_field % v01(gauss_idx)
    v02 = v_field % v02(gauss_idx)
    v03 = v_field % v03(gauss_idx)
    dv01 = v_field % d_v01_dr(gauss_idx)
    dv02 = v_field % d_v02_dr(gauss_idx)
    dv03 = v_field % d_v03_dr(gauss_idx)
    ddv01 = v_field % dd_v01_dr(gauss_idx)
    ddv02 = v_field % dd_v02_dr(gauss_idx)
    ddv03 = v_field % dd_v03_dr(gauss_idx)

    rho = rho_field % rho0(gauss_idx)
    drho = rho_field % d_rho0_dr(gauss_idx)
    T0 = T_field % T0(gauss_idx)
    dT0 = T_field % d_T0_dr(gauss_idx)

    B01 = B_field % B01
    B02 = B_field % B02(gauss_idx)
    dB02 = B_field % d_B02_dr(gauss_idx)
    B03 = B_field % B03(gauss_idx)
    dB03 = B_field % d_B03_dr(gauss_idx)

    ! Calculate derivatives eps*B02, B02/eps
    drB02 = deps * B02 + eps * dB02
    dB02_r = dB02 / eps - B02 * deps / eps**2

    Fop_plus = get_F_operator(gauss_idx, which="plus")
    Gop_plus = get_G_operator(gauss_idx, which="plus")
    Gop_min = get_G_operator(gauss_idx, which="minus")
    WVop = get_wv_operator(gauss_idx)

    eta_H = hall_field % hallfactor(gauss_idx)
    mu = settings%physics%viscosity%get_viscosity_value()
    efrac = settings%physics%hall%get_electron_fraction()

    elements = new_matrix_elements(state_vector=settings%get_state_vector())
    has_viscosity = settings%physics%viscosity%is_enabled()

    ! Hall by substitution of the momentum equation
    ! ==================== Quadratic * Cubic ====================
    call elements%add( &
      eta_H * (k2 * v02 / eps + k3 * v03), &
      "a1", &
      "v1", &
      spline1=h_quad, &
      spline2=h_cubic &
    )
    ! ==================== Quadratic * Quadratic ====================
    call elements%add( &
      -eta_H * (1.0d0 - efrac) * dT0 / rho, &
      "a1", &
      "rho", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      -2.0d0 * eta_H * deps * v02, "a1", "v2", spline1=h_quad, spline2=h_quad &
    )
    call elements%add( &
      eta_H * (1.0d0 - efrac) * drho / rho, &
      "a1", &
      "T", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    ! ==================== Cubic * Cubic ====================
    call elements%add( &
      -eta_H * (dv02 - v02 * deps / eps), &
      "a2", &
      "v1", &
      spline1=h_cubic, &
      spline2=h_cubic &
    )
    call elements%add(-eta_H * dv03, "a3", "v1", spline1=h_cubic, spline2=h_cubic)
    ! ==================== Cubic * Quadratic ====================
    call elements%add( &
      eta_H * (k2 * v02 + eps * k3 * v03), &
      "a2", &
      "v2", &
      spline1=h_cubic, &
      spline2=h_quad &
    )
    call elements%add( &
      eta_H * (k2 * v02 / eps + k3 * v03), &
      "a3", &
      "v3", &
      spline1=h_cubic, &
      spline2=h_quad &
    )

    ! elements below are only relevant when viscosity is included
    if (.not. has_viscosity) return

    ! ==================== Quadratic * Cubic ====================
    call elements%add( &
      -eta_H * ic * mu * ( &
        (drho / rho + 1.0d0 / eps) * deps / eps + (k2 / eps)**2 + k3**2 &
      ) / rho, &
      "a1", &
      "v1", &
      spline1=h_quad, &
      spline2=h_cubic &
    )
    ! ==================== Quadratic * dCubic ====================
    call elements%add( &
      eta_H * ic * mu * (4.0d0 * drho / rho - deps / eps) / (3.0d0 * rho), &
      "a1", &
      "v1", &
      spline1=h_quad, &
      spline2=dh_cubic &
    )
    ! ==================== dQuadratic * Cubic ====================
    call elements%add( &
      eta_H * ic * mu * deps / (eps * rho), &
      "a1", &
      "v1", &
      spline1=dh_quad, &
      spline2=h_cubic &
    )
    ! ==================== dQuadratic * dCubic ====================
    call elements%add( &
      -4.0d0 * eta_H * ic * mu / (3.0d0 * rho), &
      "a1", &
      "v1", &
      spline1=dh_quad, &
      spline2=dh_cubic &
    )
    ! ==================== Quadratic * Quadratic ====================
    call elements%add( &
      eta_H * mu * 4.0d0 * (ddv01 + deps * (dv01 - v01 / eps) / eps) &
      / (3.0d0 * rho**2), &
      "a1", &
      "rho", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      7.0d0 * eta_H * ic * mu * deps * k2 / (3.0d0 * eps * rho), &
      "a1", &
      "v2", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      eta_H * ic * mu * k3 * deps / (3.0d0 * eps * rho), &
      "a1", &
      "v3", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    ! ==================== Quadratic * dQuadratic ====================
    call elements%add( &
      -eta_H * ic * mu * k2 / (3.0d0 * rho), &
      "a1", &
      "v2", &
      spline1=h_quad, &
      spline2=dh_quad &
    )
    call elements%add( &
      -eta_H * ic * mu * k3 / (3.0d0 * rho), &
      "a1", &
      "v3", &
      spline1=h_quad, &
      spline2=dh_quad &
    )
    ! ==================== Cubic * Cubic ====================
    call elements%add( &
      2.0d0 * eta_H * ic * mu * k2 * deps / (eps**2 * rho), &
      "a2", &
      "v1", &
      spline1=h_cubic, &
      spline2=h_cubic &
    )
    ! ==================== Cubic * dCubic ====================
    call elements%add( &
      eta_H * ic * mu * k2 / (3.0d0 * eps * rho), &
      "a2", &
      "v1", &
      spline1=h_cubic, &
      spline2=dh_cubic &
    )
    call elements%add( &
      eta_H * ic * mu * k3 / (3.0d0 * rho), &
      "a3", &
      "v1", &
      spline1=h_cubic, &
      spline2=dh_cubic &
    )
    ! ==================== Cubic * Quadratic ====================
    call elements%add( &
      -ic * eta_H * mu * (ddv02 + deps * (dv02 - v02 / eps) / eps) / rho**2, &
      "a2", &
      "rho", &
      spline1=h_cubic, &
      spline2=h_quad &
    )
    call elements%add( &
      -eta_H * ic * mu * ( &
        4.0d0 * k2**2 / (3.0d0 * eps) + eps * k3**2 + deps / eps &
      ) / rho, &
      "a2", &
      "v2", &
      spline1=h_cubic, &
      spline2=h_quad &
    )
    call elements%add( &
      -eta_H * ic * mu * k2 * k3 / (3.0d0 * eps * rho), &
      "a2", &
      "v3", &
      spline1=h_cubic, &
      spline2=h_quad &
    )
    call elements%add( &
      -ic * eta_H * mu * (ddv03 + deps * dv03 / eps) / rho**2, &
      "a3", &
      "rho", &
      spline1=h_cubic, &
      spline2=h_quad &
    )
    call elements%add( &
      -eta_H * ic * mu * k2 * k3 / (3.0d0 * rho), &
      "a3", &
      "v2", &
      spline1=h_cubic, &
      spline2=h_quad &
    )
    call elements%add( &
      -eta_H * ic * mu * ( &
        (k2 / eps)**2 + 4.0d0 * k3**2 / 3.0d0 + drho * deps / (eps * rho) &
      ) / rho, &
      "a3", &
      "v3", &
      spline1=h_cubic, &
      spline2=h_quad &
    )
    ! ==================== Cubic * dQuadratic ====================
    call elements%add( &
      eta_H * ic * mu * eps * drho / rho**2, &
      "a2", &
      "v2", &
      spline1=h_cubic, &
      spline2=dh_quad &
    )
    call elements%add( &
      eta_H * ic * mu * drho / rho**2, &
      "a3", &
      "v3", &
      spline1=h_cubic, &
      spline2=dh_quad &
    )
    ! ==================== dCubic * dQuadratic ====================
    call elements%add( &
      -eta_H * ic * mu * eps / rho, "a2", "v2", spline1=dh_cubic, spline2=dh_quad &
    )
    call elements%add( &
      -eta_H * ic * mu / rho, "a3", "v3", spline1=dh_cubic, spline2=dh_quad &
    )
    ! ==================== dCubic * Quadratic ====================
    call elements%add( &
      eta_H * ic * mu * deps / (eps * rho), &
      "a3", &
      "v3", &
      spline1=dh_cubic, &
      spline2=h_quad &
    )

    call add_to_quadblock(quadblock, elements, current_weight, settings%dims)
    call elements%delete()
  end procedure add_hall_matrix_terms

end submodule smod_hall_matrix
