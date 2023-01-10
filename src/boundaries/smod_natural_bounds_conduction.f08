submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_conduction
  implicit none

contains

  module procedure add_natural_conduction_terms
    use mod_equilibrium, only: kappa_field
    use mod_matrix_shortcuts, only: get_Kp_operator, get_F_operator, get_G_operator

    real(dp)  :: eps, deps
    real(dp)  :: dT0
    real(dp)  :: B0, B01, B02, B03
    real(dp)  :: dkappa_para_dT
    real(dp)  :: kappa_perp
    real(dp)  :: dkappa_perp_drho, dkappa_perp_dT
    real(dp)  :: Fop, Gop_min, Kp, Kp_plusplus
    real(dp) :: gamma_1
    type(matrix_elements_t) :: elements

    if (.not. settings%physics%conduction%is_enabled()) return

    gamma_1 = settings%physics%get_gamma_1()

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

    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Quadratic * Quadratic ====================
    call elements%add( &
      ic * gamma_1 * dT0 * dkappa_perp_drho * (1.0d0 - B01**2 / B0**2), &
      "T", &
      "rho", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      gamma_1 * ( &
        -B01 * Kp * (2.0d0 * (deps / eps) * ic * B01 + 3.0d0 * Fop) &
        - deps * ic * kappa_perp / eps &
        + ic * dT0 * ( &
          B01**2 * dkappa_para_dT / B0**2 + dkappa_perp_dT * (1.0d0 - B01**2 / B0**2) &
        ) &
      ), &
      "T", &
      "T", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      2.0d0 * ic * gamma_1 * eps * dT0 * Gop_min * Kp_plusplus, &
      "T", &
      "a1", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    ! ==================== Quadratic * dQuadratic ====================
    call elements%add( &
      ic * gamma_1 * (2.0d0 * B01**2 * Kp + kappa_perp), &
      "T", &
      "T", &
      spline1=h_quad, &
      spline2=dh_quad &
    )
    ! ==================== Quadratic * Cubic ====================
    call elements%add( &
      2.0d0 * gamma_1 * k3 * dT0 * B01 * Kp_plusplus, &
      "T", &
      "a2", &
      spline1=h_quad, &
      spline2=h_cubic &
    )
    call elements%add( &
      -2.0d0 * gamma_1 * k2 * dT0 * B01 * Kp_plusplus, &
      "T", &
      "a3", &
      spline1=h_quad, &
      spline2=h_cubic &
    )
    ! ==================== Quadratic * dCubic ====================
    call elements%add( &
      2.0d0 * ic * gamma_1 * dT0 * B03 * Kp_plusplus, &
      "T", &
      "a2", &
      spline1=h_quad, &
      spline2=dh_cubic &
    )
    call elements%add( &
      -2.0d0 * ic * gamma_1 * dT0 * eps * B02 * Kp_plusplus, &
      "T", &
      "a3", &
      spline1=h_quad, &
      spline2=dh_cubic &
    )

    call add_to_quadblock(quadblock, elements, weight, settings%dims)
    call elements%delete()
  end procedure add_natural_conduction_terms

end submodule smod_natural_bounds_conduction
