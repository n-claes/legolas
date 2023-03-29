submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_conduction
  implicit none

contains

  module procedure add_natural_conduction_terms
    real(dp) :: eps, deps
    real(dp) :: dT0
    real(dp) :: dkappa_para_dT
    real(dp) :: kappa_perp
    real(dp) :: dkappa_perp_drho, dkappa_perp_dT
    real(dp) :: gamma_1
    type(matrix_elements_t) :: elements

    if (.not. settings%physics%conduction%is_enabled()) return


    gamma_1 = settings%physics%get_gamma_1()
    eps = grid%get_eps(x)
    deps = grid%get_deps()
    dT0 = background%temperature%dT0(x)
    dkappa_para_dT = physics%conduction%dtcparadT(x)
    kappa_perp = physics%conduction%tcperp(x)
    dkappa_perp_drho = physics%conduction%dtcperpdrho(x)
    dkappa_perp_dT = physics%conduction%dtcperpdT(x)

    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Quadratic * Quadratic ====================
    call elements%add( &
      ic * gamma_1 * dT0 * dkappa_perp_drho, &
      "T", &
      "rho", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      gamma_1 * (-deps * ic * kappa_perp / eps + ic * dT0 * dkappa_perp_dT), &
      "T", &
      "T", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    ! ==================== Quadratic * dQuadratic ====================
    call elements%add( &
      ic * gamma_1 * kappa_perp, "T", "T", spline1=h_quad, spline2=dh_quad &
    )

    if (settings%has_bfield()) then
      call add_natural_conduction_terms_bfield( &
        x, settings, grid, background, physics, elements &
      )
    end if

    call add_to_quadblock(quadblock, elements, weight, settings%dims)
    call elements%delete()
  end procedure add_natural_conduction_terms


  subroutine add_natural_conduction_terms_bfield( &
    x, settings, grid, background, physics, elements &
  )
    real(dp), intent(in) :: x
    type(settings_t), intent(in) :: settings
    type(grid_t), intent(in) :: grid
    type(background_t), intent(in) :: background
    type(physics_t), intent(in) :: physics
    type(matrix_elements_t), intent(inout) :: elements

    real(dp) :: eps, deps
    real(dp) :: dT0
    real(dp) :: B0, B01, B02, B03
    real(dp) :: dkappa_para_dT
    real(dp) :: kappa_perp
    real(dp) :: dkappa_perp_drho, dkappa_perp_dT, dkappa_perp_dB2
    real(dp) :: Fop, Gop_min, Kp, Kp_plus, Kp_plusplus
    real(dp) :: gamma_1

    gamma_1 = settings%physics%get_gamma_1()
    eps = grid%get_eps(x)
    deps = grid%get_deps()
    dT0 = background%temperature%dT0(x)
    dkappa_para_dT = physics%conduction%dtcparadT(x)
    kappa_perp = physics%conduction%tcperp(x)
    dkappa_perp_drho = physics%conduction%dtcperpdrho(x)
    dkappa_perp_dT = physics%conduction%dtcperpdT(x)
    B0 = background%magnetic%get_B0(x)
    B01 = background%magnetic%B01(x)
    B02 = background%magnetic%B02(x)
    B03 = background%magnetic%B03(x)
    dkappa_perp_dB2 = physics%conduction%dtcperpdB2(x)

    Gop_min = k3 * B02 - k2 * B03 / eps
    Fop = k2 * B02 / eps + k3 * B03
    Kp = physics%conduction%tcprefactor(x)
    Kp_plus = Kp + dkappa_perp_dB2
    Kp_plusplus = dkappa_perp_dB2 - (B01**2 * Kp_plus / B0**2)

    ! ==================== Quadratic * Quadratic ====================
    call elements%add( &
      -ic * gamma_1 * dT0 * dkappa_perp_drho * B01**2 / B0**2, &
      "T", &
      "rho", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      gamma_1 * ( &
        -B01 * Kp * (2.0d0 * (deps / eps) * ic * B01 + 3.0d0 * Fop) &
        + ic * dT0 * ( &
          B01**2 * dkappa_para_dT / B0**2 - dkappa_perp_dT * B01**2 / B0**2 &
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
      ic * gamma_1 * 2.0d0 * B01**2 * Kp, "T", "T", spline1=h_quad, spline2=dh_quad &
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
  end subroutine add_natural_conduction_terms_bfield

end submodule smod_natural_bounds_conduction
