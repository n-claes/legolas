submodule (mod_matrix_manager) smod_conduction_matrix
  implicit none

contains

  module procedure add_conduction_matrix_terms
    real(dp) :: eps, deps
    real(dp) :: dT0
    real(dp) :: WVop
    real(dp) :: kappa_perp, dkappa_perp_drho, dkappa_perp_dT
    real(dp) :: gamma_1
    type(matrix_elements_t) :: elements

    if (settings%physics%is_incompressible) return

    gamma_1 = settings%physics%get_gamma_1()
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    dT0 = background%temperature%dT0(x_gauss)
    kappa_perp = physics%conduction%tcperp(x_gauss)
    dkappa_perp_drho = physics%conduction%dtcperpdrho(x_gauss)
    dkappa_perp_dT = physics%conduction%dtcperpdT(x_gauss)
    ! operators
    WVop = k2**2 / eps + eps * k3**2

    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Quadratic * Quadratic ====================
    call elements%add( &
      -gamma_1 * ic * WVop * kappa_perp / eps, &
      "T", &
      "T", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    ! ==================== dQuadratic * Quadratic ====================
    call elements%add( &
      -ic * gamma_1 * dT0 * dkappa_perp_drho, &
      "T", &
      "rho", &
      spline1=dh_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      gamma_1 * (ic * deps * kappa_perp / eps - ic * dT0 * dkappa_perp_dT), &
      "T", &
      "T", &
      spline1=dh_quad, &
      spline2=h_quad &
    )
    ! ==================== dQuadratic * dQuadratic ====================
    call elements%add( &
      -ic * gamma_1 * kappa_perp, "T", "T", spline1=dh_quad, spline2=dh_quad &
    )

    if (settings%has_bfield()) then
      call add_conduction_matrix_terms_bfield( &
        gauss_idx, x_gauss, settings, background, physics, elements &
      )
    end if

    call add_to_quadblock(quadblock, elements, weight, settings%dims)
    call elements%delete()
  end procedure add_conduction_matrix_terms


  subroutine add_conduction_matrix_terms_bfield( &
    gauss_idx, x_gauss, settings, background, physics, elements &
  )
    integer, intent(in) :: gauss_idx
    real(dp), intent(in) :: x_gauss
    type(settings_t), intent(in) :: settings
    type(background_t), intent(in) :: background
    type(physics_t), intent(in) :: physics
    type(matrix_elements_t), intent(inout) :: elements

    real(dp)  :: eps, deps
    real(dp)  :: dT0, ddT0
    real(dp)  :: B0, B01, B02, B03, dB02, dB03
    real(dp)  :: diffKp, Kp, Kp_plus, Kp_plusplus
    real(dp)  :: Fop_plus, dFop_plus, Gop_min
    real(dp)  :: kappa_para, dkappa_para_dT
    real(dp)  :: kappa_perp, dkappa_perp_drho, dkappa_perp_dT, dkappa_perp_dB2
    real(dp) :: gamma_1
    complex(dp) :: Fop_B01

    gamma_1 = settings%physics%get_gamma_1()
    ! grid variables
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! temperature variables
    dT0 = background%temperature%dT0(x_gauss)
    ddT0 = background%temperature%ddT0(x_gauss)
    ! magnetic field variables
    B0 = background%magnetic%get_B0(x_gauss)
    B01 = background%magnetic%B01(x_gauss)
    B02 = background%magnetic%B02(x_gauss)
    dB02 = background%magnetic%dB02(x_gauss)
    B03 = background%magnetic%B03(x_gauss)
    dB03 = background%magnetic%dB03(x_gauss)
    ! parallel thermal conduction variables
    kappa_para = physics%conduction%tcpara(x_gauss)
    dkappa_para_dT = physics%conduction%dtcparadT(x_gauss)
    ! perpendicular thermal conduction variables
    kappa_perp = physics%conduction%tcperp(x_gauss)
    dkappa_perp_drho = physics%conduction%dtcperpdrho(x_gauss)
    dkappa_perp_dT = physics%conduction%dtcperpdT(x_gauss)
    dkappa_perp_dB2 = physics%conduction%dtcperpdB2(x_gauss)
    ! prefactors
    Kp = physics%conduction%tcprefactor(x_gauss)
    diffKp = physics%conduction%dtcprefactordr(x_gauss)
    Kp_plus = Kp + dkappa_perp_dB2
    Kp_plusplus = dkappa_perp_dB2 - (B01**2 * Kp_plus / B0**2)
    ! operators
    Fop_plus = k2 * B02 / eps + k3 * B03
    dFop_plus = (k2 / eps) * (dB02 - deps * B02 / eps) + k3 * dB03
    Gop_min = k3 * B02 - k2 * B03 / eps

    ! B01 modified F+ operator
    Fop_B01 = deps * ic * B01 / eps + Fop_plus

    ! ==================== Quadratic * Quadratic ====================
    call elements%add( &
      gamma_1 * dT0 * (B01 / B0**2) * dkappa_perp_drho * Fop_B01, &
      "T", &
      "rho", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      gamma_1 * ( &
        B01 * Fop_plus * diffKp &
        - B01 * dT0 * (dkappa_para_dT - dkappa_perp_dT) * Fop_B01 / B0**2 &
        + Kp * ( &
          B01 * deps * (2.0d0 * deps * ic * B01 / eps + 3.0d0 * Fop_plus) / eps &
          + B01 * dFop_plus &
          - ic * Fop_plus**2 &
        ) &
      ), &
      "T", &
      "T", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      2.0d0 * gamma_1 * eps * dT0 * Gop_min * Kp_plus * B01 * Fop_B01 / B0**2, &
      "T", &
      "a1", &
      spline1=h_quad, &
      spline2=h_quad &
    )
    ! ==================== dQuadratic * Quadratic ====================
    call elements%add( &
      ic * gamma_1 * dT0 * dkappa_perp_drho * B01**2 / B0**2, &
      "T", &
      "rho", &
      spline1=dh_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      gamma_1 * ( &
        B01 * Kp * (2.0d0 * deps * ic * B01 / eps + 3.0d0 * Fop_plus) &
        - ic * dT0 * ( &
          B01**2 * dkappa_para_dT / B0**2 - dkappa_perp_dT * B01**2 / B0**2 &
        ) &
      ), &
      "T", &
      "T", &
      spline1=dh_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      -2.0d0 * ic * gamma_1 * eps * dT0 * Gop_min * Kp_plusplus, &
      "T", &
      "a1", &
      spline1=dh_quad, &
      spline2=h_quad &
    )
    ! ==================== Quadratic * dQuadratic ====================
    call elements%add( &
      -2.0d0 * ic * gamma_1 * deps * B01**2 * Kp / eps, &
      "T", &
      "T", &
      spline1=h_quad, &
      spline2=dh_quad &
    )
    ! ==================== dQuadratic * dQuadratic ====================
    call elements%add( &
      -ic * gamma_1 * 2.0d0 * B01**2 * Kp, &
      "T", &
      "T", &
      spline1=dh_quad, &
      spline2=dh_quad &
    )
    ! ==================== Quadratic * Cubic ====================
    call elements%add( &
      gamma_1 * k3 * ( &
        Kp * (B01 * ddT0 + ic * dT0 * Fop_plus) &
        + B01 * dT0 * (diffKp - 2.0d0 * ic * B01 * Kp_plus * Fop_B01 / B0**2) &
      ), &
      "T", &
      "a2", &
      spline1=h_quad, &
      spline2=h_cubic &
    )
    call elements%add( &
      -gamma_1 * k2 * ( &
        Kp * (B01 * ddT0 + ic * dT0 * Fop_plus) &
        + B01 * dT0 * (diffKp - 2.0d0 * ic * B01 * Kp_plus * Fop_B01 / B0**2) &
      ), &
      "T", &
      "a3", &
      spline1=h_quad, &
      spline2=h_cubic &
    )
    ! ==================== Quadratic * dCubic ====================
    call elements%add( &
      gamma_1 * dT0 * B01 * (2.0d0 * B03 * Kp_plus * Fop_B01 / B0**2 - k3 * Kp), &
      "T", &
      "a2", &
      spline1=h_quad, &
      spline2=dh_cubic &
    )
    call elements%add( &
      gamma_1 * dT0 * B01 * ( &
        -2.0d0 * eps * B02 * Kp_plus * Fop_B01 / B0**2 + k2 * Kp &
      ), &
      "T", &
      "a3", &
      spline1=h_quad, &
      spline2=dh_cubic &
    )
    ! ==================== dQuadratic * Cubic ====================
    call elements%add( &
      -2.0d0 * gamma_1 * dT0 * B01 * k3 * Kp_plusplus, &
      "T", &
      "a2", &
      spline1=dh_quad, &
      spline2=h_cubic &
    )
    call elements%add( &
      2.0d0 * gamma_1 * dT0 * B01 * k2 * Kp_plusplus, &
      "T", &
      "a3", &
      spline1=dh_quad, &
      spline2=h_cubic &
    )
    ! ==================== dQuadratic * dCubic ====================
    call elements%add( &
      -2.0d0 * ic * gamma_1 * dT0 * B03 * Kp_plusplus, &
      "T", &
      "a2", &
      spline1=dh_quad, &
      spline2=dh_cubic &
    )
    call elements%add( &
      2.0d0 * ic * gamma_1 * dT0 * eps * B02 * Kp_plusplus, &
      "T", &
      "a3", &
      spline1=dh_quad, &
      spline2=dh_cubic &
    )
  end subroutine add_conduction_matrix_terms_bfield

end submodule smod_conduction_matrix
