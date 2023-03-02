submodule (mod_matrix_manager) smod_conduction_matrix
  use mod_equilibrium, only: kappa_field
  implicit none

contains

  module procedure add_conduction_matrix_terms
    use mod_matrix_shortcuts, only: get_Kp_operator, get_diffF_operator

    real(dp)  :: eps, deps
    real(dp)  :: dT0, ddT0
    real(dp)  :: B0, B01, B02, B03
    real(dp)  :: diffKp, Kp, Kp_plus, Kp_plusplus
    real(dp)  :: WVop, Fop_plus, dFop_plus, Gop_min
    real(dp)  :: kappa_para, dkappa_para_dT
    real(dp)  :: kappa_perp, dkappa_perp_drho, dkappa_perp_dT, dkappa_perp_dB2
    real(dp) :: gamma_1
    complex(dp) :: Fop_B01
    type(matrix_elements_t) :: elements

    if (settings%physics%is_incompressible) return

    gamma_1 = settings%physics%get_gamma_1()
    ! grid variables
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! temperature variables
    dT0 = T_field % d_T0_dr(gauss_idx)
    ddT0 = T_field % dd_T0_dr(gauss_idx)
    ! magnetic field variables
    B0 = B_field % B0(gauss_idx)
    B01 = B_field % B01
    B02 = B_field % B02(gauss_idx)
    B03 = B_field % B03(gauss_idx)
    ! prefactors
    Kp = kappa_field % prefactor(gauss_idx)
    diffKp = kappa_field % d_prefactor_dr(gauss_idx)
    Kp_plus = get_Kp_operator(gauss_idx, which="+")
    Kp_plusplus = get_Kp_operator(gauss_idx, which="++")
    ! parallel thermal conduction variables
    kappa_para = kappa_field % kappa_para(gauss_idx)
    dkappa_para_dT = kappa_field % d_kappa_para_dT(gauss_idx)
    ! perpendicular thermal conduction variables
    kappa_perp = kappa_field % kappa_perp(gauss_idx)
    dkappa_perp_drho = kappa_field % d_kappa_perp_drho(gauss_idx)
    dkappa_perp_dT = kappa_field % d_kappa_perp_dT(gauss_idx)
    dkappa_perp_dB2 = kappa_field % d_kappa_perp_dB2(gauss_idx)
    ! operators
    WVop = get_wv_operator(gauss_idx)
    Fop_plus = get_F_operator(gauss_idx, which="plus")
    dFop_plus = get_diffF_operator(gauss_idx, which="plus")
    Gop_min = get_G_operator(gauss_idx, which="minus")

    ! B01 modified F+ operator
    Fop_B01 = deps * ic * B01 / eps + Fop_plus

    elements = new_matrix_elements(state_vector=settings%get_state_vector())

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
        - ic * WVop * kappa_perp / eps &
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
      -ic * gamma_1 * dT0 * dkappa_perp_drho * (1.0d0 - B01**2 / B0**2), &
      "T", &
      "rho", &
      spline1=dh_quad, &
      spline2=h_quad &
    )
    call elements%add( &
      gamma_1 * ( &
        B01 * Kp * (2.0d0 * deps * ic * B01 / eps + 3.0d0 * Fop_plus) &
        + ic * deps * kappa_perp / eps &
        - ic * dT0 * ( &
          B01**2 * dkappa_para_dT / B0**2 + dkappa_perp_dT * (1.0d0 - B01**2 / B0**2) &
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
      -ic * gamma_1 * (2.0d0 * B01**2 * Kp + kappa_perp), &
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

    call add_to_quadblock(quadblock, elements, current_weight, settings%dims)
    call elements%delete()
  end procedure add_conduction_matrix_terms

end submodule smod_conduction_matrix
