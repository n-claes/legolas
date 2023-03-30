submodule (mod_matrix_manager) smod_resistive_matrix
  implicit none

contains

  module procedure add_resistive_matrix_terms
    real(dp)  :: eps, deps
    real(dp)  :: B02, dB02, drB02, ddB02
    real(dp)  :: B03, dB03, ddB03
    real(dp)  :: eta, detadT, deta
    real(dp)  :: WVop, Rop_pos, Rop_neg
    real(dp) :: gamma_1
    type(matrix_elements_t) :: elements

    gamma_1 = settings%physics%get_gamma_1()

    ! grid variables
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! magnetic field variables
    B02 = background%magnetic%B02(x_gauss)
    dB02 = background%magnetic%dB02(x_gauss)
    drB02 = deps * B02 + eps * dB02
    ddB02 = background%magnetic%ddB02(x_gauss)
    B03 = background%magnetic%B03(x_gauss)
    dB03 = background%magnetic%dB03(x_gauss)
    ddB03 = background%magnetic%ddB03(x_gauss)
    ! resistivity variables
    eta = physics%resistivity%eta(x_gauss)
    detadT = physics%resistivity%detadT(x_gauss)
    ! total derivative eta = deta_dr + dT0_dr * deta_dT
    deta = ( &
      physics%resistivity%detadr(x_gauss) &
      + (background%temperature%dT0(x_gauss) * detadT) &
    )

    WVop = k2**2 / eps + eps * k3**2
    Rop_pos = deps * eta / eps + deta
    Rop_neg = deps * eta / eps - deta

    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Quadratic * Quadratic ====================
    call elements%add(-ic * eta * WVop, "a1", "a1", spline1=h_quad, spline2=h_quad)

    ! ==================== Quadratic * dCubic ====================
    call elements%add(ic * eta * k2 / eps, "a1", "a2", spline1=h_quad, spline2=dh_cubic)
    call elements%add(ic * eta * eps * k3, "a1", "a3", spline1=h_quad, spline2=dh_cubic)

    ! ==================== Cubic * Quadratic ====================
    call elements%add(ic * dB03 * detadT, "a2", "T", spline1=h_cubic, spline2=h_quad)
    call elements%add(ic * k2 * Rop_pos, "a2", "a1", spline1=h_cubic, spline2=h_quad)
    call elements%add( &
      -ic * drB02 * detadT / eps, "a3", "T", spline1=h_cubic, spline2=h_quad &
    )
    call elements%add(ic * deta * eps * k3, "a3", "a1", spline1=h_cubic, spline2=h_quad)

    ! ==================== dCubic * Quadratic ====================
    call elements%add(ic * eta * k2, "a2", "a1", spline1=dh_cubic, spline2=h_quad)
    call elements%add(ic * eta * eps * k3, "a3", "a1", spline1=dh_cubic, spline2=h_quad)

    ! ==================== Cubic * Cubic ====================
    call elements%add(-ic * eta * k3**2, "a2", "a2", spline1=h_cubic, spline2=h_cubic)
    call elements%add(ic * eta * k2 * k3, "a2", "a3", spline1=h_cubic, spline2=h_cubic)
    call elements%add( &
      ic * eta * k2 * k3 / eps, "a3", "a2", spline1=h_cubic, spline2=h_cubic &
    )
    call elements%add( &
      -ic * eta * k2**2 / eps, "a3", "a3", spline1=h_cubic, spline2=h_cubic &
    )

    ! ==================== Cubic * dCubic ====================
    call elements%add(-ic * Rop_pos, "a2", "a2", spline1=h_cubic, spline2=dh_cubic)
    call elements%add(-ic * deta * eps, "a3", "a3", spline1=h_cubic, spline2=dh_cubic)

    ! ==================== dCubic * dCubic ====================
    call elements%add(-ic * eta, "a2", "a2", spline1=dh_cubic, spline2=dh_cubic)
    call elements%add(-ic * eta * eps, "a3", "a3", spline1=dh_cubic, spline2=dh_cubic)

    if (.not. settings%physics%is_incompressible) then
      call add_compressible_resistive_terms()
    end if

    call add_to_quadblock(quadblock, elements, weight, settings%dims)
    call elements%delete()

  contains

    subroutine add_compressible_resistive_terms()
      ! ==================== Quadratic * Quadratic ====================
      call elements%add( &
        ic * gamma_1 * detadT * ((drB02 / eps)**2 + dB03**2), &
        "T", &
        "T", &
        spline1=h_quad, &
        spline2=h_quad &
      )
      call elements%add( &
        2.0d0 * ic * gamma_1 * ( &
          k2 * (dB03 * Rop_pos + eta * ddB03) &
          + k3 * (drB02 * Rop_neg - eta * (2.0d0 * deps * dB02 + eps * ddB02)) &
        ), &
        "T", &
        "a1", &
        spline1=h_quad, &
        spline2=h_quad &
      )
      ! ==================== dQuadratic * Quadratic ====================
      call elements%add( &
        -2.0d0 * ic * gamma_1 * eta * (k3 * drB02 - k2 * dB03), &
        "T", &
        "a1", &
        spline1=dh_quad, &
        spline2=h_quad &
      )
      ! ==================== Quadratic * Cubic ====================
      call elements%add( &
        -2.0d0 * ic * gamma_1 * eta * (drB02 * k2 * k3 / eps**2 + k3**2 * dB03), &
        "T", &
        "a2", &
        spline1=h_quad, &
        spline2=h_cubic &
      )
      call elements%add( &
        2.0d0 * ic * gamma_1 * eta * (drB02 * k2**2 / eps**2 + k2 * k3 * dB03), &
        "T", &
        "a3", &
        spline1=h_quad, &
        spline2=h_cubic &
      )
      ! ==================== Quadratic * dCubic ====================
      call elements%add( &
        -2.0d0 * ic * gamma_1 * (dB03 * Rop_pos + ddB03 * eta), &
        "T", &
        "a2", &
        spline1=h_quad, &
        spline2=dh_cubic &
      )
      call elements%add( &
        -2.0d0 * ic * gamma_1 * ( &
          drB02 * Rop_neg - eta * (2.0d0 * deps * dB02 + eps * ddB02) &
        ), &
        "T", &
        "a3", &
        spline1=h_quad, &
        spline2=dh_cubic &
      )
      ! ==================== dQuadratic * dCubic ====================
      call elements%add( &
        -2.0d0 * ic * gamma_1 * eta * dB03, &
        "T", &
        "a2", &
        spline1=dh_quad, &
        spline2=dh_cubic &
      )
      call elements%add( &
        2.0d0 * ic * gamma_1 * drB02 * eta, &
        "T", &
        "a3", &
        spline1=dh_quad, &
        spline2=dh_cubic &
      )
    end subroutine add_compressible_resistive_terms

  end procedure add_resistive_matrix_terms

end submodule smod_resistive_matrix
