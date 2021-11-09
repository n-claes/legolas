submodule (mod_matrix_manager) smod_regular_matrix
  implicit none

  !> list containing matrix elements
  type(matrix_elements_t) :: elements

contains

  module procedure add_bmatrix_terms
    real(dp)  :: rho, eps

    rho = rho_field % rho0(gauss_idx)
    eps = eps_grid(gauss_idx)

    ! Quadratic * Quadratic
    call elements%add(1.0d0, location=["rho", "rho"])
    call elements%add(eps * rho, location=["v2", "v2"])
    call elements%add(rho, location=["v3", "v3"])
    call elements%add(rho, location=["T", "T"])
    call elements%add(eps, location=["a1", "a1"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      h_quad, &
      h_quad &
    )
    call elements%delete()

    ! Cubic * Cubic
    call elements%add(rho, location=["v1", "v1"])
    call elements%add(1.0d0, location=["a2", "a2"])
    call elements%add(eps, location=["a3", "a3"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      h_cubic, &
      h_cubic &
    )
    call elements%delete()
  end procedure add_bmatrix_terms


  module procedure add_regular_matrix_terms
    use mod_global_variables, only: external_gravity
    use mod_equilibrium, only: grav_field

    real(dp)  :: eps, deps
    real(dp)  :: rho, drho
    real(dp)  :: T0, dT0
    real(dp)  :: B01, B02, dB02, drB02, B03, db03
    real(dp)  :: Fop_plus, Gop_plus, Gop_min, WVop
    real(dp)  :: element

    ! grid variables
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! density variables
    rho = rho_field % rho0(gauss_idx)
    drho = rho_field % d_rho0_dr(gauss_idx)
    ! temperature variables
    T0 = T_field % T0(gauss_idx)
    dT0 = T_field % d_T0_dr(gauss_idx)
    ! magnetic field variables
    B01 = B_field % B01
    B02 = B_field % B02(gauss_idx)
    dB02 = B_field % d_B02_dr(gauss_idx)
    drB02 = deps * B02 + eps * dB02
    B03 = B_field % B03(gauss_idx)
    dB03 = B_field % d_B03_dr(gauss_idx)
    ! operators
    Fop_plus = get_F_operator(gauss_idx, which="plus")
    Gop_plus = get_G_operator(gauss_idx, which="plus")
    Gop_min = get_G_operator(gauss_idx, which="minus")
    WVop = get_wv_operator(gauss_idx)

    ! ==================== Quadratic * Cubic ====================
    call elements%add(-drho, location=["rho", "v1 "])
    call elements%add(k3 * (drB02 - ic * k2 * B01) / eps, location=["v2", "a2"])
    call elements%add(k2 * (ic * k2 * B01 - drB02) / eps, location=["v2", "a3"])
    call elements%add(k3 * (dB03 - ic * k3 * B01), location=["v3", "a2"])
    call elements%add(k2 * (ic * B01 * k3 - dB03), location=["v3", "a3"])
    call elements%add(-dT0 * rho, location=["T ", "v1"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      h_quad, &
      h_cubic &
    )
    call elements%delete()

    ! ==================== Quadratic * dCubic ====================
    call elements%add(-rho, location=["rho", "v1 "])
    call elements%add(k2 * B03 / eps, location=["v2", "a2"])
    call elements%add(eps * k3 * B03, location=["v2", "a3"])
    call elements%add(-(k2 * B02 + ic * deps * B01) / eps, location=["v3", "a2"])
    call elements%add(-eps * k3 * B02, location=["v3", "a3"])
    call elements%add(-gamma_1 * T0 * rho, location=["T ", "v1"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      h_quad, &
      dh_cubic &
    )
    call elements%delete()

    ! ==================== Quadratic * Quadratic ====================
    call elements%add(rho * k2, location=["rho", "v2 "])
    call elements%add(rho * k3, location=["rho", "v3 "])
    call elements%add(k2 * T0 / eps, location=["v2 ", "rho"])
    call elements%add(k2 * rho / eps, location=["v2", "T "])
    call elements%add(-WVop * B03, location=["v2", "a1"])
    call elements%add(k3 * T0, location=["v3 ", "rho"])
    call elements%add(k3 * rho, location=["v3", "T "])
    call elements%add(ic * deps * k2 * B01 / eps + B02 * WVop, location=["v3", "a1"])
    call elements%add(gamma_1 * k2 * rho * T0, location=["T ", "v2"])
    call elements%add(gamma_1 * k3 * rho * T0, location=["T ", "v3"])
    call elements%add(-eps * B03, location=["a1", "v2"])
    call elements%add(B02, location=["a1", "v3"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      h_quad, &
      h_quad &
    )
    call elements%delete()

    ! ==================== Cubic * Quadratic ====================
    element = -deps * T0 / eps
    if (external_gravity) then
      ! adds gravity term to A(2, 1) matrix element
      element = element + grav_field % grav(gauss_idx)
    end if
    call elements%add(element, location=["v1 ", "rho"])
    call elements%add(-deps * rho / eps, location=["v1", "T "])
    call elements%add(deps * Gop_plus, location=["v1", "a1"])
    call elements%add(ic * B01, location=["a2", "v3"])
    call elements%add(-ic * eps * B01, location=["a3", "v2"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      h_cubic, &
      h_quad &
    )
    call elements%delete()

    ! ==================== dCubic * Quadratic ====================
    call elements%add(-T0, location=["v1 ", "rho"])
    call elements%add(-rho, location=["v1", "T "])
    call elements%add(-eps * Gop_min, location=["v1", "a1"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      dh_cubic, &
      h_quad &
    )
    call elements%delete()

    ! ==================== Cubic * Cubic ====================
    call elements%add(-k3 * Fop_plus, location=["v1", "a2"])
    call elements%add(k2 * Fop_plus, location=["v1", "a3"])
    call elements%add(-B03, location=["a2", "v1"])
    call elements%add(B02, location=["a3", "v1"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      h_cubic, &
      h_cubic &
    )
    call elements%delete()

    ! ==================== Cubic * dCubic ====================
    call elements%add(-deps * B03 / eps, location=["v1", "a2"])
    ! A(2, 8)
    call elements%add(-deps * B02, location=["v1", "a3"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      h_cubic, &
      dh_cubic &
    )
    call elements%delete()

    ! ==================== dCubic * dCubic ====================
    call elements%add(-B03, location=["v1", "a2"])
    call elements%add(eps * B02, location=["v1", "a3"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      dh_cubic, &
      dh_cubic &
    )
    call elements%delete()

    ! ==================== dQuadratic * Quadratic ====================
    call elements%add(-ic * eps * k3 * B01, location=["v2", "a1"])
    call elements%add(ic * k2 * B01, location=["v3", "a1"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      dh_quad, &
      h_quad &
    )
    call elements%delete()

    ! ==================== dQuadratic * dCubic ====================
    call elements%add(ic * eps * B01, location=["v2", "a3"])
    call elements%add(-ic * B01, location=["v3", "a2"])
    call subblock( &
      quadblock, &
      elements%get_values(), &
      elements%get_positions(), &
      current_weight, &
      dh_quad, &
      dh_cubic &
    )
    call elements%delete()
  end procedure add_regular_matrix_terms

end submodule smod_regular_matrix
