submodule (mod_matrix_manager) smod_resistive_matrix
  use mod_equilibrium, only: eta_field
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
    B02 = B_field % B02(gauss_idx)
    dB02 = B_field % d_B02_dr(gauss_idx)
    drB02 = deps * B02 + eps * dB02
    ddB02 = eta_field % dd_B02_dr(gauss_idx)
    B03 = B_field % B03(gauss_idx)
    dB03 = B_field % d_B03_dr(gauss_idx)
    ddB03 = eta_field % dd_B03_dr(gauss_idx)
    ! resistivity variables
    eta = eta_field % eta(gauss_idx)
    detadT = eta_field % d_eta_dT(gauss_idx)
    deta = get_deta(gauss_idx)

    WVop = get_wv_operator(gauss_idx)
    Rop_pos = get_R_operator(gauss_idx, which="plus")
    Rop_neg = get_R_operator(gauss_idx, which="minus")

    elements = new_matrix_elements(state_vector=settings%get_state_vector())

    ! ==================== Quadratic * Quadratic ====================
    if (.not. settings%physics%is_incompressible) then
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
    end if
    call elements%add(-ic * eta * WVop, "a1", "a1", spline1=h_quad, spline2=h_quad)

    if (.not. settings%physics%is_incompressible) then
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
    end if

    ! ==================== Quadratic * dCubic ====================
    if (.not. settings%physics%is_incompressible) then
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
    end if
    call elements%add(ic * eta * k2 / eps, "a1", "a2", spline1=h_quad, spline2=dh_cubic)
    call elements%add(ic * eta * eps * k3, "a1", "a3", spline1=h_quad, spline2=dh_cubic)

    if (.not. settings%physics%is_incompressible) then
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
    end if

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

    call add_to_quadblock(quadblock, elements, current_weight, settings%dims)
    call elements%delete()
  end procedure add_resistive_matrix_terms


  !> Calculates the $$\boldsymbol{\mathcal{R}}$$ operator, given as
  !! $$
  !! \boldsymbol{\mathcal{R}} =
  !!      \left(\frac{\varepsilon'}{\varepsilon}\eta_0' \pm \eta_0\right)
  !! $$
  function get_R_operator(gauss_idx, which) result(Roperator)
    use mod_global_variables, only: NaN

    !> current index in the Gaussian grid
    integer, intent(in) :: gauss_idx
    !> which operator to calculate, <tt>"plus"</tt> or <tt>"minus"</tt>
    character(len=*), intent(in)  :: which
    !> the R operator on return
    real(dp)  :: Roperator

    real(dp)  :: eps, deps
    real(dp)  :: eta, deta

    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    eta = eta_field % eta(gauss_idx)
    deta = get_deta(gauss_idx)

    if (which == "plus") then
      Roperator = (deps * eta / eps + deta)
    else if (which == "minus") then
      Roperator = (deps * eta / eps - deta)
    else
      Roperator = NaN
      call logger%error("requesting invalid R-operator sign: " // trim(which))
    end if
  end function get_R_operator


  !> Calculates the total derivative of $$\eta$$, given as
  !! $$ \eta_0(r, T)' = \frac{d\eta_0}{dr} + \frac{dT_0}{dr}\frac{d\eta_0}{dT} $$
  function get_deta(gauss_idx)  result(deta)
    !> current intex in the Gaussian grid
    integer, intent(in) :: gauss_idx
    !> full eta derivative on return
    real(dp)  :: deta

    deta = eta_field % d_eta_dr(gauss_idx) + ( &
      T_field % d_T0_dr(gauss_idx) * eta_field % d_eta_dT(gauss_idx) &
    )
  end function get_deta

end submodule smod_resistive_matrix
