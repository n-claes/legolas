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

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=3)
    ! R(5, 5)
    factors(1) = 0.0d0
    positions(1, :) = [5, 5]
    ! R(5, 6)
    factors(2) = 0.0d0
    positions(2, :) = [5, 6]
    if (.not. settings%physics%is_incompressible) then
      factors(1) = ic * gamma_1 * detadT * ((drB02 / eps)**2 + dB03**2)
      factors(2) = 2.0d0 * ic * gamma_1 * ( &
        k2 * (dB03 * Rop_pos + eta * ddB03) &
        + k3 * (drB02 * Rop_neg - eta * (2.0d0 * deps * dB02 + eps * ddB02)) &
      )
    end if
    ! R(6, 6)
    factors(3) = -ic * eta * WVop
    positions(3, :) = [6, 6]
    call subblock( &
      quadblock, factors, positions, current_weight, h_quad, h_quad, settings%dims &
    )

    if (.not. settings%physics%is_incompressible) then
      ! ==================== dQuadratic * Quadratic ====================
      call reset_factor_positions(new_size=1)
      ! R(5, 6)
      factors(1) = -2.0d0 * ic * gamma_1 * eta * (k3 * drB02 - k2 * dB03)
      positions(1, :) = [5, 6]
      call subblock( &
        quadblock, factors, positions, current_weight, dh_quad, h_quad, settings%dims &
      )

      ! ==================== Quadratic * Cubic ====================
      call reset_factor_positions(new_size=2)
      ! R(5, 7)
      factors(1) = -2.0d0 * ic * gamma_1 * eta * (drB02 * k2 * k3 / eps**2 + k3**2 * dB03)
      positions(1, :) = [5, 7]
      ! R(5, 8)
      factors(2) = 2.0d0 * ic * gamma_1 * eta * (drB02 * k2**2 / eps**2 + k2 * k3 * dB03)
      positions(2, :) = [5, 8]
      call subblock( &
        quadblock, factors, positions, current_weight, h_quad, h_cubic, settings%dims &
      )
    end if

    ! ==================== Quadratic * dCubic ====================
    call reset_factor_positions(new_size=4)
    ! R(5, 7)
    factors(1) = 0.0d0
    positions(1, :) = [5, 7]
    ! R(5, 8)
    factors(2) = 0.0d0
    positions(2, :) = [5, 8]
    if (.not. settings%physics%is_incompressible) then
      factors(1) = -2.0d0 * ic * gamma_1 * (dB03 * Rop_pos + ddB03 * eta)
      factors(2) = -2.0d0 * ic * gamma_1 * ( &
        drB02 * Rop_neg - eta * (2.0d0 * deps * dB02 + eps * ddB02) &
      )
    end if
    ! R(6, 7)
    factors(3) = ic * eta * k2 / eps
    positions(3, :) = [6, 7]
    ! R(6, 8)
    factors(4) = ic * eta * eps * k3
    positions(4, :) = [6, 8]
    call subblock( &
      quadblock, factors, positions, current_weight, h_quad, dh_cubic, settings%dims &
    )

    if (.not. settings%physics%is_incompressible) then
      ! ==================== dQuadratic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! R(5, 7)
      factors(1) = -2.0d0 * ic * gamma_1 * eta * dB03
      positions(1, :) = [5, 7]
      ! R(5, 8)
      factors(2) = 2.0d0 * ic * gamma_1 * drB02 * eta
      positions(2, :) = [5, 8]
      call subblock( &
        quadblock, &
        factors, &
        positions, &
        current_weight, &
        dh_quad, &
        dh_cubic, &
        settings%dims &
      )
    end if

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=4)
    ! R(7, 5)
    factors(1) = ic * dB03 * detadT
    positions(1, :) = [7, 5]
    ! R(7, 6)
    factors(2) = ic * k2 * Rop_pos
    positions(2, :) = [7, 6]
    ! R(8, 5)
    factors(3) = -ic * drB02 * detadT / eps
    positions(3, :) = [8, 5]
    ! R(8, 6)
    factors(4) = ic * deta * eps * k3
    positions(4, :) = [8, 6]
    call subblock( &
      quadblock, factors, positions, current_weight, h_cubic, h_quad, settings%dims &
    )

    ! ==================== dCubic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! R(7, 6)
    factors(1) = ic * eta * k2
    positions(1, :) = [7, 6]
    ! R(8, 6)
    factors(2) = ic * eta * eps * k3
    positions(2, :) = [8, 6]
    call subblock( &
      quadblock, factors, positions, current_weight, dh_cubic, h_quad, settings%dims &
    )

    ! ==================== Cubic * Cubic ====================
    call reset_factor_positions(new_size=4)
    ! R(7, 7)
    factors(1) = -ic * eta * k3**2
    positions(1, :) = [7, 7]
    ! R(7, 8)
    factors(2) = ic * eta * k2 * k3
    positions(2, :) = [7, 8]
    ! R(8, 7)
    factors(3) = ic * eta * k2 * k3 / eps
    positions(3, :) = [8, 7]
    ! R(8, 8)
    factors(4) = -ic * eta * k2**2 / eps
    positions(4, :) = [8, 8]
    call subblock( &
      quadblock, factors, positions, current_weight, h_cubic, h_cubic, settings%dims &
    )

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! R(7, 7)
    factors(1) = -ic * Rop_pos
    positions(1, :) = [7, 7]
    ! R(8, 8)
    factors(2) = -ic * deta * eps
    positions(2, :) = [8, 8]
    call subblock( &
      quadblock, factors, positions, current_weight, h_cubic, dh_cubic, settings%dims &
    )

    ! ==================== dCubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! R(7, 7)
    factors(1) = -ic * eta
    positions(1, :) = [7, 7]
    ! R(8, 8)
    factors(2) = -ic * eta * eps
    positions(2, :) = [8, 8]
    call subblock( &
      quadblock, factors, positions, current_weight, dh_cubic, dh_cubic, settings%dims &
    )

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
