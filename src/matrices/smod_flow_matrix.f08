submodule (mod_matrix_manager) smod_flow_matrix
  use mod_equilibrium, only: v_field
  implicit none

contains

  module procedure add_flow_matrix_terms
    use mod_global_variables, only: incompressible

    real(dp)  :: eps, deps
    real(dp)  :: rho, drho
    real(dp)  :: T0
    real(dp)  :: v01, dv01, drv01
    real(dp)  :: v02, dv02, drv02
    real(dp)  :: v03, dv03
    real(dp)  :: Vop

    ! grid variables
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! density variables
    rho = rho_field % rho0(gauss_idx)
    drho = rho_field % d_rho0_dr(gauss_idx)
    ! temperature variables
    T0 = T_field % T0(gauss_idx)
    ! flow variables
    v01 = v_field % v01(gauss_idx)
    dv01 = v_field % d_v01_dr(gauss_idx)
    drv01 = deps * v01 + eps * dv01
    v02 = v_field % v02(gauss_idx)
    dv02 = v_field % d_v02_dr(gauss_idx)
    drv02 = deps * v02 + eps * dv02
    v03 = v_field % v03(gauss_idx)
    dv03 = v_field % d_v03_dr(gauss_idx)
    Vop = get_flow_operator(gauss_idx)

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=8)
    ! Phi(1, 1)
    factors(1) = Vop - ic * dv01
    positions(1, :) = [1, 1]
    ! Phi(3, 1)
    factors(2) = -drv02 * ic * v01 / eps
    positions(2, :) = [3, 1]
    ! Phi(3, 3)
    factors(3) = rho * (eps * Vop - ic * deps * v01)
    positions(3, :) = [3, 3]
    ! Phi(4, 1)
    factors(4) = -ic * v01 * dv03
    positions(4, :) = [4, 1]
    ! Phi(4, 4)
    factors(5) = rho * (Vop + ic * dv01) + (deps * rho / eps + drho) * ic * v01
    positions(5, :) = [4, 4]
    ! Phi(5, 1)
    factors(6) = 0.0d0
    positions(6, :) = [5, 1]
    ! Phi(5, 5)
    factors(7) = 0.0d0
    positions(7, :) = [5, 5]
    if (.not. incompressible) then
      factors(6) = -ic * gamma_1 * drv01 * T0 / eps
      factors(7) = ( &
        rho * (Vop + ic * dv01 - ic * gamma_1 * drv01 / eps) &
        + ic * v01 * (deps * rho / eps + drho) &
      )
    end if
    ! Phi(6, 6)
    factors(8) = eps * Vop
    positions(8, :) = [6, 6]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_quad)

    ! ==================== Quadratic * dQuadratic ====================
    call reset_factor_positions(new_size=2)
    ! Phi(1, 1)
    factors(1) = -ic * v01
    positions(1, :) = [1, 1]
    ! Phi(3, 3)
    factors(2) = -ic * eps * rho * v01
    positions(2, :) = [3, 3]
    call subblock(quadblock, factors, positions, current_weight, h_quad, dh_quad)

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=4)
    ! Phi(2, 1)
    factors(1) = v01 * dv01 - deps * v02**2 / eps
    positions(1, :) = [2, 1]
    ! Phi(2, 3)
    factors(2) = -2.0d0 * deps * rho * v02
    positions(2, :) = [2, 3]
    ! Phi(7, 6)
    factors(3) = ic * v01 * k2
    positions(3, :) = [7, 6]
    ! Phi(8, 6)
    factors(4) = eps * k3 * ic * v01
    positions(4, :) = [8, 6]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_quad)

    ! ==================== Cubic * Cubic ====================
    call reset_factor_positions(new_size=5)
    ! Phi(2, 2)
    factors(1) = rho * Vop + (deps * rho / eps + drho) * ic * v01
    positions(1, :) = [2, 2]
    ! Phi(7, 7)
    factors(2) = k3 * v03
    positions(2, :) = [7, 7]
    ! Phi(7, 8)
    factors(3) = -k2 * v03
    positions(3, :) = [7, 8]
    ! Phi(8, 7)
    factors(4) = - k3 * v02
    positions(4, :) = [8, 7]
    ! Phi(8, 8)
    factors(5) = k2 * v02
    positions(5, :) = [8, 8]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_cubic)

    ! ==================== dCubic * Cubic ====================
    call reset_factor_positions(new_size=1)
    ! Phi(2, 2)
    factors(1) = ic * rho * v01
    positions(1, :) = [2, 2]
    call subblock(quadblock, factors, positions, current_weight, dh_cubic, h_cubic)

    ! ==================== Quadratic * Cubic ====================
    call reset_factor_positions(new_size=2)
    ! Phi(3, 2)
    factors(1) = -drv02 * rho / eps
    positions(1, :) = [3, 2]
    ! Phi(4, 2)
    factors(2) = -rho * dv03
    positions(2, :) = [4, 2]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_cubic)

    ! ==================== dQuadratic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! Phi(4, 4)
    factors(1) = ic * rho * v01
    positions(1, :) = [4, 4]
    ! Phi(5, 5)
    factors(2) = 0.0d0
    if (.not. incompressible) then
      factors(2) = ic * rho * v01
    end if
    positions(2, :) = [5, 5]
    call subblock(quadblock, factors, positions, current_weight, dh_quad, h_quad)

    ! ==================== Quadratic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! Phi(6, 7)
    factors(1) = -v02
    positions(1, :) = [6, 7]
    ! Phi(6, 8)
    factors(2) = -eps * v03
    positions(2, :) = [6, 8]
    call subblock(quadblock, factors, positions, current_weight, h_quad, dh_cubic)

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! Phi(7, 7)
    factors(1) = -ic * v01
    positions(1, :) = [7, 7]
    ! Phi(8, 8)
    factors(2) = - ic * v01
    positions(2, :) = [8, 8]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, dh_cubic)

  end procedure add_flow_matrix_terms


  !> Calculates the $$\boldsymbol{\mathcal{V}}$$ operator, given as
  !! $$ \boldsymbol{\mathcal{V}} = \left(\frac{k_2}{\eps}v_{02} + k_3v_{03}\right) $$
  function get_flow_operator(gauss_idx) result(Voperator)
    !> current index in the Gaussian grid
    integer, intent(in) :: gauss_idx
    !> the V operator on return
    real(dp)  :: Voperator

    Voperator = ( &
      k2 * v_field % v02(gauss_idx) / eps_grid(gauss_idx) &
      + k3 * v_field % v03(gauss_idx) &
    )
  end function get_flow_operator

end submodule smod_flow_matrix
