submodule (mod_matrix_manager) smod_viscosity_matrix
  use mod_global_variables, only: viscosity_value, viscous_heating, incompressible
  use mod_equilibrium, only: v_field
  implicit none

contains

  module procedure add_viscosity_matrix_terms
    real(dp)  :: eps, deps
    real(dp)  :: v01, dv01, ddv01
    real(dp)  :: v02, dv02
    real(dp)  :: v03, dv03, ddv03
    real(dp)  :: mu
    real(dp)  :: WVop

    ! grid variables
    eps = eps_grid(gauss_idx)
    deps = d_eps_grid_dr(gauss_idx)
    ! viscous heating variables
    v01 = v_field % v01(gauss_idx)
    dv01 = v_field % d_v01_dr(gauss_idx)
    ddv01 = v_field % dd_v01_dr(gauss_idx)
    v02 = v_field % v02(gauss_idx)
    dv02 = v_field % d_v02_dr(gauss_idx)
    v03 = v_field % v03(gauss_idx)
    dv03 = v_field % d_v03_dr(gauss_idx)
    ddv03 = v_field % dd_v03_dr(gauss_idx)
    ! viscosity value
    mu = viscosity_value
    ! operators
    WVop = get_wv_operator(gauss_idx)

    ! ==================== Cubic * Cubic ====================
    call reset_factor_positions(new_size=1)
    ! Sigma(2, 2)
    factors(1) = -ic * mu * (deps / eps + WVop) / eps
    positions(1, :) = [2, 2]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_cubic)

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=1)
    ! Sigma(2, 2)
    factors(1) = -ic * mu * deps / (3.0d0 * eps)
    positions(1, :) = [2, 2]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, dh_cubic)

    ! ==================== dCubic * Cubic ====================
    call reset_factor_positions(new_size=1)
    ! Sigma(2, 2)
    factors(1) = ic * mu * deps / eps
    positions(1, :) = [2, 2]
    call subblock(quadblock, factors, positions, current_weight, dh_cubic, h_cubic)

    ! ==================== dCubic * dCubic ====================
    call reset_factor_positions(new_size=1)
    ! Sigma(2, 2)
    factors(1) = -4.0d0 * ic * mu / 3.0d0
    positions(1, :) = [2, 2]
    call subblock(quadblock, factors, positions, current_weight, dh_cubic, dh_cubic)

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! Sigma(2, 3)
    factors(1) = 7.0d0 * deps * ic * mu * k2 / (3.0d0 * eps)
    positions(1, :) = [2, 3]
    ! Sigma(2, 4)
    factors(2) = ic * mu * deps * k3 / (3.0d0 * eps)
    positions(2, :) = [2, 4]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_quad)

    ! ==================== dCubic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! Sigma(2, 3)
    factors(1) = ic * mu * k2 / 3.0d0
    positions(1, :) = [2, 3]
    ! Sigma(2, 4)
    factors(2) = ic * mu * k3 / 3.0d0
    positions(2, :) = [2, 4]
    call subblock(quadblock, factors, positions, current_weight, dh_cubic, h_quad)

    ! ==================== Quadratic * Cubic ====================
    call reset_factor_positions(new_size=2)
    ! Sigma(3, 2)
    factors(1) = ic * mu * deps * 2.0d0 * k2 / eps**2
    positions(1, :) = [3, 2]
    ! Sigma(5, 2)
    factors(2) = (0.0d0, 0.0d0)
    if (viscous_heating .and. (.not. incompressible)) then
      factors(2) = factors(2) + 2.0d0 * gamma_1 * mu * ( &
        (deps**2 * v01 - ic * deps * k2 * v02) / eps**2 - deps * dv01 / eps - ddv01 &
      )
    end if
    positions(2, :) = [5, 2]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_cubic)

    ! ==================== Quadratic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! Sigma(3, 2)
    factors(1) = ic * mu * k2 / (3.0d0 * eps)
    positions(1, :) = [3, 2]
    ! Sigma(4, 2)
    factors(2) = ic * mu * k3 / 3.0d0
    positions(2, :) = [4, 2]
    call subblock(quadblock, factors, positions, current_weight, h_quad, dh_cubic)

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=6)
    ! Sigma(3, 3)
    factors(1) = -ic * mu * (deps / eps + 4.0d0 * k2**2 / (3.0d0 * eps) + eps * k3**2)
    positions(1, :) = [3, 3]
    ! Sigma(3, 4)
    factors(2) = -ic * mu * k2 * k3 / (3.0d0 * eps)
    positions(2, :) = [3, 4]
    ! Sigma(4, 3)
    factors(3) = -ic * mu * k2 * k3 / 3.0d0
    positions(3, :) = [4, 3]
    ! Sigma(4, 4)
    factors(4) = -ic * mu * (k2**2 / eps**2 + 4.0d0 * k3**2 / 3.0d0)
    positions(4, :) = [4, 4]
    ! Sigma(5, 3)
    factors(5) = (0.0d0, 0.0d0)
    if (viscous_heating .and. (.not. incompressible)) then
      factors(5) = factors(5) + 2.0d0 * gamma_1 * mu * ( &
        deps**2 * ic * v02 - deps * k2 * v01 &
      ) / eps
    end if
    positions(5, :) = [5, 3]
    ! Sigma(5, 4)
    factors(6) = (0.0d0, 0.0d0)
    if (viscous_heating .and. (.not. incompressible)) then
      factors(6) = factors(6) - 2.0d0 * ic * gamma_1 * mu * (deps * dv03 / eps + ddv03)
    end if
    positions(6, :) = [5, 4]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_quad)

    ! ==================== dQuadratic * dQuadratic ====================
    call reset_factor_positions(new_size=2)
    ! Sigma(3, 3)
    factors(1) = -ic * mu * eps
    positions(1, :) = [3, 3]
    ! Sigma(4, 4)
    factors(2) = -ic * mu
    positions(2, :) = [4, 4]
    call subblock(quadblock, factors, positions, current_weight, dh_quad, dh_quad)

    ! ==================== dQuadratic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! Sigma(4, 4)
    factors(1) = ic * mu * deps / eps
    positions(1, :) = [4, 4]
    ! Sigma(5, 4)
    factors(2) = (0.0d0, 0.0d0)
    if (viscous_heating .and. (.not. incompressible)) then
      factors(2) = factors(2) - 2.0d0 * ic * gamma_1 * mu * dv03
    end if
    positions(2, :) = [5, 4]
    call subblock(quadblock, factors, positions, current_weight, dh_quad, h_quad)

    ! ==================== dQuadratic * Cubic ====================
    call reset_factor_positions(new_size=1)
    ! Sigma(5, 2)
    factors(1) = (0.0d0, 0.0d0)
    if (viscous_heating .and. (.not. incompressible)) then
      factors(1) = factors(1) - 2.0d0 * gamma_1 * mu * dv01
    end if
    positions(1, :) = [5, 2]
    call subblock(quadblock, factors, positions, current_weight, dh_quad, h_cubic)

    ! ==================== Quadratic * dQuadratic ====================
    call reset_factor_positions(new_size=1)
    ! Sigma(5, 3)
    factors(1) = (0.0d0, 0.0d0)
    if (viscous_heating .and. (.not. incompressible)) then
      factors(1) = factors(1) + 2.0d0 * ic * gamma_1 * mu * eps * dv02
    end if
    positions(1, :) = [5, 3]
    call subblock(quadblock, factors, positions, current_weight, h_quad, dh_quad)

  end procedure add_viscosity_matrix_terms

end submodule smod_viscosity_matrix
