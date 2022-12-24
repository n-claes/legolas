submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_hall
  implicit none

contains

  module procedure add_natural_hall_Bterms
    use mod_equilibrium, only: hall_field

    real(dp)  :: eps
    real(dp)  :: rho
    real(dp)  :: eta_e

    if ( &
      .not. (settings%physics%hall%is_enabled() &
      .and. settings%physics%hall%has_electron_inertia()) &
    ) return

    eps = eps_grid(grid_idx)
    rho = rho_field % rho0(grid_idx)
    eta_e = hall_field % inertiafactor(grid_idx)

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! B_H(7, 6)
    factors(1) = eta_e * k2 / rho
    positions(1, :) = [7, 6]
    ! B_H(8, 6)
    factors(2) = eta_e * k3 * eps / rho
    positions(2, :) = [8, 6]
    call subblock(quadblock, factors, positions, weight, h_cubic, h_quad, settings%dims)

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! B_H(7, 7)
    factors(1) = -eta_e / rho
    positions(1, :) = [7, 7]
    ! B_H(8, 8)
    factors(2) = -eta_e * eps / rho
    positions(2, :) = [8, 8]
    call subblock( &
      quadblock, factors, positions, weight, h_cubic, dh_cubic, settings%dims &
    )

  end procedure add_natural_hall_Bterms

  module procedure add_natural_hall_terms
    use mod_equilibrium, only: hall_field

    real(dp)  :: eps, deps
    real(dp)  :: rho, T0, B01, B02, B03
    real(dp)  :: eta_H, mu, efrac

    if (.not. settings%physics%hall%is_enabled()) then
      return
    end if

    eps = eps_grid(grid_idx)
    deps = d_eps_grid_dr(grid_idx)

    rho = rho_field % rho0(grid_idx)
    T0 = T_field % T0(grid_idx)
    B01 = B_field % B01
    B02 = B_field % B02(grid_idx)
    B03 = B_field % B03(grid_idx)

    eta_H = hall_field % hallfactor(grid_idx)
    mu = 0.0_dp
    efrac = settings%physics%hall%get_electron_fraction()

    ! Hall by substitution of the momentum equation
    if (settings%physics%hall%is_using_substitution()) then
      if (settings%physics%viscosity%is_enabled()) then
        mu = settings%physics%viscosity%get_viscosity_value()
        ! ==================== Quadratic * Cubic ====================
        call reset_factor_positions(new_size=1)
        ! H(6, 2)
        factors(1) = -eta_H * ic * mu * deps / (eps * rho)
        positions(1, :) = [6, 2]
        call subblock( &
          quadblock, factors, positions, weight, h_quad, h_cubic, settings%dims &
        )

        ! ==================== Quadratic * dCubic ====================
        call reset_factor_positions(new_size=1)
        ! H(6, 2)
        factors(1) = 4.0d0 * eta_H * ic * mu / (3.0d0 * rho)
        positions(1, :) = [6, 2]
        call subblock( &
          quadblock, factors, positions, weight, h_quad, dh_cubic, settings%dims &
        )

        ! ==================== Cubic * dQuadratic ====================
        call reset_factor_positions(new_size=2)
        ! H(7, 3)
        factors(1) = eta_H * ic * mu * eps / rho
        positions(1, :) = [7, 3]
        ! H(8, 4)
        factors(2) = eta_H * ic * mu / rho
        positions(2, :) = [8, 4]
        call subblock( &
          quadblock, factors, positions, weight, h_cubic, dh_quad, settings%dims &
        )

        ! ==================== Cubic * Quadratic ====================
        call reset_factor_positions(new_size=1)
        ! H(8, 4)
        factors(1) = -eta_H * ic * mu * deps / (eps * rho)
        positions(1, :) = [8, 4]
        call subblock( &
          quadblock, factors, positions, weight, h_cubic, h_quad, settings%dims &
        )
      end if

    ! Hall without substitution, only E redefinition up to a gradient
    else
      ! ==================== Quadratic * Quadratic ====================
      call reset_factor_positions(new_size=1)
      ! H(6, 6)
      factors(1) = eta_H * (k2 * B03 - eps * k3 * B02) / rho
      positions(1, :) = [6, 6]
      call subblock( &
        quadblock, factors, positions, weight, h_quad, h_quad, settings%dims &
      )

      ! ==================== Quadratic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! H(6, 7)
      factors(1) = -eta_H * B03 / rho
      positions(1, :) = [6, 7]
      ! H(6, 8)
      factors(2) = eta_H * eps * B02 / rho
      positions(2, :) = [6, 8]
      call subblock( &
        quadblock, factors, positions, weight, h_quad, dh_cubic, settings%dims &
      )

      ! ==================== Cubic * Quadratic ====================
      call reset_factor_positions(new_size=2)
      ! H(7, 6)
      factors(1) = -eta_H * ic * B01 * eps * k3 / rho
      positions(1, :) = [7, 6]
      ! H(8, 6)
      factors(2) = eta_H * ic * B01 * k2 / rho
      positions(2, :) = [8, 6]
      call subblock( &
        quadblock, factors, positions, weight, h_cubic, h_quad, settings%dims &
      )

      ! ==================== Cubic * dCubic ====================
      call reset_factor_positions(new_size=2)
      ! H(7, 8)
      factors(1) = eta_H * ic * B01 * eps / rho
      positions(1, :) = [7, 8]
      ! H(8, 7)
      factors(2) = -eta_H * ic * B01 / rho
      positions(2, :) = [8, 7]
      call subblock( &
        quadblock, factors, positions, weight, h_cubic, dh_cubic, settings%dims &
      )
    end if

  end procedure add_natural_hall_terms

end submodule
