submodule (mod_boundary_manager:smod_natural_boundaries) smod_natural_bounds_regular
  implicit none

contains

  module procedure add_natural_regular_terms
    use mod_matrix_shortcuts, only: get_G_operator

    real(dp)  :: eps
    real(dp)  :: rho, T0
    real(dp)  :: B01, B02, B03
    real(dp)  :: Gop_min

    eps = eps_grid(grid_idx)
    rho = rho_field % rho0(grid_idx)
    T0 = T_field % T0(grid_idx)
    B01 = B_field % B01
    B02 = B_field % B02(grid_idx)
    B03 = B_field % B03(grid_idx)
    Gop_min = get_G_operator(grid_idx, which="minus")

    ! ==================== Cubic * Quadratic ====================
    call reset_factor_positions(new_size=3)
    ! A(2, 1)
    factors(1) = T0
    positions(1, :) = [2, 1]
    ! A(2, 5)
    factors(2) = rho
    positions(2, :) = [2, 5]
    ! A(2, 6)
    factors(3) = eps * Gop_min
    positions(3, :) = [2, 6]
    call subblock(quadblock, factors, positions, weight, h_cubic, h_quad)

    ! ==================== Cubic * dCubic ====================
    call reset_factor_positions(new_size=2)
    ! A(2, 7)
    factors(1) = B03
    positions(1, :) = [2, 7]
    ! A(2, 8)
    factors(2) = -eps * B02
    positions(2, :) = [2, 8]
    call subblock(quadblock, factors, positions, weight, h_cubic, dh_cubic)

    ! ==================== Quadratic * Quadratic ====================
    call reset_factor_positions(new_size=2)
    ! A(3, 6)
    factors(1) = ic * eps * k3 * B01
    positions(1, :) = [3, 6]
    ! A(4, 6)
    factors(2) = -ic * k2 * B01
    positions(2, :) = [4, 6]
    call subblock(quadblock, factors, positions, weight, h_quad, h_quad)

    ! ==================== Cubic * dQuadratic ====================
    call reset_factor_positions(new_size=1)
    ! A(3, 8)
    factors(1) = -ic * eps * B01
    positions(1, :) = [3, 8]
    call subblock(quadblock, factors, positions, weight, h_cubic, dh_quad)

    ! ==================== Quadratic * dcubic ====================
    call reset_factor_positions(new_size=1)
    ! A(4, 7)
    factors(1) = ic * B01
    positions(1, :) = [4, 7]
    call subblock(quadblock, factors, positions, weight, h_quad, dh_cubic)
  end procedure add_natural_regular_terms

  module procedure add_natural_regular_terms_vacuum
    use mod_global_variables, only: geometry, dp_LIMIT, flow, radiative_cooling, &
                                    external_gravity, thermal_conduction, resistivity, &
                                    viscosity, hall_mhd
    use mod_equilibrium_params, only: cte_Bv2, cte_Bv3

    real(dp)  :: x, eps, deps
    real(dp)  :: rho, T0
    real(dp)  :: B01, B02, B03, B0
    real(dp)  :: Fv, ktot, fac, pbalance
    real(dp)  :: BesKm, BesKm_1, BesKm1
    integer   :: nz, flag

    x = grid(grid_idx)
    eps = eps_grid(grid_idx)
    deps = d_eps_grid_dr(grid_idx)
    rho = rho_field % rho0(grid_idx)
    T0 = T_field % T0(grid_idx)
    B01 = B_field % B01
    B02 = B_field % B02(grid_idx)
    B03 = B_field % B03(grid_idx)
    B0  = B_field % B0(grid_idx)
    Fv = k2 * cte_Bv2 / eps + k3 * cte_Bv3
    ktot = sqrt(k2**2 + k3**2)

    if (B01 > dp_LIMIT) then
      call log_message('B01 cannot be non-zero for a plasma-vacuum interface', level='error')
    end if

    !> At the interface, there should be total pressure balance on both sides
    pbalance = rho * T0 + (B0**2 - (cte_Bv2**2 + cte_Bv3**2)) / 2.0d0
    if (pbalance > dp_LIMIT) then
      call log_message('pressure balance not satisfied at the interface !', level='warning')
      call log_message('discrepancy = ' // str(pbalance), level='warning')
    end if

    if (flow .or. radiative_cooling .or. external_gravity .or. thermal_conduction &
        .or. resistivity .or. viscosity .or. hall_mhd) then
      call log_message('flow, conductivity, resistivity, radiative cooling, viscosity, &
                        &external gravity and Hall effects are currently not supported &
                        &for a plasma-vacuum transition', level='error')
    end if

    if (geometry == 'Cartesian') then
      ! ==================== Cubic * Cubic ====================
      call reset_factor_positions(new_size=2)
      ! A(2, 7)
      factors(1) = -eps * Fv**2 * B03 / (ktot * B0**2)
      positions(1, :) = [2, 7]
      ! A(2, 8)
      factors(2) = eps * Fv**2 * B02 / (ktot * B0**2)
      positions(2, :) = [2, 8]
      call subblock(quadblock, factors, positions, weight, h_cubic, h_cubic)

    else if (geometry == 'cylindrical') then
      if (abs(k3) > dp_LIMIT) then
        !> cbesk computes the Bessel K_m(|k3|R) value
        !> It is defined in the amos library
        call cbesk(abs(k3) * x, int(k2), 1, 1, BesKm, nz, flag)
        call cbesk(abs(k3) * x, int(k2)-1, 1, 1, BesKm_1, nz, flag)
        call cbesk(abs(k3) * x, int(k2)+1, 1, 1, BesKm1, nz, flag)
        fac = Fv**2 * BesKm / (0.5d0 * abs(k3) * (BesKm_1 + BesKm1)) + (B02**2 - cte_Bv2**2) / x
      else if (abs(k2) > dp_LIMIT) then
        fac = (k2 * cte_Bv2)**2 / (abs(k2) * x) + (B02**2 - cte_Bv2**2) / x
      else
        fac = (B02**2 - cte_Bv2**2) / x
      end if

      ! ==================== Cubic * Cubic ====================
      call reset_factor_positions(new_size=2)
      ! A(2, 7)
      factors(1) = -eps * B03 * fac / B0**2
      positions(1, :) = [2, 7]
      ! A(2, 8)
      factors(2) = eps * B02 * fac / B0**2
      positions(2, :) = [2, 8]
      call subblock(quadblock, factors, positions, weight, h_cubic, h_cubic)
    end if
  end procedure add_natural_regular_terms_vacuum

end submodule smod_natural_bounds_regular
