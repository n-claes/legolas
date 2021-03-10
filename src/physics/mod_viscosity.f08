! =============================================================================
!> This module is responsible for calculating the viscosity contributions
!! to the A-matrix based on the equilibrium configuration.
module mod_viscosity
  use mod_global_variables, only: dp, dim_quadblock, viscosity_value, use_kinematic_viscosity

  implicit none

  private

  public :: get_viscosity_terms
  public :: viscosity_boundaries

contains

  subroutine get_viscosity_terms(gauss_idx, eps, d_eps_dr, curr_weight, quadblock_viscosity, &
                            h_quadratic, dh_quadratic_dr, h_cubic, dh_cubic_dr)
    use mod_equilibrium_params, only: k2, k3
    use mod_equilibrium, only: rho_field, v_field, viscosity_field
    use mod_make_subblock, only: subblock, reset_factors, reset_positions

    !> current index in the Gaussian grid
    integer, intent(in)       :: gauss_idx
    !> current value for the scale factor epsilon
    real(dp), intent(in)      :: eps
    !> current value for the derivative of epsilon
    real(dp), intent(in)      :: d_eps_dr
    !> current weight in the Gaussian quadrature
    real(dp), intent(in)      :: curr_weight
    !> the A-quadblock for a particular grid interval
    complex(dp), intent(out)  :: quadblock_viscosity(dim_quadblock, dim_quadblock)
    !> array containing the 4 quadratic basis functions
    real(dp), intent(in)      :: h_quadratic(4)
    !> array containing the derivatives of the 4 quadratic basis functions
    real(dp), intent(in)      :: dh_quadratic_dr(4)
    !> array containing the 4 cubic basis functions
    real(dp), intent(in)      :: h_cubic(4)
    !> array containing the derivatives of the 4 cubic basis functions
    real(dp), intent(in)      :: dh_cubic_dr(4)

    complex(dp), allocatable  :: factors(:)
    integer, allocatable      :: positions(:, :)

    real(dp)                  :: eps_inv, zeta
    real(dp)                  :: drho0, v02, dv02, ddv02, v03, dv03, ddv03

    quadblock_viscosity = (0.0d0, 0.0d0)
    eps_inv = 1.0d0 / eps

    ! Equilibrium quantities
    zeta = 1.0d0
    if (use_kinematic_viscosity) then
      zeta = rho_field % rho0(gauss_idx)
    end if
    drho0 = rho_field % d_rho0_dr(gauss_idx)

    v02 = v_field % v02(gauss_idx)
    dv02 = v_field % d_v02_dr(gauss_idx)
    v03 = v_field % v03(gauss_idx)
    dv03 = v_field % d_v03_dr(gauss_idx)

    ddv02 = viscosity_field % dd_v02_dr(gauss_idx)
    ddv03 = viscosity_field % dd_v03_dr(gauss_idx)

    ! Setup of matrix elements

    ! Quadratic * Quadratic
    call reset_factors(factors, 6)
    call reset_positions(positions, 6)
    ! A(3, 1)
    factors(1) = 0.0d0
    positions(1, :) = [3, 1]
    ! A(3, 3)
    factors(2) = -zeta * ((k2 * eps_inv)**2 + k3**2) - zeta * (k2 *eps_inv)**2 / 3.0d0
    positions(2, :) = [3, 3]
    ! A(3, 4)
    factors(3) = -zeta * k2 * k3 *eps_inv**2 / 3.0d0
    positions(3, :) = [3, 4]
    ! A(4, 1)
    factors(4) = 0.0d0
    positions(4, :) = [4, 1]
    ! A(4, 3)
    factors(5) = -zeta * k2 * k3 * eps_inv / 3.0d0
    positions(5, :) = [4, 3]
    ! A(4, 4)
    factors(6) = -zeta * eps_inv * ((k2 * eps_inv)**2 + k3**2) - zeta * k3**2 * eps_inv / 3.0d0 &
                 -zeta * d_eps_dr**2 * eps_inv**3
    positions(6, :) = [4, 4]
    if (use_kinematic_viscosity) then
      factors(1) = factors(1) + eps_inv**2 * (d_eps_dr * dv02 + eps * ddv02)
      factors(4) = factors(4) + eps_inv**2 * (d_eps_dr * dv03 + eps * ddv03)
      factors(6) = factors(6) + drho0 * d_eps_dr * eps_inv**2
    end if
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_quadratic, h_quadratic)

    ! Quadratic * d(Quadratic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(3,3)
    factors(1) = zeta * d_eps_dr * eps_inv
    positions(1, :) = [3, 3]
    ! A(4,4)
    factors(2) = zeta * d_eps_dr * eps_inv**2
    positions(2, :) = [4, 4]
    if (use_kinematic_viscosity) then
      factors(1) = factors(1) - drho0
      factors(2) = factors(2) - drho0 * eps_inv
    end if
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_quadratic, dh_quadratic_dr)

    ! d(Quadratic)/dr * Quadratic
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(4, 4)
    factors(1) = zeta * d_eps_dr * eps_inv**2
    positions(1, :) = [4, 4]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, dh_quadratic_dr, h_quadratic)

    ! d(Quadratic)/dr * d(Quadratic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(3,3)
    factors(1) = -zeta
    positions(1, :) = [3, 3]
    ! A(4,4)
    factors(2) = -zeta * eps_inv
    positions(2, :) = [4, 4]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, dh_quadratic_dr, dh_quadratic_dr)

    ! Quadratic * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(3,2)
    factors(1) = zeta * k2 * eps_inv**2 / 3.0d0
    positions(1, :) = [3, 2]
    ! A(4,2)
    factors(2) = zeta * k3 * eps_inv / 3.0d0
    positions(2, :) = [4, 2]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_quadratic, dh_cubic_dr)

    ! Cubic * Quadratic
    if (use_kinematic_viscosity) then
      call reset_factors(factors, 2)
      call reset_positions(positions, 2)
      ! A(2,3)
      factors(1) = drho0 * k2 * eps_inv / 3.0d0
      positions(1, :) = [2, 3]
      ! A(2,4)
      factors(2) = drho0 * k3 * eps_inv / 3.0d0
      positions(2, :) = [2, 4]
      call subblock(quadblock_viscosity, factors, positions, curr_weight, h_cubic, h_quadratic)
    end if

    ! d(Cubic)/dr * Quadratic
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(2,3)
    factors(1) = zeta * k2 * eps_inv / 3.0d0
    positions(1, :) = [2, 3]
    ! A(2,4)
    factors(2) = zeta * k3 * eps_inv / 3.0d0
    positions(2, :) = [2, 4]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, dh_cubic_dr, h_quadratic)

    ! Cubic * Cubic
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(2,2)
    factors(1) = -zeta * eps_inv * ((k2 * eps_inv)**2 + k3**2) - zeta * d_eps_dr**2 * eps_inv**3
    positions(1, :) = [2, 2]
    if (use_kinematic_viscosity) then
      factors(1) = factors(1) + drho0 * d_eps_dr * eps_inv**2
    end if
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_cubic, h_cubic)

    ! d(Cubic)/dr * Cubic
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(2,2)
    factors(1) = zeta * d_eps_dr * eps_inv**2
    positions(1, :) = [2, 2]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, dh_cubic_dr, h_cubic)

    ! Cubic * d(Cubic)/dr
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(2,2)
    factors(1) = zeta * d_eps_dr * eps_inv**2
    positions(1, :) = [2, 2]
    if (use_kinematic_viscosity) then
      factors(1) = factors(1) - 4.0d0 * drho0 * eps_inv / 3.0d0
    end if
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_cubic, dh_cubic_dr)

    ! d(Cubic)/dr * d(Cubic)/dr
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(2,2)
    factors(1) = -4.0d0 * zeta * eps_inv / 3.0d0
    positions(1, :) = [2, 2]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, dh_cubic_dr, dh_cubic_dr)

    deallocate(factors)
    deallocate(positions)
  end subroutine get_viscosity_terms

  !> Creates a quadblock for the A matrix containing the natural boundary
  !! conditions coming from the viscosity terms, depending on the supplied edge.
  subroutine viscosity_boundaries(quadblock_viscosity, edge)
    use mod_global_variables, only: gauss_gridpts, dim_subblock, geometry, &
                                    use_kinematic_viscosity
    use mod_logging, only: log_message
    use mod_equilibrium, only: rho_field
    use mod_grid, only: eps_grid, d_eps_grid_dr
    use mod_make_subblock, only: reset_factors, reset_positions

    !> the quadblock corresponding to the left/right edge
    complex(dp), intent(out)      :: quadblock_viscosity(dim_quadblock, dim_quadblock)
    !> the edge, either "l_edge" or "r_edge"
    character(len=6), intent(in)  :: edge

    complex(dp), allocatable  :: surface_terms(:)
    real(dp)                  :: eps, eps_inv, d_eps_dr
    real(dp)                  :: zeta
    integer, allocatable      :: positions(:, :)
    integer                   :: idx, i

    quadblock_viscosity = (0.0d0, 0.0d0)

    if (edge == 'l_edge') then
      idx = 1
    else if (edge == 'r_edge') then
      idx = gauss_gridpts
    else
      call log_message('Viscous boundaries: wrong edge supplied' // edge, level='error')
    end if

    ! retrieve variables at current edge
    eps = eps_grid(idx)
    eps_inv = 1.0d0 / eps
    d_eps_dr = d_eps_grid_dr(idx)
    zeta = 1.0d0
    if (use_kinematic_viscosity) then
      zeta = rho_field % rho0(idx)
    end if

    if (geometry == 'cylindrical' .and. edge == 'l_edge') then
      ! Add to bottom-right of 2x2 block, for top-left subblock
      ! (Quadratic * d(Quadratic)/dr, Quadratic * Quadratic)
      call reset_factors(surface_terms, 2)
      call reset_positions(positions, 2)

      ! surface term for element (3, 3)
      surface_terms(1) = -3.0d0 * zeta
      positions(1, :) = [3, 3]
      ! surface term for element (4, 4)
      surface_terms(2) = -3.0d0 * zeta * eps_inv - zeta * d_eps_dr * eps_inv**2
      positions(2, :) = [4, 4]

      positions = 2 * positions
      do i = 1, size(surface_terms)
        quadblock_viscosity(positions(i, 1), positions(i, 2)) = &
        quadblock_viscosity(positions(i, 1), positions(i, 2)) + surface_terms(i)
      end do

      ! Add to bottom-left of 2x2 block, for top-right subblock
      ! (Quadratic * d(Quadratic)/dr)
      call reset_factors(surface_terms, 2)
      call reset_positions(positions, 2)

      ! surface term for element (3, 3)
      surface_terms(1) = 4.0d0 * zeta
      positions(1, :) = [3, 3]
      ! surface term for element (4, 4)
      surface_terms(2) = 4.0d0 * zeta * eps_inv
      positions(2, :) = [4, 4]

      positions = 2 * positions + spread([0, dim_subblock-1], 1, size(surface_terms))
      do i = 1, size(surface_terms)
        quadblock_viscosity(positions(i, 1), positions(i, 2)) = &
        quadblock_viscosity(positions(i, 1), positions(i, 2)) + surface_terms(i)
      end do

      ! Add to bottom-right of 2x2 block, for top-right subblock
      ! (Quadratic * d(Quadratic)/dr)
      call reset_factors(surface_terms, 2)
      call reset_positions(positions, 2)

      ! surface term for element (3, 3)
      surface_terms(1) = -zeta
      positions(1, :) = [3, 3]
      ! surface term for element (4, 4)
      surface_terms(2) = -zeta * eps_inv
      positions(2, :) = [4, 4]

      positions = 2 * positions + spread([0, dim_subblock], 1, 2)
      do i = 1, size(surface_terms)
        quadblock_viscosity(positions(i, 1), positions(i, 2)) = &
        quadblock_viscosity(positions(i, 1), positions(i, 2)) + surface_terms(i)
      end do

      deallocate(surface_terms)
      deallocate(positions)
    end if
  end subroutine viscosity_boundaries
end module mod_viscosity
