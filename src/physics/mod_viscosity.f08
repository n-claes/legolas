! =============================================================================
!> This module is responsible for calculating the viscosity contributions
!! to the A-matrix based on the equilibrium configuration.
module mod_viscosity
  use mod_global_variables, only: dp, dim_quadblock, viscosity_value

  implicit none

  private

  public :: get_viscosity_terms
  public :: viscosity_boundaries

contains

  subroutine get_viscosity_terms(gauss_idx, eps, d_eps_dr, curr_weight, quadblock_viscosity, &
                            h_quadratic, dh_quadratic_dr, h_cubic, dh_cubic_dr)
    use mod_global_variables, only: viscous_heating
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

    real(dp)                  :: eps_inv
    real(dp)                  :: drho0, v02, dv02, ddv02, v03, dv03, ddv03

    quadblock_viscosity = (0.0d0, 0.0d0)
    eps_inv = 1.0d0 / eps

    ! Equilibrium quantities
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
    ! A(3, 3)
    factors(1) = -((k2 * eps_inv)**2 + k3**2) - (k2 * eps_inv)**2 / 3.0d0 - d_eps_dr * eps_inv**2
    positions(1, :) = [3, 3]
    ! A(3, 4)
    factors(2) = -k2 * k3 * eps_inv**2 / 3.0d0
    positions(2, :) = [3, 4]
    ! A(4, 3)
    factors(3) = -k2 * k3 * eps_inv / 3.0d0
    positions(3, :) = [4, 3]
    ! A(4, 4)
    factors(4) = -eps_inv * ((k2 * eps_inv)**2 + k3**2) - k3**2 * eps_inv / 3.0d0 &
                 -d_eps_dr**2 * eps_inv**3
    positions(4, :) = [4, 4]
    ! A(5, 3)
    factors(5) = 0.0d0
    positions(5, :) = [5, 3]
    ! A(5, 4)
    factors(6) = 0.0d0
    positions(6, :) = [5, 4]
    if (viscous_heating) then
      factors(5) = factors(5) + 2.0d0 * (d_eps_dr * eps_inv)**2 * v02
      factors(6) = factors(6) - 2.0d0 * eps_inv * ddv03
    end if
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_quadratic, h_quadratic)

    ! Quadratic * d(Quadratic)/dr
    call reset_factors(factors, 3)
    call reset_positions(positions, 3)
    ! A(3,3)
    factors(1) = d_eps_dr * eps_inv
    positions(1, :) = [3, 3]
    ! A(4,4)
    factors(2) = d_eps_dr * eps_inv**2
    positions(2, :) = [4, 4]
    ! A(5,3)
    factors(3) = 0.0d0
    positions(3, :) = [5, 3]
    if (viscous_heating) then
      factors(3) = factors(3) + 2.0d0 * dv02
    end if
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_quadratic, dh_quadratic_dr)

    ! d(Quadratic)/dr * Quadratic
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(4, 4)
    factors(1) = d_eps_dr * eps_inv**2
    positions(1, :) = [4, 4]
    ! A(5, 4)
    factors(2) = 0.0d0
    positions(2, :) = [5, 4]
    if (viscous_heating) then
      factors(2) = factors(2) - 2.0d0 * eps_inv * dv03
    end if
    call subblock(quadblock_viscosity, factors, positions, curr_weight, dh_quadratic_dr, h_quadratic)

    ! d(Quadratic)/dr * d(Quadratic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(3,3)
    factors(1) = -1.0d0
    positions(1, :) = [3, 3]
    ! A(4,4)
    factors(2) = -eps_inv
    positions(2, :) = [4, 4]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, dh_quadratic_dr, dh_quadratic_dr)

    if (viscous_heating) then
      ! Quadratic * Cubic
      call reset_factors(factors, 1)
      call reset_positions(positions, 1)
      ! A(5,2)
      factors(1) = -2.0d0 * d_eps_dr * eps_inv**3 * k2 * v02
      positions(1, :) = [5, 2]
      call subblock(quadblock_viscosity, factors, positions, curr_weight, h_quadratic, h_cubic)
    end if

    ! Quadratic * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(3,2)
    factors(1) = k2 * eps_inv**2 / 3.0d0
    positions(1, :) = [3, 2]
    ! A(4,2)
    factors(2) = k3 * eps_inv / 3.0d0
    positions(2, :) = [4, 2]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_quadratic, dh_cubic_dr)

    ! Cubic * Quadratic
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(2,3)
    factors(1) = 2.0d0 * d_eps_dr * eps_inv**2 * k2
    positions(1, :) = [2, 3]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_cubic, h_quadratic)

    ! d(Cubic)/dr * Quadratic
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(2,3)
    factors(1) = k2 * eps_inv / 3.0d0
    positions(1, :) = [2, 3]
    ! A(2,4)
    factors(2) = k3 * eps_inv / 3.0d0
    positions(2, :) = [2, 4]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, dh_cubic_dr, h_quadratic)

    ! Cubic * Cubic
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(2,2)
    factors(1) = -eps_inv * ((k2 * eps_inv)**2 + k3**2) - d_eps_dr**2 * eps_inv**3 - d_eps_dr * eps_inv**3
    positions(1, :) = [2, 2]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_cubic, h_cubic)

    ! d(Cubic)/dr * Cubic
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(2,2)
    factors(1) = d_eps_dr * eps_inv**2
    positions(1, :) = [2, 2]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, dh_cubic_dr, h_cubic)

    ! Cubic * d(Cubic)/dr
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(2,2)
    factors(1) = d_eps_dr * eps_inv**2
    positions(1, :) = [2, 2]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, h_cubic, dh_cubic_dr)

    ! d(Cubic)/dr * d(Cubic)/dr
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(2,2)
    factors(1) = -4.0d0 * eps_inv / 3.0d0
    positions(1, :) = [2, 2]
    call subblock(quadblock_viscosity, factors, positions, curr_weight, dh_cubic_dr, dh_cubic_dr)

    deallocate(factors)
    deallocate(positions)
  end subroutine get_viscosity_terms

  !> Creates a quadblock for the A matrix containing the natural boundary
  !! conditions coming from the viscosity terms, depending on the supplied edge.
  subroutine viscosity_boundaries(quadblock_viscosity, edge)
    use mod_global_variables, only: dp_LIMIT, gauss_gridpts, geometry, coaxial, viscous_heating
    use mod_logging, only: log_message
    use mod_equilibrium, only: v_field, kappa_field
    use mod_grid, only: eps_grid, d_eps_grid_dr
    use mod_make_subblock, only: reset_factors, reset_positions

    !> the quadblock corresponding to the left/right edge
    complex(dp), intent(out)      :: quadblock_viscosity(dim_quadblock, dim_quadblock)
    !> the edge, either "l_edge" or "r_edge"
    character(len=6), intent(in)  :: edge

    complex(dp), allocatable  :: surface_terms(:)
    real(dp)                  :: eps, eps_inv, d_eps_dr, dv03
    integer, allocatable      :: positions(:, :)
    integer                   :: idx, i
    logical, save             :: kappa_perp_is_zero

    quadblock_viscosity = (0.0d0, 0.0d0)
    kappa_perp_is_zero = .true.
    if (any(abs(kappa_field % kappa_perp) > dp_LIMIT)) then
      kappa_perp_is_zero = .false.
    end if

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

    dv03 = v_field % d_v03_dr(idx)

    if (geometry == 'cylindrical' .and. edge == 'l_edge' .and. .not. coaxial) then
      ! Add to bottom-right of 2x2 block, for top-left subblock only
      ! (Quadratic * Quadratic, Quadratic * d(Quadratic)/dr)
      call reset_factors(surface_terms, 2)
      call reset_positions(positions, 2)

      ! surface term for element (4, 4)
      surface_terms(1) = -d_eps_dr * eps_inv**2
      positions(1, :) = [4, 4]
      ! surface term for element (5, 4)
      if (kappa_perp_is_zero .and. viscous_heating) then
        surface_terms(2) = 2.0d0 * eps_inv * dv03
      end if
      positions(2, :) = [5, 4]

      positions = 2 * positions
      do i = 1, size(surface_terms)
        quadblock_viscosity(positions(i, 1), positions(i, 2)) = &
        quadblock_viscosity(positions(i, 1), positions(i, 2)) + surface_terms(i)
      end do

      deallocate(surface_terms)
      deallocate(positions)
    end if
  end subroutine viscosity_boundaries
end module mod_viscosity
