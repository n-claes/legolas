! =============================================================================
!> This module handles everything related to the Hall term, computing
!! the normalised Hall factor, the contributions to the A matrix, and natural
!! boundary conditions.
module mod_hallmhd
  use mod_global_variables, only: dp, dim_quadblock

implicit none

private

public  :: get_hallfactor
public  :: get_hallterm
public  :: hall_boundaries

contains

  !> Retrieves the normalised Hall factor as described by Porth et al. (2014)
  subroutine get_hallfactor()
    use mod_global_variables, only: cgs_units, hallfactor
    use mod_physical_constants, only: mp_cgs, mp_si, ec_cgs, ec_si
    use mod_units, only: unit_velocity, unit_length, unit_magneticfield

    if (cgs_units) then
      hallfactor = (mp_cgs * unit_velocity) / (ec_cgs * unit_length * unit_magneticfield)
    else
      hallfactor = (mp_si * unit_velocity) / (ec_si * unit_length * unit_magneticfield)
    end if
  end subroutine get_hallfactor

  !> Retrieves the A-matrix Hall contribution for a given quadblock corresponding
  !! to a particular point in the grid and a particular Gaussian weight.
  !! If <tt>hall_mhd = True</tt>, this routine is called <tt>n_gauss</tt> times
  !! for every grid interval.
  subroutine get_hallterm(gauss_idx, eps, d_eps_dr, curr_weight, quadblock_Hall, &
                          h_quadratic, h_cubic, dh_cubic_dr)
    use mod_equilibrium_params, only: k2, k3
    use mod_equilibrium, only: rho_field, B_field
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
    complex(dp), intent(out)  :: quadblock_Hall(dim_quadblock, dim_quadblock)
    !> array containing the 4 quadratic basis functions
    real(dp), intent(in)      :: h_quadratic(4)
    !> array containing the 4 cubic basis functions
    real(dp), intent(in)      :: h_cubic(4)
    !> array containing the derivatives of the 4 cubic basis functions
    real(dp), intent(in)      :: dh_cubic_dr(4)

    complex(dp), allocatable  :: factors(:)
    integer, allocatable      :: positions(:, :)

    real(dp)                  :: eps_inv
    real(dp)                  :: rho0, B02, B03
    real(dp)                  :: drho0, drB02, dB02_r, dB03, dB02

    quadblock_Hall = (0.0d0, 0.0d0)
    eps_inv = 1.0d0 / eps

    !! Equilibrium quantities
    rho0 = rho_field % rho0(gauss_idx)
    drho0 = rho_field % d_rho0_dr(gauss_idx)

    B02 = B_field % B02(gauss_idx)
    dB02 = B_field % d_B02_dr(gauss_idx)
    B03 = B_field % B03(gauss_idx)
    dB03   = B_field % d_B03_dr(gauss_idx)

    !! Calculate derivatives eps*B02, B02/eps
    drB02 = d_eps_dr * B02 + eps * dB02
    dB02_r = eps_inv * dB02 - d_eps_dr * eps_inv**2 * B02

    !! Setup of matrix elements

    ! Quadratic * Quadratic
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(6, 6)
    factors(1) = - (eps_inv * k2 * B02 + k3 * B03) * ((k2 * eps_inv)**2 + k3**2) / rho0
    positions(1, :) = [6, 6]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_quadratic, h_quadratic)

    ! Quadratic * Cubic
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(6, 7)
    factors(1) = k3 * eps_inv**2 * (k3 * drB02 - k2 * dB03) / rho0
    positions(1, :) = [6, 7]
    ! A(6, 8)
    factors(2) = - k2 * eps_inv**2 * (k3 * drB02 - k2 * dB03) / rho0
    positions(2, :) = [6, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_quadratic, h_cubic)

    ! Quadratic * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(6, 7)
    factors(1) = k2 * eps_inv**2 * (k2 * B02 * eps_inv + k3 * B03) / rho0
    positions(1, :) = [6, 7]
    ! A(6, 8)
    factors(2) = k3 * (k2 * B02 * eps_inv + k3 * B03) / rho0
    positions(2, :) = [6, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_quadratic, dh_cubic_dr)

    ! Cubic * Quadratic
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(7, 1)
    factors(1) = - k3 * eps_inv * (B02 * eps_inv * drB02 + B03 * dB03) / rho0**2
    positions(1, :) = [7, 1]
    ! A(7, 6)
    factors(2) = k3**2 * (drB02 * eps_inv - eps * dB02_r) / rho0 &
                  + k3 * drho0 * (k3 * B02 - k2 * B03 * eps_inv) / rho0**2
    positions(2, :) = [7, 6]
    ! A(8, 1)
    factors(3) = k2 * eps_inv**2 * (B02 * eps_inv * drB02 + B03 * dB03) / rho0**2
    positions(3, :) = [8, 1]
    ! A(8, 6)
    factors(4) = k2 * k3 * eps_inv * (dB02 - drB02 * eps_inv - drho0 * B02 / rho0 &
                  - 2.0d0 * d_eps_dr * B02 * eps_inv) / rho0 &
                  + B03 * (drho0 * (k2 * eps_inv)**2 / rho0 - d_eps_dr * eps_inv * k3**2) / rho0
    positions(4, :) = [8, 6]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_cubic, h_quadratic)

    ! d(Cubic)/dr * Quadratic
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(7, 6)
    factors(1) = k2 * eps_inv * (k2 * eps_inv * B02 + k3 * B03) / rho0
    positions(1, :) = [7, 6]
    ! A(8, 6)
    factors(2) = k3 * (k2 * eps_inv * B02 + k3 * B03) / rho0
    positions(2, :) = [8, 6]
    call subblock(quadblock_Hall, factors, positions, curr_weight, dh_cubic_dr, h_quadratic)

    ! Cubic * Cubic
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(7, 7)
    factors(1) = - k3**2 * eps_inv * (k2 * eps_inv * B02 + k3 * B03) / rho0
    positions(1, :) = [7, 7]
    ! A(7, 8)
    factors(2) = k2 * k3 * eps_inv * (k2 * eps_inv * B02 + k3 * B03) / rho0
    positions(2, :) = [7, 8]
    ! A(8, 7)
    factors(3) = k2 * k3 * eps_inv**2 * (k2 * eps_inv * B02 + k3 * B03) / rho0 &
                  + d_eps_dr * eps_inv**3 * k3 * drB02 / rho0
    positions(3, :) = [8, 7]
    ! A(8, 8)
    factors(4) = - (k2 * eps_inv)**2 * (k2 * eps_inv * B02 + k3 * B03) / rho0 &
                  - d_eps_dr * eps_inv**3 * k2 * drB02 / rho0
    positions(4, :) = [8, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_cubic, h_cubic)

    ! Cubic * d(Cubic)/dr
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(7, 7)
    factors(1) = k3 * B03 * eps_inv * drho0 / rho0**2
    positions(1, :) = [7, 7]
    ! A(7, 8)
    factors(2) = k3 * (eps * dB02_r - drB02 * eps_inv - B02 * drho0 / rho0) / rho0
    positions(2, :) = [7, 8]
    ! A(8, 7)
    factors(3) = - k2 * B03 * eps_inv**2 * drho0 / rho0**2
    positions(3, :) = [8, 7]
    ! A(8, 8)
    factors(4) = k2 * eps_inv * (2.0d0 * d_eps_dr * eps_inv * B02 &
                  + drB02 * eps_inv - dB02 + B02 * drho0 / rho0) / rho0 &
                  + k3 * B03 * d_eps_dr * eps_inv / rho0
    positions(4, :) = [8, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_cubic, dh_cubic_dr)

    ! d(Cubic)/dr * Cubic
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(7, 7)
    factors(1) = k3 * dB03 * eps_inv / rho0
    positions(1, :) = [7, 7]
    ! A(7, 8)
    factors(2) = - k2 * dB03 * eps_inv / rho0
    positions(2, :) = [7, 8]
    ! A(8, 7)
    factors(3) = - k3 * drB02 * eps_inv**2 / rho0
    positions(3, :) = [8, 7]
    ! A(8, 8)
    factors(4) = k2 * drB02 * eps_inv**2 / rho0
    positions(4, :) = [8, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, dh_cubic_dr, h_cubic)

    ! d(Cubic)/dr * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(7, 7)
    factors(1) = - eps_inv * (k2 * B02 * eps_inv + k3 * B03) / rho0
    positions(1, :) = [7, 7]
    ! A(8, 8)
    factors(2) = - (k2 * B02 * eps_inv + k3 * B03) / rho0
    positions(2, :) = [8, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, dh_cubic_dr, dh_cubic_dr)
  end subroutine get_hallterm

  !> Creates a quadblock for the A matrix containing the natural boundary
  !! conditions coming from the Hall term.
  !! @note    This subroutine is currently unused since \(a2 = a3 = 0\) for a
  !!          fixed wall boundary. For potential, different boundary types, the
  !!          quadblock_Hall created by this subroutine should be added in the
  !!          <tt>apply_boundary_conditions</tt> function to the existinng
  !!          quadblock from the essential and natural boundary conditions as
  !!          <tt>quadblock = quadblock - hallfactor * quadblock_Hall</tt>
  !!          for both the left and right edge. @endnote
  subroutine hall_boundaries(quadblock_Hall, edge)
    use mod_global_variables, only: gauss_gridpts, dim_subblock
    use mod_logging, only: log_message
    use mod_equilibrium, only: rho_field, B_field
    use mod_grid, only: eps_grid, d_eps_grid_dr
    use mod_equilibrium_params, only: k2, k3
    use mod_make_subblock, only: reset_factors, reset_positions

    !> the quadblock corresponding to the left/right edge
    complex(dp), intent(out)      :: quadblock_Hall(dim_quadblock, dim_quadblock)
    !> the edge, either "l_edge" or "r_edge"
    character(len=6), intent(in)  :: edge

    complex(dp), allocatable  :: surface_terms(:)
    real(dp)                  :: eps, d_eps_dr, eps_inv
    real(dp)                  :: rho0, B02, dB02, B03, dB03, drB02
    integer, allocatable      :: positions(:, :)
    integer                   :: idx, i

    quadblock_Hall = (0.0d0, 0.0d0)

    if (edge == 'l_edge') then
      idx = 1
    else if (edge == 'r_edge') then
      idx = gauss_gridpts
    else
      call log_message('Hall boundaries: wrong edge supplied' // edge, level='error')
    end if

    ! retrieve variables at current edge
    eps = eps_grid(idx)
    eps_inv = 1.0d0 / eps
    d_eps_dr = d_eps_grid_dr(idx)
    rho0 = rho_field % rho0(idx)
    B02 = B_field % B02(idx)
    dB02 = B_field % d_B02_dr(idx)
    B03 = B_field % B03(idx)
    dB03 = B_field % d_B03_dr(idx)
    drB02 = d_eps_dr * B02 + eps * dB02

    ! Cubic * Quadratic
    call reset_factors(surface_terms, 2)
    call reset_positions(positions, 2)

    ! surface term for element (7, 6)
    surface_terms(1) = - k2 * eps_inv * (k2 * eps_inv * B02 + k3 * B03) / rho0
    positions(1, :) = [7, 6]
    ! surface term for element (8, 6)
    surface_terms(2) = - k3 * (k2 * eps_inv * B02 + k3 * B03) / rho0
    positions(2, :) = [8, 6]

    ! l_edge: add to top-right of 2x2 block, for top-left subblock only
    ! r_edge: add to top-right of 2x2 block, for bottom-right subblock only
    if (edge == 'l_edge') then
      positions = 2 * positions - spread([1, 0], 1, 2)
    else if (edge == 'r_edge') then
      positions = 2 * positions - spread([1, 0], 1, 2) + dim_subblock
    end if

    do i = 1, size(surface_terms)
      quadblock_Hall(positions(i, 1), positions(i, 2)) = quadblock_Hall(positions(i, 1), positions(i, 2)) + surface_terms(i)
    end do

    ! Cubic * Cubic
    call reset_factors(surface_terms, 4)
    call reset_positions(positions, 4)

    ! surface term for element (7, 7)
    surface_terms(1) = - k3 * eps_inv * dB03 / rho0
    positions(1, :) = [7, 7]
    ! surface term for element (7, 8)
    surface_terms(2) = k2 * dB03 * eps_inv / rho0
    positions(2, :) = [7, 8]
    ! surface term for element (8, 7)
    surface_terms(3) = k3 * drB02 * eps_inv**2 / rho0
    positions(3, :) = [8, 7]
    ! surface term for element (8, 8)
    surface_terms(4) = - k2 * eps_inv**2 * drB02 / rho0
    positions(4, :) = [8, 8]

    ! l_edge: add to top-left of 2x2 block, for top-left subblock only
    ! r_edge: add to top-left of 2x2 block, for bottom-right subblock only
    if (edge == 'l_edge') then
      positions = 2 * positions - spread([1, 1], 1, 4)
    else if (edge == 'r_edge') then
      positions = 2 * positions - spread([1, 1], 1, 4) + dim_subblock
    end if

    do i = 1, size(surface_terms)
      quadblock_Hall(positions(i, 1), positions(i, 2)) = quadblock_Hall(positions(i, 1), positions(i, 2)) + surface_terms(i)
    end do

    ! Cubic * d(Cubic)/dr
    call reset_factors(surface_terms, 2)
    call reset_positions(positions, 2)

    ! surface term for element (7, 7)
    surface_terms(1) = eps_inv * (k2 * B02 * eps_inv + k3 * B03) / rho0
    positions(1, :) = [7, 7]
    ! surface term for element (8, 8)
    surface_terms(2) = (k2 * B02 * eps_inv + k3 * B03) / rho0
    positions(2, :) = [8, 8]

    ! l_edge: add to top-right of 2x2 block, for top-left subblock only
    ! r_edge: add to top-right of 2x2 block, for bottom-right subblock only
    if (edge == 'l_edge') then
      positions = 2 * positions - spread([1, 0], 1, 2)
    else if (edge == 'r_edge') then
      positions = 2 * positions - spread([1, 0], 1, 2) + dim_subblock
    end if

    do i = 1, size(surface_terms)
      quadblock_Hall(positions(i, 1), positions(i, 2)) = quadblock_Hall(positions(i, 1), positions(i, 2)) + surface_terms(i)
    end do

    deallocate(positions)
    deallocate(surface_terms)
  end subroutine hall_boundaries

end module mod_hallmhd
