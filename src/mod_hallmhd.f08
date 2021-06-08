! =============================================================================
!> This module handles everything related to the Hall term, computing
!! the normalised Hall factor, the contributions to the A matrix, and natural
!! boundary conditions, when necessary.
module mod_hallmhd
  use mod_global_variables, only: dp, dim_quadblock

implicit none

private

public  :: set_hallfactor
public  :: get_hallterm
public  :: get_hallterm_Bmat
public  :: hall_boundaries

contains

  !> Retrieves the normalised Hall factor as described by Porth et al. (2014),
  !! with a dropoff at the boundary
  subroutine set_hallfactor(hall_field)
    use mod_grid, only: grid_gauss
    use mod_physical_constants, only: dpi
    use mod_global_variables, only: cgs_units, gauss_gridpts, dropoff_edge_dist, &
                                    dropoff_width, hall_dropoff, inertia_dropoff
    use mod_physical_constants, only: mp_cgs, mp_si, ec_cgs, ec_si, me_cgs, me_si
    use mod_units, only: unit_velocity, unit_length, unit_magneticfield
    use mod_types, only: hall_type

    type (hall_type), intent(inout)  :: hall_field

    real(dp)  :: sleft, sright, width, hallval, inertiaval
    real(dp)  :: x, shift, stretch, shift2, stretch2
    integer   :: i

    width = dropoff_width
    if (cgs_units) then
      hallval = (mp_cgs * unit_velocity) / (ec_cgs * unit_length * unit_magneticfield)
      inertiaval = (mp_cgs * me_cgs * unit_velocity**2) / (ec_cgs * unit_length * unit_magneticfield)**2
    else
      hallval = (mp_si * unit_velocity) / (ec_si * unit_length * unit_magneticfield)
      inertiaval = (mp_si * me_si * unit_velocity**2) / (ec_si * unit_length * unit_magneticfield)**2
    end if

    sleft = grid_gauss(1) + 0.5d0 * width + dropoff_edge_dist
    sright = grid_gauss(gauss_gridpts) - dropoff_edge_dist - 0.5d0 * width

    if (hall_dropoff) then
      shift = hallval * tanh(-dpi) / (tanh(-dpi) - tanh(dpi))
      stretch = hallval / (tanh(dpi) - tanh(-dpi))

      do i = 1, gauss_gridpts
        x = grid_gauss(i)
        if (sleft - 0.5d0*width <= x .and. x <= sleft + 0.5d0*width) then
          hall_field % hallfactor(i) = shift + stretch * tanh(2.0d0 * dpi * (x - sleft) / width)
        else if (sleft + 0.5d0*width < x .and. x < sright - 0.5d0*width) then
          hall_field % hallfactor(i) = hallval
        else if (sright - 0.5d0*width <= x .and. x <= sright + 0.5d0*width) then
          hall_field % hallfactor(i) = shift + stretch * tanh(2.0d0 * dpi * (sright - x) / width)
        else
          hall_field % hallfactor(i) = 0.0d0
        end if
      end do
    else
      hall_field % hallfactor = hallval
    end if

    if (inertia_dropoff) then
      shift2 = inertiaval * tanh(-dpi) / (tanh(-dpi) - tanh(dpi))
      stretch2 = inertiaval / (tanh(dpi) - tanh(-dpi))

      do i = 1, gauss_gridpts
        x = grid_gauss(i)
        if (sleft - 0.5d0*width <= x .and. x <= sleft + 0.5d0*width) then
          hall_field % inertiafactor(i) = shift2 + stretch2 * tanh(2.0d0 * dpi * (x - sleft) / width)
        else if (sleft + 0.5d0*width < x .and. x < sright - 0.5d0*width) then
          hall_field % inertiafactor(i) = inertiaval
        else if (sright - 0.5d0*width <= x .and. x <= sright + 0.5d0*width) then
          hall_field % inertiafactor(i) = shift2 + stretch2 * tanh(2.0d0 * dpi * (sright - x) / width)
        else
          hall_field % inertiafactor(i) = 0.0d0
        end if
      end do
    else
      hall_field % inertiafactor = inertiaval
    end if
  end subroutine set_hallfactor

  !> Retrieves the A-matrix Hall contribution for a given quadblock corresponding
  !! to a particular point in the grid and a particular Gaussian weight.
  !! If <tt>hall_mhd = True</tt>, this routine is called <tt>n_gauss</tt> times
  !! for every grid interval.
  subroutine get_hallterm(gauss_idx, eps, d_eps_dr, curr_weight, quadblock_Hall, &
                          h_quadratic, dh_quadratic_dr, h_cubic, dh_cubic_dr)
    use mod_global_variables, only: dp_LIMIT, elec_pressure
    use mod_equilibrium_params, only: k2, k3
    use mod_equilibrium, only: rho_field, B_field, T_field
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
    !> array containing the derivatives of the 4 quadratic basis functions
    real(dp), intent(in)      :: dh_quadratic_dr(4)
    !> array containing the 4 cubic basis functions
    real(dp), intent(in)      :: h_cubic(4)
    !> array containing the derivatives of the 4 cubic basis functions
    real(dp), intent(in)      :: dh_cubic_dr(4)

    complex(dp), allocatable  :: factors(:)
    integer, allocatable      :: positions(:, :)

    real(dp)                  :: eps_inv, rho0_inv
    real(dp)                  :: rho0, B01, B02, B03, T0
    real(dp)                  :: drho0, drB02, dB02_r, dB03, dB02, dT0

    quadblock_Hall = (0.0d0, 0.0d0)
    eps_inv = 1.0d0 / eps

    ! Equilibrium quantities
    rho0 = rho_field % rho0(gauss_idx)
    drho0 = rho_field % d_rho0_dr(gauss_idx)
    rho0_inv = 1.0d0 / rho0

    T0 = T_field % T0(gauss_idx)
    dT0 = T_field % d_T0_dr(gauss_idx)

    ! To be uncommented when B01 is implemented
    B01 = 0.0d0 !B_field % B01(gauss_idx)
    B02 = B_field % B02(gauss_idx)
    dB02 = B_field % d_B02_dr(gauss_idx)
    B03 = B_field % B03(gauss_idx)
    dB03 = B_field % d_B03_dr(gauss_idx)

    ! Calculate derivatives eps*B02, B02/eps
    drB02 = d_eps_dr * B02 + eps * dB02
    dB02_r = eps_inv * dB02 - d_eps_dr * eps_inv**2 * B02

    ! Setup of matrix elements

    ! Quadratic * Quadratic
    call reset_factors(factors, 3)
    call reset_positions(positions, 3)
    ! A(6, 1)
    factors(1) = - eps_inv * rho0_inv**2 * (B02 * drB02 * eps_inv + B03 * dB03)
    positions(1, :) = [6, 1]
    ! A(6, 5)
    factors(2) = 0.0d0
    positions(2, :) = [6, 5]
    ! A(6, 6)
    factors(3) = - drho0 * rho0_inv**2 * (k2 * B03 * eps_inv - k3 * B02) &
                  + k3 * rho0_inv * (drB02 * eps_inv - eps * dB02_r)
    positions(3, :) = [6, 6]
    if (elec_pressure) then
      factors(1) = factors(1) - dT0 * eps_inv * rho0_inv
      factors(2) = factors(2) + drho0 * eps_inv * rho0_inv
    end if
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_quadratic, h_quadratic)

    ! d(Quadratic)/dr * Quadratic
    call reset_factors(factors, 3)
    call reset_positions(positions, 3)
    ! A(6, 1)
    factors(1) = 0.0d0
    positions(1, :) = [6, 1]
    ! A(6, 5)
    factors(2) = 0.0d0
    positions(2, :) = [6, 5]
    ! A(6, 6)
    factors(3) = rho0_inv * (k2 * B03 * eps_inv - k3 * B02)
    positions(3, :) = [6, 6]
    if (elec_pressure) then
      factors(1) = factors(1) - T0 * eps_inv * rho0_inv
      factors(2) = factors(2) - eps_inv
    end if
    call subblock(quadblock_Hall, factors, positions, curr_weight, dh_quadratic_dr, h_quadratic)

    ! Quadratic * Cubic
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(6, 7)
    factors(1) = - k3 * eps_inv * rho0_inv * (k2 * B02 * eps_inv + k3 * B03)
    positions(1, :) = [6, 7]
    ! A(6, 8)
    factors(2) = k2 * eps_inv * rho0_inv * (k2 * B02 * eps_inv + k3 * B03)
    positions(2, :) = [6, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_quadratic, h_cubic)

    ! Quadratic * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(6, 7)
    factors(1) = drho0 * rho0_inv**2 * B03 * eps_inv
    positions(1, :) = [6, 7]
    ! A(6, 8)
    factors(2) = rho0_inv * (eps * dB02_r - drB02 * eps_inv - drho0 * rho0_inv * B02)
    positions(2, :) = [6, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_quadratic, dh_cubic_dr)

    ! d(Quadratic)/dr * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(6, 7)
    factors(1) = - B03 * eps_inv * rho0_inv
    positions(1, :) = [6, 7]
    ! A(6, 8)
    factors(2) = B02 * rho0_inv
    positions(2, :) = [6, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, dh_quadratic_dr, dh_cubic_dr)

    ! Cubic * Quadratic
    call reset_factors(factors, 6)
    call reset_positions(positions, 6)
    ! A(7, 1)
    factors(1) = 0.0d0
    positions(1, :) = [7, 1]
    ! A(7, 5)
    factors(2) = 0.0d0
    positions(2, :) = [7, 5]
    ! A(7, 6)
    factors(3) = - B03 * rho0_inv * ((k2 * eps_inv)**2 + k3**2)
    positions(3, :) = [7, 6]
    ! A(8, 1)
    factors(4) = 0.0d0
    positions(4, :) = [8, 1]
    ! A(8, 5)
    factors(5) = 0.0d0
    positions(5, :) = [8, 5]
    ! A(8, 6)
    factors(6) = B02 * rho0_inv * ((k2 * eps_inv)**2 + k3**2)
    positions(6, :) = [8, 6]
    if (elec_pressure) then
      factors(1) = factors(1) + k2 * T0 * rho0_inv * eps_inv**2
      factors(2) = factors(2) + k2 * eps_inv**2
      factors(4) = factors(4) + k3 * T0 * rho0_inv * eps_inv
      factors(5) = factors(5) + k3 * eps_inv
    end if
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_cubic, h_quadratic)

    ! Cubic * Cubic
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(7, 7)
    factors(1) = k3 * drB02 * eps_inv**2 * rho0_inv
    positions(1, :) = [7, 7]
    ! A(7, 8)
    factors(2) = - k2 * drB02 * eps_inv**2 * rho0_inv
    positions(2, :) = [7, 8]
    ! A(8, 7)
    factors(3) = k3 * dB03 * eps_inv * rho0_inv
    positions(3, :) = [8, 7]
    ! A(8, 8)
    factors(4) = - k2 * dB03 * eps_inv * rho0_inv
    positions(4, :) = [8, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_cubic, h_cubic)

    ! Cubic * d(Cubic)/dr
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(7, 7)
    factors(1) = k2 * B03 * eps_inv**2 * rho0_inv
    positions(1, :) = [7, 7]
    ! A(7, 8)
    factors(2) = k3 * B03 * rho0_inv
    positions(2, :) = [7, 8]
    ! A(8, 7)
    factors(3) = - k2 * B02 * eps_inv**2 * rho0_inv
    positions(3, :) = [8, 7]
    ! A(8, 8)
    factors(4) = - k3 * B02 * rho0_inv
    positions(4, :) = [8, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_cubic, dh_cubic_dr)

    deallocate(factors)
    deallocate(positions)

    if (B01 > dp_LIMIT) then
      call get_hallterm_B01(gauss_idx, eps, d_eps_dr, curr_weight, &
                            quadblock_Hall, h_quadratic, h_cubic, dh_cubic_dr)
    end if
  end subroutine get_hallterm

  !> Retrieves the A-matrix Hall contribution for a given quadblock corresponding
  !! to a particular point in the grid and a particular Gaussian weight.
  !! If <tt>hall_mhd = True</tt>, this routine is called <tt>n_gauss</tt> times
  !! for every grid interval.
  subroutine get_hallterm_Bmat(gauss_idx, eps, d_eps_dr, curr_weight, quadblock_HallB, &
                          h_quadratic, h_cubic, dh_cubic_dr)
    use mod_equilibrium_params, only: k2, k3
    use mod_equilibrium, only: rho_field
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
    complex(dp), intent(out)  :: quadblock_HallB(dim_quadblock, dim_quadblock)
    !> array containing the 4 quadratic basis functions
    real(dp), intent(in)      :: h_quadratic(4)
    !> array containing the 4 cubic basis functions
    real(dp), intent(in)      :: h_cubic(4)
    !> array containing the derivatives of the 4 cubic basis functions
    real(dp), intent(in)      :: dh_cubic_dr(4)

    complex(dp), allocatable  :: factors(:)
    integer, allocatable      :: positions(:, :)

    real(dp)                  :: eps_inv, rho0_inv
    real(dp)                  :: rho0, drho0

    quadblock_HallB = (0.0d0, 0.0d0)
    eps_inv = 1.0d0 / eps

    ! Equilibrium quantities
    rho0 = rho_field % rho0(gauss_idx)
    drho0 = rho_field % d_rho0_dr(gauss_idx)
    rho0_inv = 1.0d0 / rho0

    ! Setup of matrix elements

    ! Quadratic * Quadratic
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(6, 6)
    factors(1) = rho0_inv * ((k2 * eps_inv)**2 + k3**2)
    positions(1, :) = [6, 6]
    call subblock(quadblock_HallB, factors, positions, curr_weight, h_quadratic, h_quadratic)

    ! Quadratic * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(6, 7)
    factors(1) = - k2 * rho0_inv * eps_inv**2
    positions(1, :) = [6, 7]
    ! A(6, 8)
    factors(2) = - k3 * rho0_inv
    positions(2, :) = [6, 8]
    call subblock(quadblock_HallB, factors, positions, curr_weight, h_quadratic, dh_cubic_dr)

    ! Cubic * Quadratic
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(7, 6)
    factors(1) = drho0 * k2 * eps_inv * rho0_inv**2
    positions(1, :) = [7, 6]
    ! A(8, 6)
    factors(2) = drho0 * k3 * rho0_inv**2 + d_eps_dr * k3 * eps_inv * rho0_inv
    positions(2, :) = [8, 6]
    call subblock(quadblock_HallB, factors, positions, curr_weight, h_cubic, h_quadratic)

    ! d(Cubic)/dr * Quadratic
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(7, 6)
    factors(1) = - k2 * eps_inv * rho0_inv
    positions(1, :) = [7, 6]
    ! A(8, 6)
    factors(2) = - k3 * rho0_inv
    positions(2, :) = [8, 6]
    call subblock(quadblock_HallB, factors, positions, curr_weight, dh_cubic_dr, h_quadratic)

    ! Cubic * Cubic
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(7, 7)
    factors(1) = k3**2 * eps_inv * rho0_inv
    positions(1, :) = [7, 7]
    ! A(7, 8)
    factors(2) = - k2 * k3 * eps_inv * rho0_inv
    positions(2, :) = [7, 8]
    ! A(8, 7)
    factors(3) = - k2 * k3 * eps_inv**2 * rho0_inv
    positions(3, :) = [8, 7]
    ! A(8, 8)
    factors(4) = (k2 * eps_inv)**2 * rho0_inv
    positions(4, :) = [8, 8]
    call subblock(quadblock_HallB, factors, positions, curr_weight, h_cubic, h_cubic)

    ! Cubic * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(7, 7)
    factors(1) = - drho0 * eps_inv * rho0_inv**2
    positions(1, :) = [7, 7]
    ! A(8, 8)
    factors(2) = - rho0_inv * (d_eps_dr * eps_inv + drho0 * rho0_inv)
    positions(2, :) = [8, 8]
    call subblock(quadblock_HallB, factors, positions, curr_weight, h_cubic, dh_cubic_dr)

    ! d(Cubic)/dr * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(7, 7)
    factors(1) = eps_inv * rho0_inv
    positions(1, :) = [7, 7]
    ! A(8, 8)
    factors(2) = rho0_inv
    positions(2, :) = [8, 8]
    call subblock(quadblock_HallB, factors, positions, curr_weight, dh_cubic_dr, dh_cubic_dr)

    deallocate(factors)
    deallocate(positions)
  end subroutine get_hallterm_Bmat

  !> Creates a quadblock for the A matrix containing the natural boundary
  !! conditions coming from the Hall term, depending on the supplied edge.
  !! @note  Currently, no boundary conditions are applied due to the use of a
  !!        dropoff profile for the Hall parameter, which is zero at the
  !!        edges. @endnote
  subroutine hall_boundaries(quadblock_Hall, kappa_perp_is_zero, edge)
    use mod_global_variables, only: gauss_gridpts, dim_subblock, elec_pressure
    use mod_logging, only: log_message
    use mod_equilibrium, only: rho_field, B_field, T_field
    use mod_grid, only: eps_grid
    use mod_equilibrium_params, only: k2, k3

    !> the quadblock corresponding to the left/right edge
    complex(dp), intent(out)      :: quadblock_Hall(dim_quadblock, dim_quadblock)
    !> the edge, either "l_edge" or "r_edge"
    character(len=6), intent(in)  :: edge
    !> boolean indicating if perpendicular thermal conduction is included or not
    logical, intent(in)           :: kappa_perp_is_zero

    complex(dp), allocatable  :: surface_terms(:)
    real(dp)                  :: eps, eps_inv, rho0_inv
    real(dp)                  :: rho0, B02, B03, T0
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
    rho0 = rho_field % rho0(idx)
    rho0_inv = 1.0d0 / rho0
    T0 = T_field % T0(idx)
    B02 = B_field % B02(idx)
    B03 = B_field % B03(idx)

    allocate(positions(4, 2))
    allocate(surface_terms(4))

    ! surface term for element (6, 5)
    surface_terms(1) = 0.0d0
    positions(1, :) = [6, 5]
    ! surface term for element (6, 6)
    surface_terms(2) = - rho0_inv * (k2 * B03 * eps_inv - k3 * B02)
    positions(2, :) = [6, 6]
    ! surface term for element (6, 7)
    surface_terms(3) = B03 * eps_inv * rho0_inv
    positions(3, :) = [6, 7]
    ! surface term for element (6, 8)
    surface_terms(4) = - B02 * rho0_inv
    positions(4, :) = [6, 8]

    ! T1 is zero at the wall if perpendicular thermal conduction is included
    if (elec_pressure .and. kappa_perp_is_zero) then
      surface_terms(2) = surface_terms(2) + eps_inv
    end if

    ! l_edge: add to bottom-right of 2x2 block, for top-left subblock only
    ! r_edge: add to bottom-right of 2x2 block, for bottom-right subblock only
    if (edge == 'l_edge') then
      positions = 2 * positions
    else if (edge == 'r_edge') then
      positions = 2 * positions + dim_subblock
    end if

    do i = 1, 2
      quadblock_Hall(positions(i, 1), positions(i, 2)) = quadblock_Hall(positions(i, 1), positions(i, 2)) + surface_terms(i)
    end do
    do i = 3, size(surface_terms)
      quadblock_Hall(positions(i, 1), positions(i, 2)) = quadblock_Hall(positions(i, 1), positions(i, 2)) + surface_terms(i)
    end do

    deallocate(positions)
    deallocate(surface_terms)
  end subroutine hall_boundaries

  !> Retrieves the B01 (constant) contributions to the A-matrix Hall
  !! contribution for a given quadblock corresponding
  !! to a particular point in the grid and a particular Gaussian weight.
  !! If <tt>hall_mhd = True</tt>, this routine is called <tt>n_gauss</tt> times
  !! for every grid interval.
  !! @note  Usage of a constant B01 component is currently not implemented. This
  !!        subroutine supports a future extension of the code. @endnote
  subroutine get_hallterm_B01(gauss_idx, eps, d_eps_dr, curr_weight, &
                              quadblock_Hall, h_quadratic, h_cubic, dh_cubic_dr)
    use mod_global_variables, only: ic
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
    complex(dp), intent(inout):: quadblock_Hall(dim_quadblock, dim_quadblock)
    !> array containing the 4 quadratic basis functions
    real(dp), intent(in)      :: h_quadratic(4)
    !> array containing the 4 cubic basis functions
    real(dp), intent(in)      :: h_cubic(4)
    !> array containing the derivatives of the 4 cubic basis functions
    real(dp), intent(in)      :: dh_cubic_dr(4)

    complex(dp), allocatable  :: factors(:)
    integer, allocatable      :: positions(:, :)

    real(dp)                  :: eps_inv, rho0_inv
    real(dp)                  :: rho0, B01, B02, B03
    real(dp)                  :: drho0, drB02, dB02_r, dB03, dB02

    eps_inv = 1.0d0 / eps

    ! Equilibrium quantities
    rho0 = rho_field % rho0(gauss_idx)
    drho0 = rho_field % d_rho0_dr(gauss_idx)
    rho0_inv = 1.0d0 / rho0

    ! To be uncommented when B01 is implemented
    B01 = 0.0d0 !B_field % B01(gauss_idx)
    B02 = B_field % B02(gauss_idx)
    dB02 = B_field % d_B02_dr(gauss_idx)
    B03 = B_field % B03(gauss_idx)
    dB03 = B_field % d_B03_dr(gauss_idx)

    ! Calculate derivatives eps*B02, B02/eps
    drB02 = d_eps_dr * B02 + eps * dB02
    dB02_r = eps_inv * dB02 - d_eps_dr * eps_inv**2 * B02

    ! Setup of matrix elements

    ! Cubic * Quadratic
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(7, 1)
    factors(1) = - ic * eps_inv**2 * rho0_inv**2 * B01 * drB02
    positions(1, :) = [7, 1]
    ! A(7, 6)
    factors(2) = ic * k3 * B01 * rho0_inv * (d_eps_dr * eps_inv + drho0 * rho0_inv)
    positions(2, :) = [7, 6]
    ! A(8, 1)
    factors(3) = - ic * eps_inv * rho0_inv**2 * B01 * dB03
    positions(3, :) = [8, 1]
    ! A(8, 6)
    factors(4) = - ic * drho0 * rho0_inv**2 * k2 * B01 * eps_inv
    positions(4, :) = [8, 6]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_cubic, h_quadratic)

    ! d(Cubic)/dr * Quadratic
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(7, 6)
    factors(1) = - ic * k3 * B01 * rho0_inv
    positions(1, :) = [7, 6]
    ! A(8, 6)
    factors(2) = ic * k2 * B01 * eps_inv * rho0_inv
    positions(2, :) = [8, 6]
    call subblock(quadblock_Hall, factors, positions, curr_weight, dh_cubic_dr, h_quadratic)

    ! Cubic * Cubic
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(7, 7)
    factors(1) = - ic * k2 * k3 * B01 * eps_inv**2 * rho0_inv
    positions(1, :) = [7, 7]
    ! A(7, 8)
    factors(2) = ic * k2**2 * B01 * eps_inv**2 * rho0_inv
    positions(2, :) = [7, 8]
    ! A(8, 7)
    factors(3) = - ic * k3**2 * B01 * eps_inv * rho0_inv
    positions(3, :) = [8, 7]
    ! A(8, 8)
    factors(4) = ic * k2 * k3 * B01 * eps_inv * rho0_inv
    positions(4, :) = [8, 8]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_cubic, h_cubic)

    ! Cubic * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(7, 8)
    factors(1) = ic * B01 * rho0_inv * (d_eps_dr * eps_inv + drho0 * rho0_inv)
    positions(1, :) = [7, 8]
    ! A(8, 7)
    factors(2) = ic * B01 * eps_inv * drho0 * rho0_inv**2
    positions(2, :) = [8, 7]
    call subblock(quadblock_Hall, factors, positions, curr_weight, h_cubic, dh_cubic_dr)

    ! d(Cubic)/dr * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(7, 8)
    factors(1) = - ic * B01 * rho0_inv
    positions(1, :) = [7, 8]
    ! A(8, 7)
    factors(2) = - ic * B01 * eps_inv * rho0_inv
    positions(2, :) = [8, 7]
    call subblock(quadblock_Hall, factors, positions, curr_weight, dh_cubic_dr, dh_cubic_dr)

    deallocate(factors)
    deallocate(positions)
  end subroutine get_hallterm_B01

end module mod_hallmhd
