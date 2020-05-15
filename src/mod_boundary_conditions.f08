!
! MODULE: mod_boundary_conditions
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module to calculate the boundary conditions for the eigenvalue problem.
!
module mod_boundary_conditions
  use mod_global_variables, only: dp, matrix_gridpts, dim_quadblock, dim_subblock
  implicit none

  private

  public :: apply_boundary_conditions

contains

  subroutine apply_boundary_conditions(matrix_A, matrix_B)
    complex(dp), intent(inout)  :: matrix_A(matrix_gridpts, matrix_gridpts)
    real(dp), intent(inout)     :: matrix_B(matrix_gridpts, matrix_gridpts)
    complex(dp)                 :: quadblock(dim_quadblock, dim_quadblock)
    integer                     :: idx_end_left, idx_start_right

    ! end of index first-gridpoint quadblock
    idx_end_left = dim_quadblock
    ! start of index last-gridpoint quadblock
    idx_start_right = matrix_gridpts - dim_quadblock + 1

    ! matrix B left-edge quadblock
    quadblock = matrix_B(1:idx_end_left, 1:idx_end_left)
    call essential_boundaries(quadblock, edge='l_edge', matrix='B')
    matrix_B(1:idx_end_left, 1:idx_end_left) = real(quadblock)
    ! matrix B right-edge quadblock
    quadblock = matrix_B(idx_start_right:matrix_gridpts, idx_start_right:matrix_gridpts)
    call essential_boundaries(quadblock, edge='r_edge', matrix='B')
    matrix_B(idx_start_right:matrix_gridpts, idx_start_right:matrix_gridpts) = real(quadblock)

    ! matrix A left-edge quadblock
    quadblock = matrix_A(1:idx_end_left, 1:idx_end_left)
    call natural_boundaries(quadblock, edge='l_edge')
    call essential_boundaries(quadblock, edge='l_edge', matrix='A')
    matrix_A(1:idx_end_left, 1:idx_end_left) = quadblock
    ! matrix A right-edge quadblock
    quadblock = matrix_A(idx_start_right:matrix_gridpts, idx_start_right:matrix_gridpts)
    call natural_boundaries(quadblock, edge='r_edge')
    call essential_boundaries(quadblock, edge='r_edge', matrix='A')
    matrix_A(idx_start_right:matrix_gridpts, idx_start_right:matrix_gridpts) = quadblock
  end subroutine apply_boundary_conditions

  subroutine essential_boundaries(quadblock, edge, matrix)
    use mod_global_variables, only: boundary_type, use_fixed_tc_perp, &
                                    fixed_tc_perp_value, dp_LIMIT
    use mod_logging, only: log_message

    complex(dp), intent(inout)    :: quadblock(dim_quadblock, dim_quadblock)
    character(len=6), intent(in)  :: edge
    character, intent(in)         :: matrix

    complex(dp)                   :: diagonal_factor
    integer                       :: i, j, qua_zeroes(5), wall_idx_left(4), wall_idx_right(4)
    logical                       :: use_Tbound

    if (matrix == 'B') then
      diagonal_factor = (1.0d0, 0.0d0)
    else if (matrix == 'A') then
      diagonal_factor = (1.0d20, 0.0d0)
    else
      call log_message("essential boundaries: invalid matrix argument", level='error')
    end if

    ! determine if additional boundary conditions on temperature must be applied.
    ! this is only done if there is perpendicular thermal conduction
    use_Tbound = .false.
    if (use_fixed_tc_perp .and. abs(fixed_tc_perp_value) > dp_LIMIT) then
      use_Tbound = .true.
    end if

    ! Always: the contribution from the 0 basis function automatically
    ! zeroes out the odd rows/columns for the quadratic variables on the left edge
    ! so we handle those indices explicitly
    qua_zeroes = [1, 5, 7, 9, 11]
    if (edge == 'l_edge') then
      do i = 1, size(qua_zeroes)
        j = qua_zeroes(i)
        quadblock(j, j) = diagonal_factor
      end do
    end if

    ! Wall/regularity conditions: handling of v1, a2 and a3 (and T if conduction).
    ! v1, a2 and a3 are cubic elements, so omit non-zero basis functions (odd rows/columns)
    ! T is a quadratic element, so omit even row/columns
    wall_idx_left = [3, 13, 15, 9]
    wall_idx_right = [19, 29, 31, 26]

    select case(boundary_type)
    case('wall')
      if (edge == 'l_edge') then
        ! left regularity/wall conditions
        do i = 1, size(wall_idx_left)
          j = wall_idx_left(i)
          if (j == 9 .and. .not. use_Tbound) then
            cycle
          end if
          quadblock(j, :) = (0.0d0, 0.0d0)
          quadblock(:, j) = (0.0d0, 0.0d0)
          quadblock(j, j) = diagonal_factor
        end do
      else if (edge == 'r_edge') then
        do i = 1, size(wall_idx_right)
          j = wall_idx_right(i)
          if ((j == 26) .and. .not. use_Tbound) then
            cycle
          end if
          quadblock(j, :) = (0.0d0, 0.0d0)
          quadblock(:, j) = (0.0d0, 0.0d0)
          quadblock(j, j) = diagonal_factor
        end do
      else
        call log_message("essential boundaries: invalid edge argument", level='error')
      end if

    case default
      call log_message( "essential boundaries: invalid boundary_type", level='error')
    end select

  end subroutine essential_boundaries

  subroutine natural_boundaries(quadblock, edge)
    use mod_global_variables, only: ic, gamma_1, gauss_gridpts, gridpts
    use mod_equilibrium_params, only: k2, k3
    use mod_spline_functions, only: quadratic_factors, quadratic_factors_deriv, &
                                    cubic_factors, cubic_factors_deriv
    use mod_grid, only: grid, eps_grid, d_eps_grid_dr
    use mod_equilibrium, only: rho_field, T_field, B_field, kappa_field, eta_field
    use mod_make_subblock, only: reset_factors, reset_positions
    use mod_logging, only: log_message

    complex(dp), intent(inout)    :: quadblock(dim_quadblock, dim_quadblock)
    character(len=6), intent(in)  :: edge

    complex(dp), allocatable      :: factors(:)
    integer, allocatable          :: positions(:, :)
    real(dp)  :: h_quadratic(4), h_cubic(4), dh_quadratic_dr(4), dh_cubic_dr(4)
    real(dp)  :: r, r_lo, r_hi, eps_inv, eps, d_eps_dr
    real(dp)  :: rho0, T0, B02, B03, dB02, dB03, drB02, eta, dT0, &
                 tc_perp, dtc_perp_dT, dtc_perp_dB2, dtc_perp_drho
    integer   :: idx

    if (edge == 'l_edge') then
      r_lo = grid(1)
      r_hi = grid(2)
      r = r_lo
      idx = 1
    else if (edge == 'r_edge') then
      r_lo = grid(gridpts - 1)
      r_hi = grid(gridpts)
      r = r_hi
      idx = gauss_gridpts
    else
      call log_message("natural boundaries: invalid edge argument", level='error')
    end if

    eps = eps_grid(idx)
    d_eps_dr = d_eps_grid_dr(idx)
    eps_inv = 1.0d0 / eps

    ! Equilibrium quantities at the boundary
    rho0 = rho_field % rho0(idx)
    T0 = T_field % T0(idx)
    dT0 = T_field % d_T0_dr(idx)
    B02 = B_field % B02(idx)
    B03 = B_field % B03(idx)
    dB02 = B_field % d_B02_dr(idx)
    dB03 = B_field % d_B03_dr(idx)
    eta = eta_field % eta(idx)
    tc_perp = kappa_field % kappa_perp(idx)
    dtc_perp_dT = kappa_field % d_kappa_perp_dT(idx)
    dtc_perp_drho = kappa_field % d_kappa_perp_drho(idx)
    dtc_perp_dB2 = kappa_field % d_kappa_perp_dB2(idx)
    drB02 = d_eps_dr * B02 + eps * dB02
    ! Spline functions at the boundaries
    call quadratic_factors(r, r_lo, r_hi, h_quadratic)
    call quadratic_factors_deriv(r, r_lo, r_hi, dh_quadratic_dr)
    call cubic_factors(r, r_lo, r_hi, h_cubic)
    call cubic_factors_deriv(r, r_lo, r_hi, dh_cubic_dr)

    ! Quadratic * Quadratic
    call reset_factors(factors, 3)
    call reset_positions(positions, 3)
    ! A(5, 1)
    factors(1) = ic * gamma_1 * eps_inv * dT0 * dtc_perp_drho
    positions(1, :) = [5, 1]
    ! A(5, 5)
    factors(2) = ic * gamma_1 * eps_inv * (dT0 * dtc_perp_dT - d_eps_dr * eps_inv * tc_perp)
    positions(2, :) = [5, 5]
    ! A(5, 6)
    factors(3) = 2.0d0 * ic * gamma_1 * eps_inv * ( dT0 * (eps * B02 * k3 - B03 * k2) * dtc_perp_dB2 &
                                                    - eta * dB03 * k2 + eta * k3 * drB02)
    positions(3, :) = [5, 6]
    call add_factors_quadblock(quadblock, factors, positions, h_quadratic, h_quadratic, edge)

    ! Quadratic *  d(Quadratic)/dr
    call reset_factors(factors, 1)
    call reset_positions(positions, 1)
    ! A(5, 5)
    factors(1) = ic * gamma_1 * eps_inv * tc_perp
    positions(1, :) = [5, 5]
    call add_factors_quadblock(quadblock, factors, positions, h_quadratic, dh_quadratic_dr, edge)

    ! Quadratic * d(Cubic)/dr
    call reset_factors(factors, 2)
    call reset_positions(positions, 2)
    ! A(5, 7)
    factors(1) = 2.0d0 * ic * gamma_1 * eps_inv * (dT0 * B03 * dtc_perp_dB2 + eta * dB03)
    positions(1, :) = [5, 7]
    ! A(5, 8)
    factors(2) = -2.0d0 * ic * gamma_1 * (dT0 * B02 * dtc_perp_dB2 + eta * eps_inv * drB02)
    positions(2, :) = [5, 8]
    call add_factors_quadblock(quadblock, factors, positions, h_quadratic, dh_cubic_dr, edge)

    ! Cubic * Quadratic
    call reset_factors(factors, 5)
    call reset_positions(positions, 5)
    ! A(2, 1)
    factors(1) = eps_inv * T0
    positions(1, :) = [2, 1]
    ! A(2, 5)
    factors(2) = eps_inv * rho0
    positions(2, :) = [2, 5]
    ! A(2, 6)
    factors(3) = B02 * k3 - eps_inv * B03 * k2
    positions(3, :) = [2, 6]
    ! A(7, 6)
    factors(4) = -ic * eta * eps_inv * k2
    positions(4, :) = [7, 6]
    ! A(8, 6)
    factors(5) = -ic * eta * k3
    positions(5, :) = [8, 6]
    call add_factors_quadblock(quadblock, factors, positions, h_cubic, h_quadratic, edge)

    ! Cubic * d(Cubic)/dr
    call reset_factors(factors, 4)
    call reset_positions(positions, 4)
    ! A(2, 7)
    factors(1) = eps_inv * B03
    positions(1, :) = [2, 7]
    ! A(2, 8)
    factors(2) = -B02
    positions(2, :) = [2, 8]
    ! A(7, 7)
    factors(3) = ic * eta * eps_inv
    positions(3, :) = [7, 7]
    ! A(8, 8)
    factors(4) = ic * eta
    positions(4, :) = [8, 8]
    call add_factors_quadblock(quadblock, factors, positions, h_cubic, dh_cubic_dr, edge)

    deallocate(factors)
    deallocate(positions)
  end subroutine natural_boundaries

  !> Adds the factors to the quadblock depending on their positions and the
  !! value of 'edge'.
  !! @param[in, out]  quadblock   The 32x32 quadblock. Out: factors added to
  !!                              their respective blocks
  !! @param[in] factors   The boundary conditions for each element
  !! @param[in] positions The positions of each boundary condition
  !! @param[in] spline1   The left spline to apply, evaluated at the edge
  !! @param[in] spline2   The right spline to apply, evaluated at the edge
  !! @param[in] edge  'l_edge' for left boundary, 'r_edge' for right boundary
  !!                  - left edge : do not fill the bottom-right subblock
  !!                  - right edge: do not fill the top-left subblock
  subroutine add_factors_quadblock(quadblock, factors, positions, spline1, spline2, edge)
    complex(dp), intent(inout)  :: quadblock(dim_quadblock, dim_quadblock)
    complex(dp), intent(in)     :: factors(:)
    integer, intent(in)         :: positions(:, :)
    real(dp), intent(in)        :: spline1(4), spline2(4)
    character(6), intent(in)    :: edge

    integer                     :: i, len_factors
    integer                     :: curr_position(2), idx(2)

     idx(:) = 0
     len_factors = size(factors)

     do i = 1, len_factors
       curr_position = positions(i, :)

       ! Top-left corner subblock, only for left boundary conditions
       if (edge == 'l_edge') then
         idx(:) = 2*curr_position(:)

         quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                         spline1(2) * factors(i) * spline2(2)
         quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                         spline1(2) * factors(i) * spline2(4)
         quadblock(idx(1),   idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                         spline1(4) * factors(i) * spline2(2)
         quadblock(idx(1),   idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                         spline1(4) * factors(i) * spline2(4)
       end if

       ! Subblock top-right corner, filled both for left and right bounds
       idx(:) = [2*curr_position(1), 2*curr_position(2) + dim_subblock]
       quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                       spline1(2) * factors(i) * spline2(1)
       quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                       spline1(2) * factors(i) * spline2(3)
       quadblock(idx(1)  , idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                       spline1(4) * factors(i) * spline2(1)
       quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                       spline1(4) * factors(i) * spline2(3)

       ! Subblock bottom-left corner, filled both for left and right bounds
       idx(:) = [2*curr_position(1) + dim_subblock, 2*curr_position(2)]
       quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                       spline1(1) * factors(i) * spline2(2)
       quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                       spline1(1) * factors(i) * spline2(4)
       quadblock(idx(1)  , idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                       spline1(3) * factors(i) * spline2(2)
       quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                       spline1(3) * factors(i) * spline2(4)

       if (edge == 'r_edge') then
         ! Subblock bottom-right corner, only for right boundary conditions
         idx(:) = 2*curr_position(:) + dim_subblock
         quadblock(idx(1)-1, idx(2)-1) = quadblock(idx(1)-1, idx(2)-1) + &
                                         spline1(1) * factors(i) * spline2(1)
         quadblock(idx(1)-1, idx(2)  ) = quadblock(idx(1)-1, idx(2)  ) + &
                                         spline1(1) * factors(i) * spline2(3)
         quadblock(idx(1)  , idx(2)-1) = quadblock(idx(1)  , idx(2)-1) + &
                                         spline1(3) * factors(i) * spline2(1)
         quadblock(idx(1)  , idx(2)  ) = quadblock(idx(1)  , idx(2)  ) + &
                                         spline1(3) * factors(i) * spline2(3)
       end if
    end do
  end subroutine add_factors_quadblock

end module mod_boundary_conditions
