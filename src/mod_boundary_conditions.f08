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
  use mod_global_variables, only: dp, dim_quadblock
  implicit none

  public :: boundaries_B_left_edge
  public :: boundaries_B_right_edge
  public :: boundaries_A_left_edge
  public :: boundaries_A_right_edge

contains

  !> Boundary conditions matrix B, left edge (first iteration over quadblock).
  !! B-matrix has no natural boundary conditions, only essential ones.
  !! Dirichlet (fixed) boundary conditions are implemented at the left edge.
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: boundary conditions applied.
  subroutine boundaries_B_left_edge(quadblock)
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)

    call fixed_boundaries(quadblock, "l_edge", "B")

  end subroutine boundaries_B_left_edge


  !> Boundary conditions matrix B, right edge (final iteration over quadblock).
  !! B-matrix has no natural boundary conditions, only essential ones.
  !! Dirichlet (fixed) boundary conditions are implemented at the right edge.
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: boundary conditions applied.
  subroutine boundaries_B_right_edge(quadblock)
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)

    call fixed_boundaries(quadblock, "r_edge", "B")

  end subroutine boundaries_B_right_edge


  !> Boundary conditions matrix A, left edge. A-matrix has both natural
  !! (from partial integration) and essential (fixed) boundary conditions.
  !! @param[in] eps   The value for epsilon: 1 for Cartesian, r for cylindrical
  !! @param[in] d_eps_dr  Epsilon derivative: 0 for Cartesian, 1 for cylindrical
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: boundary conditions applied.
  subroutine boundaries_A_left_edge(quadblock, eps, d_eps_dr)
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)
    real(dp), intent(in)       :: eps, d_eps_dr

    call natural_boundaries(eps, d_eps_dr, quadblock, "l_edge")
    call fixed_boundaries(quadblock, "l_edge", "A")

  end subroutine boundaries_A_left_edge


  !> Boundary conditions matrix A, right edge. A-matrix has both natural
  !! (from partial integration) and essential (fixed) boundary conditions.
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: boundary conditions applied.
  subroutine boundaries_A_right_edge(quadblock, eps, d_eps_dr)
    complex(dp), intent(inout) :: quadblock(dim_quadblock, dim_quadblock)
    real(dp), intent(in)       :: eps, d_eps_dr

    call natural_boundaries(eps, d_eps_dr, quadblock, "r_edge")
    call fixed_boundaries(quadblock, "r_edge", "A")

  end subroutine boundaries_A_right_edge


  !> Boundaries originating from the fixed conditions. When edge = l_edge,
  !! set regularity conditions for the left gridpoint. If edge = r_edge,
  !! set them for the right gridpoint.
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: fixed boundary conditions applied.
  !! @param[in] edge    'r_edge' for right boundary, 'l_edge' for left boundary
  !! @param[in] matrix  'A' for A-matrix boundaries, 'B' for B-matrix boundaries
  subroutine fixed_boundaries(quadblock, edge, matrix)
    use mod_global_variables, only: dim_subblock, boundary_type, thermal_conduction

    complex(dp), intent(inout)  :: quadblock(dim_quadblock, dim_quadblock)
    character(6), intent(in)    :: edge
    character, intent(in)       :: matrix
    integer                     :: i
    complex(dp)                 :: unity

    !! LEDA sets diagonal elements to ZBIG instead of unity. This is
    !! done to clearly indicate which additional eigenvalues are introduced
    !! by forcing the boundary conditions. For now, these are set to zero
    !! (so we introduce eigenvalues at 0).
    if (matrix == "B") then
      unity = (1.0d0, 0.0d0)
    else if (matrix == "A") then
      unity = (0.0d0, 0.0d0)
    else
      write(*, *) "Wrong matrix argument passed to boundaries."
      write(*, *) "Currently set on: ", matrix
      stop
    end if

    !! === LEFT EDGE ===
    !! For the two subbblocks on the left edge of the quadblock, set the
    !! odd columns to zero. For the two subblocks on the top edge of the
    !! quadblock, set the odd rows to zero (such that the right-bottom
    !! subblock is not modified as this one is not at the boundary).
    !! The first diagonal elements of the top-left quadblock are set to unity.
    if (edge == "l_edge") then
      do i = 1, dim_subblock, 2
        ! odd rows to zero
        quadblock(i, :) = (0.0d0, 0.0d0)
        ! odd columns to zero
        quadblock(:, i) = (0.0d0, 0.0d0)
        ! Unity on first diagonal elements for B, zero for A
        quadblock(i, i) = unity
      end do
      !! Additional requirements on T1 if thermal conduction is included
      if (thermal_conduction) then
        !! second row + column of T1 set to zero
        quadblock(10, :) = (0.0d0, 0.0d0)
        quadblock(:, 10) = (0.0d0, 0.0d0)
        !! Unity on second diagonal element for B, zero for A
        quadblock(10, 10) = unity
      end if

    !! === RIGHT EDGE ===
    !! For the rigid wall on the right edge, we require v1 = a2 = a3 = 0.
    else if (edge == "r_edge") then

      select case(boundary_type)
        case("wall")
          !! Rows and columns for v1 (index 19 in bottom-right subblock)
          quadblock(19, :) = (0.0d0, 0.0d0)
          quadblock(:, 19) = (0.0d0, 0.0d0)
          !! Insert unity at first diagonal element for B, zero for A
          quadblock(19, 19) = unity

          !! Rows and columns for a2 (index 29 in bottom-right subblock)
          quadblock(29, :) = (0.0d0, 0.0d0)
          quadblock(:, 29) = (0.0d0, 0.0d0)
          !! Insert unity at first diagonal element for B, zero for A
          quadblock(29, 29) = unity

          !! Rows and columns for a3 (index 31 in bottom-right subblock)
          quadblock(31, :) = (0.0d0, 0.0d0)
          quadblock(:, 31) = (0.0d0, 0.0d0)
          !! Insert unity at first diagonal element for B, zero for A
          quadblock(31, 31) = unity

          !! Additional requiremens on T1 if thermal conduction is included
          !! Rows and columns for T1 (2nd) (index 26 in bottom-right subblock)
          if (thermal_conduction) then
            quadblock(26, :) = (0.0d0, 0.0d0)
            quadblock(:, 26) = (0.0d0, 0.0d0)
            !! Insert unity at second diagonal element for B, zero for A
            quadblock(26, 26) = unity
          end if

        case default
          write(*, *) "No proper boundary condition selected."
          stop
      end select

    else    ! if-case edge
      write(*, *) "Wrong edge passed when calling boundaries."
      write(*, *) "Currently set on: ", edge
      stop
    end if

  end subroutine fixed_boundaries


  !> Boundary conditions originating from the partially integrated terms of the
  !! A-matrix.
  !! @param[in] eps   The value for epsilon: 1 for Cartesian, r for cylindrical
  !! @param[in] d_eps_dr  Epsilon derivative: 0 for Cartesian, 1 for cylindrical
  !! @param[in, out] quadblock    The 32x32 quadblock, containing 4 subblocks.
  !!                              Out: natural boundary conditions applied.
  !! @param[in] edge  'r_edge' for right boundary, 'l_edge' for left boundary
  subroutine natural_boundaries(eps, d_eps_dr, quadblock, edge)
    use mod_global_variables, only: ic, gamma_1, gauss_gridpts, gridpts
    use mod_equilibrium_params, only: k2, k3
    use mod_spline_functions, only: quadratic_factors, quadratic_factors_deriv, cubic_factors, cubic_factors_deriv
    use mod_grid, only: grid
    use mod_equilibrium, only: rho_field, T_field, B_field, kappa_field, eta_field
    use mod_make_subblock, only: subblock, reset_factors, reset_positions

    real(dp), intent(in)        :: eps, d_eps_dr
    complex(dp), intent(inout)  :: quadblock(dim_quadblock, dim_quadblock)
    character(6), intent(in)    :: edge

    complex(dp), allocatable    :: factors(:)
    integer, allocatable        :: positions(:, :)
    real(dp)                    :: h_quadratic(4), h_cubic(4)
    real(dp)                    :: dh_quadratic_dr(4), dh_cubic_dr(4)
    real(dp)                    :: r, r_lo, r_hi, eps_inv
    integer                     :: idx

    real(dp)                    :: rho0, T0, B02, B03
    real(dp)                    :: tc_perp, eta
    real(dp)                    :: drB02, dT0, dB02, dB03
    real(dp)                    :: dtc_perp_dT, dtc_perp_drho, dtc_perp_dB2

    !! \note: for the boundaries we require the GRID, not grid_gauss.
    if (edge == 'l_edge') then
     r_lo = grid(1)
     r_hi = grid(2)
     r    = r_lo
     idx  = 1
    else if (edge == "r_edge") then
     r_lo = grid(gridpts-1)
     r_hi = grid(gridpts)
     r    = r_hi
     idx  = gauss_gridpts
    else
     write(*, *) "Wrong edge passed to natural boundaries routine"
     write(*, *) "Currently set on: ", edge
     stop
    end if

    eps_inv = 1.0d0 / eps

    !! Equilibrium quantities at the boundary
    rho0    = rho_field % rho0(idx)
    T0      = T_field % T0(idx)
    dT0     = T_field % d_T0_dr(idx)
    B02     = B_field % B02(idx)
    B03     = B_field % B03(idx)
    dB02    = B_field % d_B02_dr(idx)
    dB03    = B_field % d_B03_dr(idx)
    eta     = eta_field % eta(idx)
    tc_perp = kappa_field % kappa_perp(idx)
    dtc_perp_dT   = kappa_field % d_kappa_perp_dT(idx)
    dtc_perp_drho = kappa_field % d_kappa_perp_drho(idx)
    dtc_perp_dB2  = kappa_field % d_kappa_perp_dB2(idx)

    drB02 = d_eps_dr * B02 + eps * dB02

    !! Spline functions for the boundaries. Interval is [grid(1), grid(2)]
    !! for the left edge, [grid(N-1), grid(N)] for the right edge.
    !! r is equal to grid(1) for left, grid(N) for right
    call quadratic_factors(r, r_lo, r_hi, h_quadratic)
    call quadratic_factors_deriv(r, r_lo, r_hi, dh_quadratic_dr)
    call cubic_factors(r, r_lo, r_hi, h_cubic)
    call cubic_factors_deriv(r, r_lo, r_hi, dh_cubic_dr)

    !! Calculating boundary terms

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
    use mod_global_variables, only: dim_subblock

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
