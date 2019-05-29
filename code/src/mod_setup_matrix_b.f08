module mod_setup_matrix_b
  implicit none

  ! Sets up the B-matrix for the eigenvalue problem wBX = AX
  double precision                 :: h_cubic(4), h_quadratic(4)

contains

  subroutine construct_B(grid, grid_gauss, matrix_B)
    use mod_global_variables
    use mod_setup_equilibrium
    use mod_spline_functions

    double precision, intent(in)   :: grid(gridpts), grid_gauss(4*gridpts)
    double precision, intent(inout):: matrix_B(matrix_gridpts, matrix_gridpts)
    double precision               :: quadblock(dim_quadblock, dim_quadblock)
    double precision               :: r_lo, r_ce, r_hi, eps, d_eps_dr
    double precision               :: r
    complex                        :: factors_B(5)
    integer                        :: i, j, gauss_idx

    ! Initialize matrix to zero
    matrix_B = (0.0d0, 0.0d0)

    ! Iterate over gridpoints to calculate blocks
    do i = 2, gridpts-1
      ! Set block to zero
      quadblock = (0.0d0, 0.0d0)

      ! This integrates in the interval grid_gauss(i) to grid_gauss(i + 1)
      r_lo = grid(i - 1)
      r_ce = grid(i)
      r_hi = grid(i + 1)

      do j = 1, n_gauss
        ! Current grid index (from 4*n_gauss points)
        idx = (i - 1)*n_gauss + j
        r   = grid_gauss(gauss_idx)
        if (geometry .eq. "cartesian") then
          eps      = 1.0d0
          d_eps_dr = 0.0d0
        else if (geometry .eq. "cylindrical") then
          eps      = r
          d_eps_dr = 1.0d0
        endif

        ! TODO: FACTORS ARE PROBABLY WRONG
        call quadratic_factors(r, r_lo, r_ce, r_hi, h_quadratic)
        call cubic_factors(r, r_lo, r_ce, r_hi, h_cubic)

        call get_B_elements(gauss_idx, eps, d_eps_dr, quadblock, factors_B)

      end do
    end do

  end subroutine construct_B

  subroutine get_B_elements(gauss_idx, eps, d_eps_dr, quadblock, factors_B)
    use mod_global_variables
    use mod_setup_equilibrium
    use mod_make_subblock

    integer, intent(in)             :: gauss_idx
    double precision, intent(out)   :: quadblock(dim_quadblock, dim_quadblock)
    complex, intent(out)            :: factors_B(5)

    double precision                :: rho

    rho = rho_0(idx)

    ! Quadratic elements
    !> B(1,1)
    factors_B(1) = 1.0d0 / eps
    !> B(3,3)
    factors_B(2) = rho_0
    !> B(4,4)
    factors_B(3) = rho_0 / eps
    !> B(5,5)
    factors_B(4) = rho_0 / eps
    !> B(6,6)
    factors_B(5) = 1.0d0






end module mod_setup_matrix_b
