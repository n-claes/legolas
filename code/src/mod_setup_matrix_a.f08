module mod_setup_matrix_a
  use mod_global_variables
  implicit none

  ! Sets up the A-matrix for the eigenvalue problem wBX = AX
  real(dp)                 :: h_cubic(4), h_quadratic(4)
  real(dp)                 :: dh_cubic_dr(4), dh_quadratic_dr(4)
  ! Factors and positions are allocatable so they are dynamic, in case we add
  ! additional equations (self-gravity)
  complex(dp), allocatable :: factors_A(:)
  integer, allocatable     :: positions(:, :)

contains

  subroutine construct_A(grid, grid_gauss, matrix_A)
    use mod_setup_equilibrium
    use mod_spline_functions

    real(dp), intent(in)      :: grid(gridpts), grid_gauss(4*gridpts)
    complex(dp), intent(out)  :: matrix_A(matrix_gridpts, matrix_gridpts)
    complex(dp)               :: quadblock(dim_quadblock, dim_quadblock)
    real(dp)                  :: r_lo, r_hi, eps, d_eps_dr, curr_weight
    real(dp)                  :: r
    integer                   :: i, j, gauss_idx

    ! Initialise matrix to zero (A is complex)
    matrix_A = (0.0d0, 0.0d0)

    do i = 2, gridpts-1
      ! Set quadblock to zero
      quadblock = (0.0d0, 0.0d0)

      r_lo = grid(i)
      r_hi = grid(i + 1)

      do j = 1, n_gauss
        ! Current grid index (from 4*n_gauss points)
        gauss_idx = (i - 1)*n_gauss + j
        r = grid_gauss(gauss_idx)
        if ((geometry .eq. "cartesian") .or. (geometry .eq. "Cartesian")) then
          eps      = 1.0d0
          d_eps_dr = 0.0d0
        else if (geometry .eq. "cylindrical") then
          eps      = r
          d_eps_dr = 1.0d0
        else
          write(*,*) "Geometry not defined correctly."
          write(*,*) "Currently set on:    ", geometry
          stop
        end if

        curr_weight = gaussian_weights(j)

        call quadratic_factors(r, r_lo, r_hi, h_quadratic)
        call quadratic_factors_deriv(r, r_lo, r_hi, dh_quadratic_dr)
        call cubic_factors(r, r_lo, r_hi, h_cubic)
        call cubic_factors_deriv(r, r_lo, r_hi, dh_cubic_dr)

      end do
    end do

  end subroutine construct_A

  subroutine get_A_elements(gauss_idx, eps, d_eps_dr, curr_weight, quadblock)
    use mod_setup_equilibrium
    use mod_make_subblock

    integer, intent(in)       :: gauss_idx
    real(dp), intent(in)      :: eps, d_eps_dr, curr_weight
    complex(dp), intent(out)  :: quadblock(dim_quadblock, dim_quadblock)

    ! TODO factors

    return

  end subroutine get_A_elements




end module mod_setup_matrix_a
