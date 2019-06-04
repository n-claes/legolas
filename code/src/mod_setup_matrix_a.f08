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
    integer                   :: i, j, gauss_idx, k, l

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

        call get_A_elements(gauss_idx, eps, d_eps_dr, curr_weight, quadblock)

      end do  ! end do iteration Gaussian points

      ! Gridpoint i = 1, place quadblock in upper left corner of A
      if (i == 1) then
        do k = 1, dim_quadblock
          do l = 1, dim_quadblock
            matrix_A(k, l) = quadblock(k, l)
          end do
        end do
      end if

      ! TODO: fill matrix A with quadblocks

    end do    ! end do iteration gridpoints

  end subroutine construct_A

  subroutine get_A_elements(gauss_idx, eps, d_eps_dr, curr_weight, quadblock)
    use mod_setup_equilibrium
    use mod_make_subblock

    integer, intent(in)       :: gauss_idx
    real(dp), intent(in)      :: eps, d_eps_dr, curr_weight
    complex(dp), intent(out)  :: quadblock(dim_quadblock, dim_quadblock)

    real(dp)                  :: eps_inv
    real(dp)                  :: rho, v01, v02, v03, T, B01, B02, B03

    eps_inv = 1.0d0 / eps
    rho = rho0_eq(gauss_idx)
    v01 = v01_eq(gauss_idx)
    v02 = v02_eq(gauss_idx)
    v03 = v03_eq(gauss_idx)
    T   = T0_eq(gauss_idx)
    B01 = B01_eq(gauss_idx)
    B02 = B02_eq(gauss_idx)
    B03 = B03_eq(gauss_idx)

    ! Quadratic * Quadratic
    call reset_factors_A(19)
    call reset_positions(19)

    ! A(1, 1)
    factors_A(1) = eps_inv * (eps_inv * v02 * k2 + v03 * k3) * d_eps_dr
    positions(1, :) = [1, 1]
    ! A(1, 3)
    factors_A(2) = eps_inv * rho * k2
    positions(2, :) = [1, 3]
    ! A(1, 4)
    factors_A(3) = eps_inv * rho * k3

    ! TODO: others

    call subblock(quadblock, factors_A, positions, curr_weight, &
                  h_quadratic, h_quadratic)



    return

  end subroutine get_A_elements



  subroutine reset_factors_A(size_factors_A)
    integer, intent(in)       :: size_factors_A

    if (allocated(factors_A)) then
      deallocate(factors_A)
    end if

    allocate(factors_A(size_factors_A))
    factors_A = (0.0d0, 0.0d0)

  end subroutine reset_factors_A



  subroutine reset_positions(size_positions)
    integer, intent(in)       :: size_positions

    if (allocated(positions)) then
      deallocate(positions)
    end if

    allocate(positions(size_positions, 2))
    positions = 0

  end subroutine reset_positions


  subroutine matrix_A_clean()
    if (allocated(positions)) then
      deallocate(positions)
    end if
    if (allocated(factors_A)) then
      deallocate(factors_A)
    end if
  end subroutine matrix_A_clean

end module mod_setup_matrix_a
