module mod_setup_matrix_b
  use mod_global_variables
  implicit none

  ! Sets up the B-matrix for the eigenvalue problem wBX = AX
  real(dp)                 :: h_cubic(4), h_quadratic(4)
  real(dp) :: Ht(4)
  ! Factors and positions are allocatable so they are dynamic, in case we add
  ! additional equations (self-gravity)
  complex(dp), allocatable :: factors_B(:)
  integer, allocatable     :: positions(:, :)

contains

  subroutine construct_B(grid, grid_gauss, matrix_B)
    use mod_setup_equilibrium
    use mod_spline_functions

    real(dp), intent(in)   :: grid(gridpts), grid_gauss(4*gridpts)
    real(dp), intent(inout):: matrix_B(matrix_gridpts, matrix_gridpts)
    complex(dp)            :: quadblock(dim_quadblock, dim_quadblock)
    real(dp)               :: r_lo, r_hi, eps, curr_weight
    real(dp)               :: r
    integer                :: i, j, gauss_idx, k, l

    ! Initialize matrix to zero
    matrix_B = 0.0d0

    ! Iterate over gridpoints to calculate blocks
    do i = 1, gridpts-1
      ! Set quadblock to zero (quadblock is complex, so tuple)
      quadblock = (0.0d0, 0.0d0)

      ! This integrates in the interval grid_gauss(i) to grid_gauss(i + 1)
      r_lo = grid(i)
      r_hi = grid(i + 1)

      do j = 1, n_gauss
        ! Current grid index (from 4*n_gauss points)
        gauss_idx = (i - 1)*n_gauss + j
        r   = grid_gauss(gauss_idx)
        if ((geometry .eq. "cartesian") .or. (geometry .eq. "Cartesian")) then
          eps      = 1.0d0
        else if (geometry .eq. "cylindrical") then
          eps      = r
        else
          write(*, *) "Geometry not defined correctly."
          write(*, *) "Currently set on:   ", geometry
          stop
        endif

        curr_weight = gaussian_weights(j)

        call quadratic_factors(r, r_lo, r_hi, h_quadratic)
        call cubic_factors(r, r_lo, r_hi, h_cubic)

        call get_B_elements(gauss_idx, eps, curr_weight, quadblock)

      end do    ! end do iteration gaussian points

      ! Gridpoint i = 1, place quadblock in upper left corner of B
      if (i == 1) then
        do k = 1, dim_quadblock
          do l = 1, dim_quadblock
            matrix_B(k, l) = quadblock(k, l)
          end do
        end do
      end if

    end do      ! end do iteration grid points

  end subroutine construct_B

  subroutine get_B_elements(gauss_idx, eps, curr_weight, quadblock)
    use mod_setup_equilibrium
    use mod_make_subblock

    integer, intent(in)          :: gauss_idx
    real(dp), intent(in)         :: eps, curr_weight
    complex(dp), intent(out)     :: quadblock(dim_quadblock, dim_quadblock)

    real(dp)                     :: rho

    rho = rho_0(gauss_idx)

    ! Quadratic elements
    call reset_factors_B(5)
    call reset_positions(5)

    ! B(1,1)
    factors_B(1) = 1.0d0 / eps
    positions(1, :) = [1, 1]
    ! B(3,3)
    factors_B(2) = rho
    positions(2, :) = [3, 3]
    ! B(4,4)
    factors_B(3) = rho / eps
    positions(3, :) = [4, 4]
    ! B(5,5)
    factors_B(4) = rho / eps
    positions(4, :) = [5, 5]
    ! B(6,6)
    factors_B(5) = 1.0d0
    positions(5, :) = [6, 6]

    call subblock(quadblock, factors_B, positions, curr_weight, h_quadratic, h_quadratic)

    ! Cubic elements
    call reset_factors_B(3)
    call reset_positions(3)

    ! B(2,2)
    factors_B(1) = rho / eps
    positions(1, :) = [2, 2]
    ! B(7,7)
    factors_B(2) = 1.0d0 / eps
    positions(2, :) = [7, 7]
    ! B(8,8)
    factors_B(3) = 1.0d0
    positions(3, :) = [8, 8]

    call subblock(quadblock, factors_B, positions, curr_weight, h_cubic, h_cubic)


  end subroutine get_B_elements



  subroutine reset_factors_B(size_factors_B)
    integer, intent(in)   :: size_factors_B

    if (allocated(factors_B)) then
      deallocate(factors_B)
    end if

    allocate(factors_B(size_factors_B))
    factors_B = (0.0d0, 0.0d0)

  end subroutine reset_factors_B



  subroutine reset_positions(size_positions)
    integer, intent(in)   :: size_positions

    if (allocated(positions)) then
      deallocate(positions)
    end if

    allocate(positions(size_positions, 2))
    positions = 0

  end subroutine reset_positions


  subroutine matrix_B_clean()
    if (allocated(positions)) then
      deallocate(positions)
    end if
    if (allocated(factors_B)) then
      deallocate(factors_B)
    end if
  end subroutine matrix_B_clean

end module mod_setup_matrix_b
