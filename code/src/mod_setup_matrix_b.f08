!
! MODULE: mod_setup_matrix_b
!
!> @author
!> Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Module to create the B-matrix in the eigenvalue problem AX = wBX.
!
module mod_setup_matrix_b
  use mod_global_variables
  implicit none

  private

  !> @Note: Factors and positions are allocatable so they are dynamic,
  !!        in case we add additional equations (eg. self-gravity)

  !> Array containing the 4 quadratic basic functions
  real(dp)                 :: h_quadratic(4)
  !> Array containing the 4 cubic basic functions
  real(dp)                 :: h_cubic(4)

  public construct_B

contains

  !> Constructs the B-matrix.
  !! @param[in, out] matrix_B   The B-matrix in AX = wBX
  subroutine construct_B(matrix_B)
    use mod_grid
    use mod_equilibrium
    use mod_spline_functions
    use mod_boundary_conditions

    real(dp), intent(inout):: matrix_B(matrix_gridpts, matrix_gridpts)
    complex(dp)            :: quadblock(dim_quadblock, dim_quadblock)
    real(dp)               :: r_lo, r_hi, eps, curr_weight
    real(dp)               :: r
    integer                :: i, j, gauss_idx, k, l
    integer                :: quadblock_idx, idx1, idx2

    ! Initialize matrix to zero (B is real)
    matrix_B = 0.0d0

    ! Initialise quadblock index to zero (used to shift the block)
    quadblock_idx = 0

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

      !! Solving in [0, 1], so scale with difference
      quadblock = quadblock * (r_hi - r_lo)

      !! Apply boundary conditions on edges
      if (i == 1) then
        call boundaries_B_left_edge(quadblock)
      end if
      if (i == gridpts - 1) then
        call boundaries_B_right_edge(quadblock)
      end if

      ! Fill B-matrix with quadblock entries
      do k = 1, dim_quadblock
        do l = 1, dim_quadblock
          idx1 = k + quadblock_idx
          idx2 = l + quadblock_idx
          matrix_B(idx1, idx2) = matrix_B(idx1, idx2) + real(quadblock(k, l))
        end do
      end do
      !! This ensures that the quadblock is shifted on the main diagonal using
      !! idx1 and idx2. Dimension dim_subblock (= 16) is added instead of
      !! dim_quadblock (= 32), as the bottom-right part of the quadblock
      !! overlaps with the top-left part when shifting along the main diagonal.
      quadblock_idx = quadblock_idx + dim_subblock

    end do      ! end do iteration grid points

  end subroutine construct_B

  !> Calculates the different integral elements for the B-matrix.
  !! @param[in] gauss_idx   The current index in the Gaussian grid (r-position)
  !! @param[in] eps         The value for the scale factor epsilon
  !! @param[in] curr_weight The current weight for the Gaussian quadrature
  !! @param[in, out] quadblock  The quadblock, used to calculate the B-matrix.
  !!                            This block is shifted on the main diagonal
  subroutine get_B_elements(gauss_idx, eps, curr_weight, quadblock)
    use mod_equilibrium
    use mod_make_subblock

    integer, intent(in)          :: gauss_idx
    real(dp), intent(in)         :: eps, curr_weight
    complex(dp), intent(inout)   :: quadblock(dim_quadblock, dim_quadblock)

    !> Array containing the integral expressions for the B-matrix
    complex(dp), allocatable     :: factors(:)
    !> Array containing the position of each integral, eg. [1, 3] for B(1, 3)
    integer, allocatable         :: positions(:, :)
    real(dp)                     :: rho

    rho = rho0_eq(gauss_idx)

    ! Quadratic * Quadratic
    call reset_factors(factors, 5)
    call reset_positions(positions, 5)

    ! B(1,1)
    factors(1) = 1.0d0 / eps
    positions(1, :) = [1, 1]
    ! B(3,3)
    factors(2) = rho
    positions(2, :) = [3, 3]
    ! B(4,4)
    factors(3) = rho / eps
    positions(3, :) = [4, 4]
    ! B(5,5)
    factors(4) = rho / eps
    positions(4, :) = [5, 5]
    ! B(6,6)
    factors(5) = 1.0d0
    positions(5, :) = [6, 6]

    call subblock(quadblock, factors, positions, curr_weight, &
                  h_quadratic, h_quadratic)

    ! Cubic * Cubic
    call reset_factors(factors, 3)
    call reset_positions(positions, 3)

    ! B(2,2)
    factors(1) = rho / eps
    positions(1, :) = [2, 2]
    ! B(7,7)
    factors(2) = 1.0d0 / eps
    positions(2, :) = [7, 7]
    ! B(8,8)
    factors(3) = 1.0d0
    positions(3, :) = [8, 8]

    call subblock(quadblock, factors, positions, curr_weight, h_cubic, h_cubic)

    if (allocated(factors)) then
      deallocate(factors)
    end if
    if (allocated(positions)) then
      deallocate(positions)
    end if

  end subroutine get_B_elements

end module mod_setup_matrix_b
