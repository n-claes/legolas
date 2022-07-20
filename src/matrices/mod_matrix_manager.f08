module mod_matrix_manager
  use mod_global_variables, only: dp, ir, ic, gamma_1
  use mod_grid, only: grid, grid_gauss, eps_grid, d_eps_grid_dr
  use mod_make_subblock, only: subblock
  use mod_equilibrium, only: rho_field, T_field, B_field
  use mod_equilibrium_params, only: k2, k3
  use mod_logging, only: log_message, str
  use mod_matrix_shortcuts, only: get_G_operator, get_F_operator, get_wv_operator
  implicit none

  !> quadratic basis functions
  real(dp)  :: h_quad(4)
  !> derivative of quadratic basis functions
  real(dp)  :: dh_quad(4)
  !> cubic basis functions
  real(dp)  :: h_cubic(4)
  !> derivative of cubic basis functions
  real(dp)  :: dh_cubic(4)
  !> array of factors, these are the integrands at every gaussian point
  complex(dp), allocatable  :: factors(:)
  !> array of positions, governs where factors are placed in the matrices
  integer, allocatable  :: positions(:, :)

  interface
    module subroutine add_bmatrix_terms(gauss_idx, current_weight, quadblock)
      integer, intent(in)   :: gauss_idx
      real(dp), intent(in)  :: current_weight
      complex(dp), intent(inout)  :: quadblock(:, :)
    end subroutine add_bmatrix_terms

    module subroutine add_regular_matrix_terms(gauss_idx, current_weight, quadblock)
      integer, intent(in)   :: gauss_idx
      real(dp), intent(in)  :: current_weight
      complex(dp), intent(inout)  :: quadblock(:, :)
    end subroutine add_regular_matrix_terms

    module subroutine add_flow_matrix_terms(gauss_idx, current_weight, quadblock)
      integer, intent(in)   :: gauss_idx
      real(dp), intent(in)  :: current_weight
      complex(dp), intent(inout)  :: quadblock(:, :)
    end subroutine add_flow_matrix_terms

    module subroutine add_resistive_matrix_terms(gauss_idx, current_weight, quadblock)
      integer, intent(in)   :: gauss_idx
      real(dp), intent(in)  :: current_weight
      complex(dp), intent(inout)  :: quadblock(:, :)
    end subroutine add_resistive_matrix_terms

    module subroutine add_cooling_matrix_terms(gauss_idx, current_weight, quadblock)
      integer, intent(in)   :: gauss_idx
      real(dp), intent(in)  :: current_weight
      complex(dp), intent(inout)  :: quadblock(:, :)
    end subroutine add_cooling_matrix_terms

    module subroutine add_conduction_matrix_terms(gauss_idx, current_weight, quadblock)
      integer, intent(in)   :: gauss_idx
      real(dp), intent(in)  :: current_weight
      complex(dp), intent(inout)  :: quadblock(:, :)
    end subroutine add_conduction_matrix_terms

    module subroutine add_viscosity_matrix_terms(gauss_idx, current_weight, quadblock)
      integer, intent(in)   :: gauss_idx
      real(dp), intent(in)  :: current_weight
      complex(dp), intent(inout)  :: quadblock(:, :)
    end subroutine add_viscosity_matrix_terms

    module subroutine add_hall_matrix_terms(gauss_idx, current_weight, quadblock)
      integer, intent(in)   :: gauss_idx
      real(dp), intent(in)  :: current_weight
      complex(dp), intent(inout)  :: quadblock(:, :)
    end subroutine add_hall_matrix_terms

    module subroutine add_hall_bmatrix_terms(gauss_idx, current_weight, quadblock)
      integer, intent(in)   :: gauss_idx
      real(dp), intent(in)  :: current_weight
      complex(dp), intent(inout)  :: quadblock(:, :)
    end subroutine add_hall_bmatrix_terms
  end interface

  private

  public  :: build_matrices
  public  :: reset_factor_positions

contains

  subroutine build_matrices(matrix_B, matrix_A)
    use mod_global_variables, only: gridpts, dim_quadblock, dim_subblock, &
        n_gauss, gaussian_weights, flow, resistivity, radiative_cooling, &
        thermal_conduction, viscosity, hall_mhd
    use mod_spline_functions, only: quadratic_factors, quadratic_factors_deriv, &
      cubic_factors, cubic_factors_deriv
    use mod_matrix_structure, only: matrix_t
    use mod_boundary_manager, only: apply_boundary_conditions

    !> the B-matrix
    type(matrix_t), intent(inout) :: matrix_B
    !> the A-matrix
    type(matrix_t), intent(inout) :: matrix_A

    !> quadblock for the A-matrix
    complex(dp) :: quadblock_A(dim_quadblock, dim_quadblock)
    !> quadblock for the B-matrix
    complex(dp) :: quadblock_B(dim_quadblock, dim_quadblock)
    !> left side of current interval
    real(dp)  :: x_left
    !> right side of current interval
    real(dp)  :: x_right
    !> current position in the Gaussian grid
    real(dp)  :: current_x_gauss
    !> current weight
    real(dp)  :: current_weight

    integer :: i, j, k, l, idx1, idx2
    integer :: quadblock_idx, gauss_idx

    ! used to shift the quadblock along the main diagonal
    quadblock_idx = 0

    do i = 1, gridpts - 1
      ! reset quadblocks
      quadblock_A = (0.0d0, 0.0d0)
      quadblock_B = (0.0d0, 0.0d0)

      ! get interval boundaries
      x_left = grid(i)
      x_right = grid(i + 1)

      ! loop over Gaussian points to calculate integral
      do j = 1, n_gauss
        ! current grid index of Gaussian grid
        gauss_idx = (i - 1) * n_gauss + j
        current_x_gauss = grid_gauss(gauss_idx)
        current_weight = gaussian_weights(j)

        ! calculate spline functions for this point in the Gaussian grid
        call quadratic_factors(current_x_gauss, x_left, x_right, h_quad)
        call quadratic_factors_deriv(current_x_gauss, x_left, x_right, dh_quad)
        call cubic_factors(current_x_gauss, x_left, x_right, h_cubic)
        call cubic_factors_deriv(current_x_gauss, x_left, x_right, dh_cubic)

        ! get matrix elements
        call add_bmatrix_terms(gauss_idx, current_weight, quadblock_B)
        call add_regular_matrix_terms(gauss_idx, current_weight, quadblock_A)
        if (flow) then
          call add_flow_matrix_terms(gauss_idx, current_weight, quadblock_A)
        end if
        if (resistivity) then
          call add_resistive_matrix_terms(gauss_idx, current_weight, quadblock_A)
        end if
        if (radiative_cooling) then
          call add_cooling_matrix_terms(gauss_idx, current_weight, quadblock_A)
        end if
        if (thermal_conduction) then
          call add_conduction_matrix_terms(gauss_idx, current_weight, quadblock_A)
        end if
        if (viscosity) then
          call add_viscosity_matrix_terms(gauss_idx, current_weight, quadblock_A)
        end if
        if (hall_mhd) then
          call add_hall_matrix_terms(gauss_idx, current_weight, quadblock_A)
          call add_hall_bmatrix_terms(gauss_idx, current_weight, quadblock_B)
        end if
      end do

      ! dx from integral
      quadblock_B = quadblock_B * (x_right - x_left)
      quadblock_A = quadblock_A * (x_right - x_left)

      ! fill matrices
      !> @note  The quadblock is shifted along the main (tri)diagonal.
      !!        We add `dim_subblock` instead of `dim_quadblock` to the indices,
      !!        due to overlapping of the bottom-right part of the quadblock with the
      !!        top-left part of the next grid interval.
      do k = 1, dim_quadblock
        do l = 1, dim_quadblock
          idx1 = k + quadblock_idx
          idx2 = l + quadblock_idx
          call matrix_B%add_element( &
            row=idx1, column=idx2, element=real(quadblock_B(k, l)) &
          )
          call matrix_A%add_element(row=idx1, column=idx2, element=quadblock_A(k, l))
        end do
      end do
      quadblock_idx = quadblock_idx + dim_subblock
    end do

    call apply_boundary_conditions(matrix_A, matrix_B)
  end subroutine build_matrices


  !> Resets the <tt>factors</tt> and <tt>positions</tt> arrays to a given
  !! new size.
  subroutine reset_factor_positions(new_size)
    !> the new size for the <tt>factors</tt> and <tt>positions</tt> arrays
    integer, intent(in) :: new_size

    if (allocated(factors)) then
      deallocate(factors)
    end if
    allocate(factors(new_size))

    if (allocated(positions)) then
      deallocate(positions)
    end if
    allocate(positions(new_size, 2))
  end subroutine reset_factor_positions

end module mod_matrix_manager
