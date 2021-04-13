submodule (mod_matrix_manager) smod_regular_matrix
  implicit none

contains

  module subroutine add_bmatrix_terms(gauss_idx, current_weight, quadblock)
    integer, intent(in)   :: gauss_idx
    real(dp), intent(in)  :: current_weight
    complex(dp), intent(inout)  :: quadblock(:, :)

    real(dp)  :: rho, eps

    rho = rho_field % rho0(gauss_idx)
    eps = eps_grid(gauss_idx)

    ! Quadratic * Quadratic
    call reset_factor_positions(new_size=5)
    ! B(1, 1)
    factors(1) = 1.0d0
    positions(1, :) = [1, 1]
    ! B(3, 3)
    factors(2) = eps * rho
    positions(2, :) = [3, 3]
    ! B(4, 4)
    factors(3) = rho
    positions(3, :) = [4, 4]
    ! B(5, 5)
    factors(4) = rho
    positions(4, :) = [5, 5]
    ! B(6, 6)
    factors(5) = eps
    positions(5, :) = [6, 6]
    call subblock(quadblock, factors, positions, current_weight, h_quad, h_quad)

    ! Cubic * Cubic
    call reset_factor_positions(new_size=3)
    ! B(2, 2)
    factors(1) = rho
    positions(1, :) = [2, 2]
    ! B(7, 7)
    factors(2) = 1.0d0
    positions(2, :) = [7, 7]
    ! B(8, 8)
    factors(3) = eps
    positions(3, :) = [8, 8]
    call subblock(quadblock, factors, positions, current_weight, h_cubic, h_cubic)
  end subroutine add_bmatrix_terms


  module subroutine add_regular_matrix_terms(gauss_idx, current_weight, quadblock)
    integer, intent(in)   :: gauss_idx
    real(dp), intent(in)  :: current_weight
    complex(dp), intent(inout)  :: quadblock(:, :)


  end subroutine add_regular_matrix_terms

end submodule smod_regular_matrix