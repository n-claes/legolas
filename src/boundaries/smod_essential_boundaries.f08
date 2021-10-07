submodule (mod_boundary_manager) smod_essential_boundaries
  implicit none

contains

  module procedure apply_essential_boundaries_left
    complex(dp) :: diagonal_factor

    diagonal_factor = get_diagonal_factor(matrix)

    ! on the left side we have a zero in a quadratic basis function (number 2), which
    ! zeroes out the odd rows/columns. We explicitly handle this by introducing an
    ! element on the diagonal for these indices as well. The odd numbered indices for
    ! the quadratic elements correspond to rho, v2, v2, T, a1 = 1, 5, 6, 9, 11.
    call set_row_col_to_value(quadblock, val=diagonal_factor, idxs=[1, 5, 7, 9, 11])
    ! wall or regularity conditions: v1, a2 and a3 have to be zero. These are cubic
    ! elements, so we force non-zero basis functions (odd rows & columns) to zero.
    call set_row_col_to_value(quadblock, val=diagonal_factor, idxs=[3, 13, 15])
    ! if T boundary conditions are needed, set even row/colum (quadratic)
    if (apply_T_bounds) then
      call set_row_col_to_value(quadblock, val=diagonal_factor, idxs=[10])
    end if
    ! if no-slip boundary conditions are needed, then v2 and v3 should equal the
    ! wall's tangential velocities, here zero in the even rows/columns (both quadratic)
    if (apply_noslip_bounds_left) then
      call set_row_col_to_value(quadblock, val=diagonal_factor, idxs=[6, 8])
    end if
  end procedure apply_essential_boundaries_left


  module procedure apply_essential_boundaries_right
    complex(dp) :: diagonal_factor

    diagonal_factor = get_diagonal_factor(matrix)

    ! for a fixed wall v1, a2 and a3 should be zero, same reasoning as for left side.
    call set_row_col_to_value(quadblock, val=diagonal_factor, idxs=[19, 29, 31])
    ! T condition
    if (apply_T_bounds) then
      call set_row_col_to_value(quadblock, val=diagonal_factor, idxs=[26])
    end if
    ! no-slip condition
    if (apply_noslip_bounds_right) then
      call set_row_col_to_value(quadblock, val=diagonal_factor, idxs=[22, 24])
    end if
  end procedure apply_essential_boundaries_right


  subroutine set_row_col_to_value(quadblock, val, idxs)
    complex(dp), intent(inout) :: quadblock(:, :)
    complex(dp), intent(in) :: val
    integer, intent(in) :: idxs(:)
    integer :: i, j

    do i = 1, size(idxs)
      j = idxs(i)
      quadblock(j, :) = (0.0d0, 0.0d0)
      quadblock(:, j) = (0.0d0, 0.0d0)
      quadblock(j, j) = val
    end do
  end subroutine set_row_col_to_value


  function get_diagonal_factor(matrix) result(diagonal_factor)
    use mod_global_variables, only: solver, NaN

    character, intent(in) :: matrix
    complex(dp) :: diagonal_factor

    if (matrix == "B") then
      diagonal_factor = (1.0d0, 0.0d0)
    else if (matrix == "A") then
      if (solver == "arnoldi") then
        diagonal_factor = (0.0d0, 0.0d0)
      else
        diagonal_factor = (1.0d20, 0.0d0)
      end if
    else
      diagonal_factor = NaN
      call log_message( &
        "get_diagonal_factor: invalid matrix argument " // matrix, level="error" &
      )
    end if
  end function get_diagonal_factor


end submodule smod_essential_boundaries