submodule (mod_boundary_manager) smod_essential_boundaries
  implicit none

contains

  module procedure apply_essential_boundaries_left
    use mod_get_indices, only: get_subblock_index

    ! on the left side we have a zero in a quadratic basis function (number 2), which
    ! zeroes out the odd rows/columns. We explicitly handle this by introducing an
    ! element on the diagonal for these indices as well.
    call zero_out_row_and_col( &
      quadblock=quadblock, &
      diagonal_factor=get_diagonal_factor(matrix), &
      idxs=get_subblock_index( &
        [character(len=5) :: "rho", "v2", "v3", "T", "a1"], odd=.true., edge="left" &
      ) &
    )
    ! wall or regularity conditions: v1, a2 and a3 have to be zero. These are cubic
    ! elements, so we force non-zero basis functions (odd rows & columns) to zero.
    call zero_out_row_and_col( &
      quadblock=quadblock, &
      diagonal_factor=get_diagonal_factor(matrix), &
      idxs=get_subblock_index(["v1", "a2", "a3"], odd=.true., edge="left") &
    )
    ! if T boundary conditions are needed, set even row/colum (quadratic)
    if (apply_T_bounds) then
      call zero_out_row_and_col( &
        quadblock=quadblock, &
        diagonal_factor=get_diagonal_factor(matrix), &
        idxs=get_subblock_index(["T"], odd=.false., edge="left") &
      )
    end if
    ! if no-slip boundary conditions are needed, then v2 and v3 should equal the
    ! wall's tangential velocities, here zero in the even rows/columns (both quadratic)
    if (apply_noslip_bounds_left) then
      call zero_out_row_and_col( &
        quadblock=quadblock, &
        diagonal_factor=get_diagonal_factor(matrix), &
        idxs=get_subblock_index(["v2", "v3"], odd=.false., edge="left") &
      )
    end if
  end procedure apply_essential_boundaries_left


  module procedure apply_essential_boundaries_right
    use mod_get_indices, only: get_subblock_index

    ! for a fixed wall v1, a2 and a3 should be zero, same reasoning as for left side.
    call zero_out_row_and_col( &
      quadblock=quadblock, &
      diagonal_factor=get_diagonal_factor(matrix), &
      idxs=get_subblock_index(["v1", "a2", "a3"], odd=.true., edge="right") &
    )
    ! T condition
    if (apply_T_bounds) then
      call zero_out_row_and_col( &
        quadblock=quadblock, &
        diagonal_factor=get_diagonal_factor(matrix), &
        idxs=get_subblock_index(["T"], odd=.false., edge="right") &
      )
    end if
    ! no-slip condition
    if (apply_noslip_bounds_right) then
      call zero_out_row_and_col( &
        quadblock=quadblock, &
        diagonal_factor=get_diagonal_factor(matrix), &
        idxs=get_subblock_index(["v2", "v3"], odd=.false., edge="right") &
      )
    end if
  end procedure apply_essential_boundaries_right


  !> Zeroes out the row and column corresponding to the given indices.
  !! Afterwards `diagonal_factor` is introduced in that row/column on the main diagonal.
  subroutine zero_out_row_and_col(quadblock, diagonal_factor, idxs)
    !> quadblock at the left or right edge, modified on exit
    complex(dp), intent(inout)  :: quadblock(:, :)
    !> value to introduce on the main diagonal
    complex(dp), intent(in) :: diagonal_factor
    !> indices of the row and column to zero out
    integer, intent(in) :: idxs(:)
    integer :: i, j

    do i = 1, size(idxs)
      j = idxs(i)
      quadblock(j, :) = (0.0d0, 0.0d0)
      quadblock(:, j) = (0.0d0, 0.0d0)
      quadblock(j, j) = diagonal_factor
    end do
  end subroutine zero_out_row_and_col


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