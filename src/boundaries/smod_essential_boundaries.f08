submodule (mod_boundary_manager) smod_essential_boundaries
  use mod_global_variables, only: dp_LIMIT
  use mod_equilibrium_params, only: k2, k3
  use mod_check_values, only: is_zero

  implicit none

  character(len=3), allocatable :: cubic_vars_to_zero_out(:)

contains

  module procedure apply_essential_boundaries_left
    use mod_get_indices, only: get_subblock_index

    integer :: limits(2)

    ! left side quadblock limits are (1, 1) -> (dim_quadblock, dim_quadblock)
    limits = [1, settings%dims%get_dim_quadblock()]

    ! on the left side we have a zero in a quadratic basis function (number 2), which
    ! zeroes out the odd rows/columns. We explicitly handle this by introducing an
    ! element on the diagonal for these indices as well.
    call zero_out_row_and_col( &
      matrix=matrix, &
      idxs=get_subblock_index( &
        variables=[character(len=3) :: "rho", "v2", "v3", "T", "a1"], &
        state_vector=settings%get_state_vector(), &
        dim_subblock=settings%dims%get_dim_subblock(), &
        odd=.true., &
        edge="left" &
      ), &
      limits=limits &
    )

    ! wall/regularity conditions: v1 has to be zero. Cubic, odd row/col to zero
    cubic_vars_to_zero_out = ["v1"]
    ! wall/regularity conditions: (k3 * a2 - k2 * a3) has to be zero.
    select case(settings%equilibrium%get_boundary_type())
    case("wall")
      ! for "wall" we force a2 = a3 = 0 regardless of k2 and k3
      cubic_vars_to_zero_out = [cubic_vars_to_zero_out, "a2", "a3"]
    case("wall_weak")
      ! for "wall_weak" we only force a2 resp. a3 to zero if k3 resp. k2 is nonzero
      if (.not. is_zero(k2)) cubic_vars_to_zero_out = [cubic_vars_to_zero_out, "a3"]
      if (.not. is_zero(k3)) cubic_vars_to_zero_out = [cubic_vars_to_zero_out, "a2"]
    end select
    ! apply wall/regularity conditions
    call zero_out_row_and_col( &
      matrix=matrix, &
      idxs=get_subblock_index( &
        cubic_vars_to_zero_out, &
        settings%get_state_vector(), &
        settings%dims%get_dim_subblock(), &
        odd=.true., &
        edge="left" &
      ), &
      limits=limits &
    )

    ! if T boundary conditions are needed, set even row/colum (quadratic) to zero
    if (apply_T_bounds) then
      call zero_out_row_and_col( &
        matrix=matrix, &
        idxs=get_subblock_index( &
          ["T"], &
          settings%get_state_vector(), &
          settings%dims%get_dim_subblock(), &
          odd=.false., &
          edge="left" &
        ), &
        limits=limits &
      )
    end if
    ! if no-slip boundary conditions are needed, then v2 and v3 should equal the
    ! wall's tangential velocities, here zero in the even rows/columns (both quadratic)
    if (apply_noslip_bounds_left) then
      call zero_out_row_and_col( &
        matrix=matrix, &
        idxs=get_subblock_index( &
          ["v2", "v3"], &
          settings%get_state_vector(), &
          settings%dims%get_dim_subblock(), &
          odd=.false., &
          edge="left" &
        ), &
        limits=limits &
      )
    end if

    if (allocated(cubic_vars_to_zero_out)) deallocate(cubic_vars_to_zero_out)
  end procedure apply_essential_boundaries_left


  module procedure apply_essential_boundaries_right
    use mod_get_indices, only: get_subblock_index

    integer :: ishift
    integer :: limits(2)

    ! index shift, even number so end of previous quadblock
    ishift = matrix%matrix_dim - settings%dims%get_dim_quadblock()
    ! last quadblock indices hence run from ishift + 1 to matrix dimension
    limits = [ishift + 1, matrix%matrix_dim]

    ! fixed wall: v1 should be zero
    cubic_vars_to_zero_out = ["v1"]
    select case(settings%equilibrium%get_boundary_type())
    case("wall")
      cubic_vars_to_zero_out = [cubic_vars_to_zero_out, "a2", "a3"]
    case("wall_weak")
      ! (k3 * a2 - k2 * a3) should be zero
      if (.not. is_zero(k2)) cubic_vars_to_zero_out = [cubic_vars_to_zero_out, "a3"]
      if (.not. is_zero(k3)) cubic_vars_to_zero_out = [cubic_vars_to_zero_out, "a2"]
    end select
    ! apply wall/regularity conditions
    call zero_out_row_and_col( &
      matrix=matrix, &
      idxs=ishift + get_subblock_index( &
        cubic_vars_to_zero_out, &
        settings%get_state_vector(), &
        settings%dims%get_dim_subblock(), &
        odd=.true., &
        edge="right" &
      ), &
      limits=limits &
    )

    ! T condition
    if (apply_T_bounds) then
      call zero_out_row_and_col( &
        matrix=matrix, &
        idxs=ishift + get_subblock_index( &
          ["T"], &
          settings%get_state_vector(), &
          settings%dims%get_dim_subblock(), &
          odd=.false., &
          edge="right" &
        ), &
        limits=limits &
      )
    end if
    ! no-slip condition
    if (apply_noslip_bounds_right) then
      call zero_out_row_and_col( &
        matrix=matrix, &
        idxs=ishift + get_subblock_index( &
          ["v2", "v3"], &
          settings%get_state_vector(), &
          settings%dims%get_dim_subblock(), &
          odd=.false., &
          edge="right" &
        ), &
        limits=limits &
      )
    end if

    if (allocated(cubic_vars_to_zero_out)) deallocate(cubic_vars_to_zero_out)
  end procedure apply_essential_boundaries_right


  !> Zeroes out the row and column corresponding to the given indices.
  !! Afterwards `diagonal_factor` is introduced in that row/column on the main diagonal.
  subroutine zero_out_row_and_col(matrix, idxs, limits)
    !> the matrix under consideration
    type(matrix_t), intent(inout) :: matrix
    !> indices of the row and column to zero out
    integer, intent(in) :: idxs(:)
    !> (start, end) limits of quadblock corresponding to (start, start):(end, end)
    integer, intent(in) :: limits(2)
    integer :: i, k, idx

    do i = 1, size(idxs)
      idx = idxs(i)
      do k = limits(1), limits(2)
        ! row at idx, columns within limits
        call matrix%rows(idx)%delete_node_from_row(column=k)
        ! column at idx, rows within limits
        call matrix%rows(k)%delete_node_from_row(column=idx)
      end do
      ! add diagonal factor to main diagonal
      call matrix%add_element(row=idx, column=idx, element=get_diagonal_factor(matrix))
    end do
  end subroutine zero_out_row_and_col


  !> Returns the value that is introduced on the main block diagonal after
  !! zeroing out the corresponding row and column.
  !! Depends on the matrix that is used.
  function get_diagonal_factor(matrix) result(diagonal_factor)
    use mod_global_variables, only: NaN

    type(matrix_t), intent(in) :: matrix
    complex(dp) :: diagonal_factor

    if (matrix%get_label() == "B") then
      diagonal_factor = (1.0d0, 0.0d0)
    else if (matrix%get_label() == "A") then
      diagonal_factor = (0.0d0, 0.0d0)
    else
      diagonal_factor = NaN
      call logger%error( &
        "get_diagonal_factor: invalid or empty matrix label: " // matrix%get_label() &
      )
    end if
  end function get_diagonal_factor


end submodule smod_essential_boundaries
