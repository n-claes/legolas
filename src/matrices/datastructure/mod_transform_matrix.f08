! =============================================================================
!> Contains various subroutines and functions to switch between linked-list
!! matrix representations, banded matrix representations, and full array matrices.
module mod_transform_matrix
  use mod_global_variables, only: dp, NaN
  use mod_matrix_structure, only: matrix_t, new_matrix
  use mod_banded_matrix, only: banded_matrix_t, new_banded_matrix
  implicit none

  private

  interface matrix_to_array
    module procedure matrix_to_complex_array
  end interface matrix_to_array

  interface matrix_to_banded
    module procedure matrix_to_complex_banded
  end interface matrix_to_banded

  interface array_to_matrix
    module procedure general_array_to_matrix
  end interface array_to_matrix

  interface array_to_banded
    module procedure array_to_complex_banded
  end interface array_to_banded

  interface banded_to_array
    module procedure banded_to_complex_array
  end interface banded_to_array

  public :: matrix_to_array
  public :: matrix_to_banded

  public :: array_to_matrix
  public :: array_to_banded

  public :: banded_to_array


contains


  !> Converts a given matrix data structure with complex nodes to a 2D complex array.
  subroutine matrix_to_complex_array(matrix, array)
    !> the original matrix datastructure
    type(matrix_t), intent(in) :: matrix
    !> the resulting complex 2D array
    complex(dp), intent(out) :: array(matrix%matrix_dim, matrix%matrix_dim)
    integer :: irow, icol

    do icol = 1, matrix%matrix_dim
      do irow = 1, matrix%matrix_dim
        array(irow, icol) = matrix%get_complex_element(row=irow, column=icol)
      end do
    end do
  end subroutine matrix_to_complex_array


  !> Converts a matrix data structure into a complex banded matrix.
  subroutine matrix_to_complex_banded(matrix, subdiags, superdiags, banded)
    !> the original matrix datastructure
    type(matrix_t), intent(in) :: matrix
    !> number of subdiagonals
    integer, intent(in) :: subdiags
    !> number of superdiagonals
    integer, intent(in) :: superdiags
    !> the resulting banded datastructure
    type(banded_matrix_t), intent(out) :: banded
    integer :: irow, icol

    banded = new_banded_matrix( &
      rows=matrix%matrix_dim, &
      cols=matrix%matrix_dim, &
      subdiags=subdiags, &
      superdiags=superdiags &
    )
    do icol = 1, matrix%matrix_dim
      do irow = max(1, icol - superdiags), min(matrix%matrix_dim, icol + subdiags)
        call banded%set_element( &
          row=irow, &
          col=icol, &
          element=matrix%get_complex_element(row=irow, column=icol) &
        )
      end do
    end do
  end subroutine matrix_to_complex_banded


  !> Converts a given 2D array to the linked-list matrix datastructure.
  function general_array_to_matrix(array, label) result(matrix)
    !> the original array
    class(*), intent(in) :: array(:, :)
    !> optional label for matrix datastructure
    character(len=*), intent(in), optional :: label
    !> the resulting matrix datastructure
    type(matrix_t) :: matrix
    integer :: irow, icol

    matrix = new_matrix(nb_rows=size(array, dim=1), label=label)
    do icol = 1, size(array, dim=2)
      do irow = 1, size(array, dim=1)
        call matrix%add_element(row=irow, column=icol, element=array(irow, icol))
      end do
    end do
  end function general_array_to_matrix


  !> Converts a given real array to a banded datastructure.
  subroutine array_to_complex_banded(array, subdiags, superdiags, banded)
    !> the original array
    class(*), intent(in) :: array(:, :)
    !> the number of subdiagonals
    integer, intent(in) :: subdiags
    !> the number of superdiagonals
    integer, intent(in) :: superdiags
    !> the resulting banded datastructure
    type(banded_matrix_t), intent(out) :: banded
    integer :: nrows, ncols, irow, icol

    nrows = size(array, dim=1)
    ncols = size(array, dim=2)
    banded = new_banded_matrix( &
      rows=nrows, cols=ncols, subdiags=subdiags, superdiags=superdiags &
    )
    do icol = 1, ncols
      do irow = max(1, icol - superdiags), min(nrows, icol + subdiags)
        call banded%set_element( &
          row=irow, &
          col=icol, &
          element=get_array_element(array, irow, icol) &
        )
      end do
    end do
  end subroutine array_to_complex_banded


  !> Converts a banded datastructure to a full complex array.
  pure function banded_to_complex_array(banded) result(array)
    !> the original banded datastructure
    type(banded_matrix_t), intent(in) :: banded
    !> the resulting complex array
    complex(dp) :: array(banded%m, banded%n)
    integer :: irow, icol

    do icol = 1, banded%n
      do irow = 1, banded%m
        array(irow, icol) = banded%get_element(row=irow, col=icol)
      end do
    end do
  end function banded_to_complex_array


  !> Retrieves the element at index (i, j) for an array of general type.
  !! Returns the element as a (casted) complex type.
  pure function get_array_element(array, irow, icol) result(element)
    !> the general array
    class(*), intent(in) :: array(:, :)
    !> row index of element
    integer, intent(in) :: irow
    !> column index of element
    integer, intent(in) :: icol
    !> the element at position (irow, icol), cast to complex
    complex(dp) :: element

    select type(array)
      type is (complex(dp))
        element = array(irow, icol)
      type is (real(dp))
        element = cmplx(array(irow, icol), kind=dp)
      class default
        element = cmplx(NaN, NaN, kind=dp)
    end select
  end function get_array_element

end module mod_transform_matrix
