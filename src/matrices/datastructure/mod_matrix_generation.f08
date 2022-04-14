module mod_matrix_generation
  use mod_global_variables, only: dp
  use mod_matrix_structure, only: matrix_t, new_matrix
  implicit none

  private

  interface generate_array_from_matrix
    module procedure generate_real_array_from_matrix
    module procedure generate_complex_array_from_matrix
  end interface generate_array_from_matrix

  public :: generate_matrix_from_array
  public :: generate_array_from_matrix

contains


  !> Converts a given 2D array to the linked-list matrix datastructure.
  function generate_matrix_from_array(array, label) result(matrix)
    !> the array used to generate the matrix datastructure
    class(*), intent(in) :: array(:, :)
    !> optional label for matrix structure
    character(len=*), intent(in), optional :: label
    !> matrix datastructure corresponding to the array
    type(matrix_t) :: matrix
    integer :: irow, icol

    matrix = new_matrix(nb_rows=size(array, dim=1), label=label)
    do icol = 1, size(array, dim=2)
      do irow = 1, size(array, dim=1)
        call matrix%add_element(row=irow, column=icol, element=array(irow, icol))
      end do
    end do
  end function generate_matrix_from_array


  !> Converts a given matrix data structure with real nodes to a 2D real array.
  subroutine generate_real_array_from_matrix(matrix, array)
    !> the matrix datastructure used to generate the array
    type(matrix_t), intent(in) :: matrix
    !> the real 2D array generated from the matrix
    real(dp), intent(out) :: array(matrix%matrix_dim, matrix%matrix_dim)
    integer :: irow, icol

    do icol = 1, matrix%matrix_dim
      do irow = 1, matrix%matrix_dim
        array(irow, icol) = matrix%get_real_element(row=irow, column=icol)
      end do
    end do
  end subroutine generate_real_array_from_matrix


  !> Converts a given matrix data structure with complex nodes to a 2D complex array.
  subroutine generate_complex_array_from_matrix(matrix, array)
    !> the matrix datastructure used to generate the array
    type(matrix_t), intent(in) :: matrix
    !> the complex 2D array generated from the matrix
    complex(dp), intent(out) :: array(matrix%matrix_dim, matrix%matrix_dim)
    integer :: irow, icol

    do icol = 1, matrix%matrix_dim
      do irow = 1, matrix%matrix_dim
        array(irow, icol) = matrix%get_complex_element(row=irow, column=icol)
      end do
    end do
  end subroutine generate_complex_array_from_matrix

end module mod_matrix_generation