module mod_matrix_generation
  use mod_global_variables, only: dp
  use mod_matrix_structure, only: matrix_t, new_matrix
  implicit none

  private

  interface generate_array_from_matrix
    module procedure generate_real_array_from_matrix
  end interface generate_array_from_matrix

  public :: generate_matrix_from_array
  public :: generate_array_from_matrix

contains

  function generate_matrix_from_array(array) result(matrix)
    !> the array used to generate the matrix datastructure
    class(*), intent(in) :: array(:, :)
    !> matrix datastructure corresponding to the array
    type(matrix_t) :: matrix
    integer :: irow, icol

    matrix = new_matrix(nb_rows=size(array, dim=1))
    do icol = 1, size(array, dim=2)
      do irow = 1, size(array, dim=1)
        call matrix%add_element(row=irow, column=icol, element=array(irow, icol))
      end do
    end do
  end function generate_matrix_from_array


  function generate_real_array_from_matrix(matrix) result(array)
    !> the matrix datastructure used to generate the array
    type(matrix_t), intent(in) :: matrix
    !> the real 2D array generated from the matrix
    real(dp) :: array(matrix%matrix_dim, matrix%matrix_dim)

    array = 0.0d0
  end function generate_real_array_from_matrix

end module mod_matrix_generation