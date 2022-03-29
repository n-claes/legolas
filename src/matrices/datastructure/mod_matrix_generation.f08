module mod_matrix_generation
  use mod_global_variables, only: dp
  use mod_matrix_structure, only: matrix_t, new_matrix
  implicit none

  private

  public :: generate_matrix_from_array

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

end module mod_matrix_generation