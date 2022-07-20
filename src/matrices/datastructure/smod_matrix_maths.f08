! =============================================================================
!> Contains the implementation of various operator overloading methods between
!! linked-list matrix datastructures.
submodule (mod_matrix_structure) smod_matrix_maths
  implicit none

contains

  !> Overloads the addition operator between two matrix datastructures.
  !! Raises an error if the two matrices are not compatible.
  module procedure add_matrices
    integer :: irow, inode
    type(node_t), pointer :: current_node

    if (matrix1%matrix_dim /= matrix2%matrix_dim) then
      call raise_incompatibility_error(matrix1%matrix_dim, matrix2%matrix_dim)
      return
    end if

    matrix = matrix1%copy()
    do irow = 1, matrix2%matrix_dim
      current_node => matrix2%rows(irow)%head
      do inode = 1, matrix2%rows(irow)%nb_elements
        call matrix%add_element( &
          row=irow, &
          column=current_node%column, &
          element=current_node%get_node_element() &
        )
        current_node => current_node%next
      end do
    end do
    nullify(current_node)
  end procedure add_matrices


  !> Overloads the subtraction operator between two matrix datastructures.
  !! Raises an error if the two matrices are not compatible.
  module procedure subtract_matrices
    integer :: irow, inode
    type(node_t), pointer :: current_node

    if (matrix1%matrix_dim /= matrix2%matrix_dim) then
      call raise_incompatibility_error(matrix1%matrix_dim, matrix2%matrix_dim)
      return
    end if

    matrix = matrix1%copy()
    do irow = 1, matrix2%matrix_dim
      current_node => matrix2%rows(irow)%head
      do inode = 1, matrix2%rows(irow)%nb_elements
        call matrix%add_element( &
          row=irow, &
          column=current_node%column, &
          element = -current_node%get_node_element() &
        )
        current_node => current_node%next
      end do
    end do
  end procedure subtract_matrices


  !> Overloads the multiplication operator between a matrix and a real vector.
  module procedure matrix_x_real_vector
    integer :: irow, inode
    complex(dp) :: array_value
    complex(dp) :: element
    type(node_t), pointer :: current_node

    do irow = 1, matrix%matrix_dim
      array_value = 0.0d0
      current_node => matrix%rows(irow)%head
      do inode = 1, matrix%rows(irow)%nb_elements
        element = current_node%get_node_element()
        array_value = array_value + element * vector(current_node%column)
        current_node => current_node%next
      end do
      array(irow) = array_value
    end do
    nullify(current_node)
  end procedure matrix_x_real_vector


  !> Overloads the multiplication operator between a matrix and a complex vector.
  module procedure matrix_x_complex_vector
    integer :: irow, inode
    complex(dp) :: array_value
    complex(dp) :: element
    type(node_t), pointer :: current_node

    do irow = 1, matrix%matrix_dim
      array_value = (0.0d0, 0.0d0)
      current_node => matrix%rows(irow)%head
      do inode = 1, matrix%rows(irow)%nb_elements
        element = current_node%get_node_element()
        array_value = array_value + element * vector(current_node%column)
        current_node => current_node%next
      end do
      array(irow) = array_value
    end do
    nullify(current_node)
  end procedure matrix_x_complex_vector


  !> Overloads the multiplication operator between a complex number and a matrix
  module procedure matrix_x_number
    integer :: irow, inode
    type(node_t), pointer :: current_node

    matrix_out = matrix%copy()
    do irow = 1, matrix_out%matrix_dim
      current_node => matrix_out%rows(irow)%head
      do inode = 1, matrix_out%rows(irow)%nb_elements
        call current_node%multiply_node_element_with(number)
        current_node => current_node%next
      end do
    end do
    nullify(current_node)
  end procedure matrix_x_number


  !> Raises an error if two square matrices are not compatible.
  subroutine raise_incompatibility_error(dim1, dim2)
    !> dimension of the first matrix
    integer, intent(in) :: dim1
    !> dimension of the second matrix
    integer, intent(in) :: dim2

    call log_message( &
      "incompatible ranks! unable to add matrix of dimension " &
      // str(dim1) // " with matrix of dimension " // str(dim2), &
      level="error" &
    )
  end subroutine raise_incompatibility_error

end submodule smod_matrix_maths
