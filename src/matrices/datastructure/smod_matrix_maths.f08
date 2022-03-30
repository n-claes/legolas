submodule (mod_matrix_structure) mod_matrix_maths
  implicit none

contains

  module procedure add_matrices
    integer :: irow, inode
    type(node_t), pointer :: current_node

    matrix = new_matrix(nb_rows=matrix1%matrix_dim)

    do irow = 1, matrix1%matrix_dim
      current_node => matrix1%rows(irow)%head
      do inode = 1, matrix1%rows(irow)%nb_elements
        call matrix%add_element( &
          row=irow, column=current_node%column, element=current_node%element &
        )
        current_node => current_node%next
      end do
      current_node => matrix2%rows(irow)%head
      do inode = 1, matrix2%rows(irow)%nb_elements
        call matrix%add_element( &
          row=irow, column=current_node%column, element=current_node%element &
        )
        current_node => current_node%next
      end do
    end do
    nullify(current_node)
  end procedure add_matrices


  module procedure matrix_x_real_vector
    integer :: irow, inode
    real(dp) :: array_value
    real(dp) :: element
    type(node_t), pointer :: current_node

    do irow = 1, matrix%matrix_dim
      array_value = 0.0d0
      current_node => matrix%rows(irow)%head
      do inode = 1, matrix%rows(irow)%nb_elements
        call current_node%get_node_element(element)
        array_value = array_value + element * vector(current_node%column)
        current_node => current_node%next
      end do
      array(irow) = array_value
    end do
    nullify(current_node)
  end procedure matrix_x_real_vector


  module procedure matrix_x_complex_vector
    integer :: irow, inode
    complex(dp) :: array_value
    complex(dp) :: element
    type(node_t), pointer :: current_node

    do irow = 1, matrix%matrix_dim
      array_value = (0.0d0, 0.0d0)
      current_node => matrix%rows(irow)%head
      do inode = 1, matrix%rows(irow)%nb_elements
        call current_node%get_node_element(element)
        array_value = array_value + element * vector(current_node%column)
        current_node => current_node%next
      end do
      array(irow) = array_value
    end do
    nullify(current_node)
  end procedure matrix_x_complex_vector

end submodule mod_matrix_maths
