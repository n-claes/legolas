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
        if (.not. associated(current_node)) exit
        call matrix%add_element( &
          row=irow, column=current_node%column, element=current_node%element &
        )
        current_node => current_node%next
      end do
      current_node => matrix2%rows(irow)%head
      do inode = 1, matrix2%rows(irow)%nb_elements
        if (.not. associated(current_node)) exit
        call matrix%add_element( &
          row=irow, column=current_node%column, element=current_node%element &
        )
        current_node => current_node%next
      end do
    end do
    nullify(current_node)
  end procedure add_matrices

end submodule mod_matrix_maths