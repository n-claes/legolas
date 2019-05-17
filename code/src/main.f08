program esonas
  implicit none

  complex, allocatable          :: matrix_A(:, :)
  real, allocatable             :: matrix_B(:, :)


  call initialisation
  call create_matrices

  ! First find a way to set global parameters?
  ! similar to use mod_global_parameters in amrvac

  call cleanup



contains

  subroutine initialisation()
    use mod_global_variables
    call init_variables

    allocate(matrix_A(matrix_gridpts, matrix_gridpts))
    allocate(matrix_B(matrix_gridpts, matrix_gridpts))

  end subroutine initialisation

  subroutine create_matrices()
    use mod_global_variables
    use mod_setup_matrix_b

    call construct_B(matrix_B)


    return

  end subroutine create_matrices


  subroutine cleanup()
    use mod_global_variables
    deallocate(matrix_A)
    deallocate(matrix_B)
    deallocate(geometry)

  end subroutine cleanup


end program esonas
