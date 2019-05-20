program esonas
  implicit none

  complex, allocatable          :: matrix_A(:, :)
  real, allocatable             :: matrix_B(:, :)
  real, allocatable             :: grid(:)


  call initialisation
  call create_matrices


  call cleanup



contains

  subroutine initialisation()
    use mod_global_variables
    use mod_setup_grid
    use mod_setup_equilibrium

    ! Initialises global variables
    call init_variables

    allocate(matrix_A(matrix_gridpts, matrix_gridpts))
    allocate(matrix_B(matrix_gridpts, matrix_gridpts))
    allocate(grid(integral_gridpts))

    ! Initialise grid
    call initialise_grid(grid)

    ! Initialise equilibrium
    call initialise_equilibrium

  end subroutine initialisation

  subroutine create_matrices()
    use mod_global_variables
    use mod_setup_matrix_b

    call construct_B(grid, matrix_B)


    return

  end subroutine create_matrices


  subroutine cleanup()
    use mod_global_variables
    use mod_setup_equilibrium
    deallocate(matrix_A)
    deallocate(matrix_B)
    deallocate(grid)
    deallocate(geometry)

    deallocate(rho_0)
    deallocate(v_0)
    deallocate(T_0)
    deallocate(B_0)

  end subroutine cleanup


end program esonas
