program esonas
  implicit none

  complex, allocatable            :: matrix_A(:, :)
  real, allocatable               :: matrix_B(:, :)
  double precision, allocatable   :: grid(:), grid_gauss(:)


  call initialisation
  call create_matrices


  call cleanup



contains

  !> Initialises the grid and equilibrium configuration
  subroutine initialisation()
    use mod_global_variables
    use mod_setup_grid
    use mod_setup_equilibrium

    ! Initialises global variables
    call init_variables

    allocate(matrix_A(matrix_gridpts, matrix_gridpts))
    allocate(matrix_B(matrix_gridpts, matrix_gridpts))
    allocate(grid(gridpts))
    allocate(grid_gauss(4*gridpts))

    ! Initialise grid
    call initialise_grid(grid)

    ! Initialise equilibrium
    call initialise_equilibrium(grid, grid_gauss)

  end subroutine initialisation

  !> Creates A and B matrices for the wBX = AX eigenvalue problem
  subroutine create_matrices()
    use mod_global_variables
    use mod_setup_matrix_b

    call construct_B(grid, matrix_B)


    return

  end subroutine create_matrices

  !> Performs cleanup, deallocates variables
  subroutine cleanup()
    use mod_global_variables
    use mod_setup_grid
    use mod_setup_equilibrium
    deallocate(matrix_A)
    deallocate(matrix_B)
    deallocate(grid)

    call variables_clean
    call grid_clean
    call equilibrium_clean

  end subroutine cleanup


end program esonas
