program esonas
  use mod_global_variables
  implicit none

  complex(dp), allocatable  :: matrix_A(:, :)
  real(dp), allocatable     :: matrix_B(:, :)
  real(dp), allocatable     :: grid(:), grid_gauss(:)


  call initialisation

  call create_matrices

  call cleanup



contains

  !> Initialises the grid and equilibrium configuration
  subroutine initialisation()
    use mod_setup_grid
    use mod_setup_equilibrium
    use mod_radiative_cooling

    ! Initialises global variables
    call initialise_variables

    if (radiative_cooling) then
      call initialise_radiative_cooling
    end if


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
    use mod_setup_matrix_b
    use mod_setup_matrix_a

    call construct_B(grid, grid_gauss, matrix_B)
    call construct_A(grid, grid_gauss, matrix_A)

    return

  end subroutine create_matrices

  !> Performs cleanup, deallocates variables
  subroutine cleanup()
    use mod_setup_grid
    use mod_setup_equilibrium
    use mod_setup_matrix_b
    use mod_radiative_cooling
    deallocate(matrix_A)
    deallocate(matrix_B)
    deallocate(grid)

    call variables_clean
    call grid_clean
    call equilibrium_clean
    call matrix_B_clean

    if (radiative_cooling) then
      call radiative_cooling_clean
    end if



  end subroutine cleanup


end program esonas
