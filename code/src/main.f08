program esonas
  use mod_global_variables
  implicit none

  !> A matrix in eigenvalue problem wBX = AX
  complex(dp), allocatable  :: matrix_A(:, :)
  !> B matrix in eigenvalue problem wBX = AX
  real(dp), allocatable     :: matrix_B(:, :)

  call initialisation

  call create_matrices

  call cleanup



contains

  !> Initialises the grid and equilibrium configuration
  subroutine initialisation()
    use mod_grid
    use mod_equilibrium
    use mod_equilibrium_derivatives


    ! Initialises global variables
    call initialise_variables()

    allocate(matrix_A(matrix_gridpts, matrix_gridpts))
    allocate(matrix_B(matrix_gridpts, matrix_gridpts))

    ! Initialise grid
    call initialise_grid()

    ! Initialise equilibrium
    call initialise_equilibrium()

    ! Initialise equilibrium derivatives
    call initialise_equilibrium_derivatives()

  end subroutine initialisation

  !> Creates A and B matrices for the wBX = AX eigenvalue problem
  subroutine create_matrices()
    use mod_setup_matrix_b
    use mod_setup_matrix_a

    call construct_B(matrix_B)
    call construct_A(matrix_A)

    return

  end subroutine create_matrices

  !> Performs cleanup, deallocates variables
  subroutine cleanup()
    use mod_grid
    use mod_equilibrium
    use mod_equilibrium_derivatives
    use mod_setup_matrix_b
    use mod_setup_matrix_a
    use mod_radiative_cooling
    deallocate(matrix_A)
    deallocate(matrix_B)

    call variables_clean
    call grid_clean
    call equilibrium_clean
    call equilibrium_derivatives_clean
    call matrix_B_clean
    call matrix_A_clean

    if (radiative_cooling) then
      call radiative_cooling_clean
    end if



  end subroutine cleanup


end program esonas
