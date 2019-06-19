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
    use mod_grid, only: initialise_grid
    use mod_equilibrium, only: initialise_equilibrium
    use mod_equilibrium_derivatives, only: initialise_equilibrium_derivatives


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

    !! TODO: checks for interpolating cooling curves beyond table values!
    !! TODO: normalise values AFTER getting derivatives

  end subroutine initialisation

  !> Creates A and B matrices for the wBX = AX eigenvalue problem
  subroutine create_matrices()
    use mod_setup_matrix_b, only: construct_B
    use mod_setup_matrix_a, only: construct_A
    use mod_solvers

    call construct_B(matrix_B)
    call construct_A(matrix_A)

    call solve_QR(matrix_A, matrix_B)

    return

  end subroutine create_matrices

  !> Performs cleanup, deallocates variables
  subroutine cleanup()
    use mod_grid, only: grid_clean
    use mod_equilibrium, only: equilibrium_clean
    use mod_equilibrium_derivatives, only: equilibrium_derivatives_clean
    use mod_setup_matrix_b, only: matrix_B_clean
    use mod_setup_matrix_a, only: matrix_A_clean
    use mod_radiative_cooling, only: radiative_cooling_clean

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
