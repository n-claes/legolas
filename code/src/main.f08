!
! PROGRAM: esonas
!
!> @author: Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Main program. Finite element code to calculate eigenvalues
!! and eigenvectors of the complete non-adiabatic MHD spectrum.
!! Included physics: flow, radiative cooling, thermal conduction, resistivity
!

program esonas
  use mod_global_variables
  implicit none

  !> A matrix in eigenvalue problem wBX = AX
  complex(dp), allocatable  :: matrix_A(:, :)
  !> B matrix in eigenvalue problem wBX = AX
  real(dp), allocatable     :: matrix_B(:, :)
  !> Solutions to the eigenvalue problem
  complex(dp), allocatable  :: omega(:)
  !> Right eigenvectors
  complex(dp), allocatable  :: eigenvecs_right(:, :)
  !> Left eigenvectors
  complex(dp), allocatable  :: eigenvecs_left(:, :)

  call initialisation()

  call create_matrices()

  call solve_eigenvalue_problem()

  call get_eigenfunctions()

  call save_solutions()

  call cleanup()

  write(*, *) "Program finished."



contains

  !> Initialises the grid and equilibrium configuration.
  subroutine initialisation()
    use mod_input, only: read_parfile
    use mod_grid, only: initialise_grid
    use mod_equilibrium, only: initialise_equilibrium, set_equilibrium
    use mod_equilibrium_derivatives, only: initialise_equilibrium_derivatives
    use mod_eigenfunctions, only: initialise_eigenfunctions

    call read_parfile()

    ! Allocate matrices
    allocate(matrix_A(matrix_gridpts, matrix_gridpts))
    allocate(matrix_B(matrix_gridpts, matrix_gridpts))
    allocate(omega(matrix_gridpts))
    allocate(eigenvecs_right(matrix_gridpts, matrix_gridpts))
    allocate(eigenvecs_left(matrix_gridpts, matrix_gridpts))

    ! Initialise grid
    call initialise_grid()
    ! Initialise equilibrium
    call initialise_equilibrium()
    ! Initialise equilibrium derivatives
    call initialise_equilibrium_derivatives()
    ! Initialise eigenfunction arrays
    call initialise_eigenfunctions()

    ! Set equilibrium
    call set_equilibrium()

  end subroutine initialisation

  !> Creates A and B matrices for the wBX = AX eigenvalue problem.
  subroutine create_matrices()
    use mod_setup_matrix_b, only: construct_B
    use mod_setup_matrix_a, only: construct_A

    call construct_B(matrix_B)
    call construct_A(matrix_A)

  end subroutine create_matrices

  !> Calls the solver.
  subroutine solve_eigenvalue_problem()
    use mod_solvers

    call solve_QR(matrix_A, matrix_B, omega, eigenvecs_left, eigenvecs_right)

  end subroutine solve_eigenvalue_problem

  !> Calculates the eigenfunctions
  subroutine get_eigenfunctions()
    use mod_eigenfunctions

    call get_all_eigenfunctions(omega, eigenvecs_right)

  end subroutine get_eigenfunctions

  !> Saves the solutions
  subroutine save_solutions()
    use mod_io

    write(*, *) "Writing eigenvalues to file..."
    call save_eigenvalues(omega, "eigenvalues", append=.false., stream=.true.)
    write(*, *) "Writing configuration to file..."
    call save_config("config")
    if (write_AB) then
      write(*, *) "Writing matrices to file..."
      call save_matrices(matrix_A, matrix_B, "matrix_A", "matrix_B")
    end if
    if (write_eigenvectors) then
      write(*, *) "Writing eigenvectors to file..."
      call save_eigenvectors(eigenvecs_left, eigenvecs_right, &
                             "v_left", "v_right")
    end if
    call plot_results()
  end subroutine save_solutions

  !> Performs cleanup, deallocates variables.
  subroutine cleanup()
    use mod_grid, only: grid_clean
    use mod_equilibrium, only: equilibrium_clean
    use mod_equilibrium_derivatives, only: equilibrium_derivatives_clean
    use mod_radiative_cooling, only: radiative_cooling_clean
    use mod_eigenfunctions, only: eigenfunctions_clean

    deallocate(matrix_A)
    deallocate(matrix_B)
    deallocate(omega)
    deallocate(eigenvecs_left)
    deallocate(eigenvecs_right)

    call grid_clean()
    call equilibrium_clean()
    call equilibrium_derivatives_clean()

    if (radiative_cooling) then
      call radiative_cooling_clean
    end if
    call eigenfunctions_clean()



  end subroutine cleanup


end program esonas
