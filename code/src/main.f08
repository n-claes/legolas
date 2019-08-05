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
  use mod_info
  use mod_timer
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

  ! Initialisations
  call set_time_start()
  call initialisation()
  call set_time_initialisation
  call show_startup_info()

  ! Matrix creation
  call set_time_start()
  call create_matrices()
  call set_time_matrices()

  ! Solving eigenvalue problem (timing is done inside module)
  call solve_eigenvalue_problem()

  ! Save solutions (timing is done inside module)
  call save_solutions()
  call show_final_info()

  ! Wrap up
  call cleanup()

  ! Plot solutions if needed
  call plot_solutions()


contains

  !> Initialises the grid and equilibrium configuration.
  subroutine initialisation()
    use mod_input, only: read_parfile, get_parfile
    use mod_grid, only: initialise_grid
    use mod_equilibrium, only: initialise_equilibrium, set_equilibrium
    use mod_equilibrium_derivatives, only: initialise_equilibrium_derivatives
    use mod_eigenfunctions, only: initialise_eigenfunctions

    character(len=str_len)  :: parfile

    call get_parfile(parfile)
    call read_parfile(parfile)

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

    write(*, *) "Creating matrices..."

    call construct_B(matrix_B)
    call construct_A(matrix_A)

  end subroutine create_matrices

  !> Calls the solver.
  subroutine solve_eigenvalue_problem()
    use mod_solvers

    write(*, *) "Solving eigenvalue problem..."

    call solve_QR(matrix_A, matrix_B, omega, eigenvecs_left, eigenvecs_right)

  end subroutine solve_eigenvalue_problem

  !> Saves the solutions
  subroutine save_solutions()
    use mod_eigenfunctions
    use mod_io

    write(*, *)
    write(*, *) "Writing configuration to file..."
    call save_config("config")

    write(*, *) "Writing eigenvalues to file..."
    call set_time_start()
    call save_eigenvalues(omega, "eigenvalues", append=.false., stream=.true.)
    call save_eigenvalues(omega, "eigenvalues_text", append=.false., stream=.false.)
    call set_time_write_omegas()

    if (write_AB) then
      write(*, *) "Writing matrices to file..."
      call set_time_start()
      call save_matrices(matrix_A, matrix_B, "matrix_A", "matrix_B")
      call set_time_write_matrices()
    end if

    if (write_eigenvectors) then
      write(*, *) "Writing eigenvectors to file..."
      call set_time_start()
      call save_eigenvectors(eigenvecs_left, eigenvecs_right, &
                             "v_left", "v_right")
      call set_time_write_eigenvectors()
    end if

    if (write_eigenfunctions) then
      write(*, *) "Writing eigenfunctions to file..."
      call set_time_start()
      call get_all_eigenfunctions(eigenvecs_right)
      call set_time_write_eigenfunctions()
    end if
  end subroutine save_solutions

  !> Routine to call Python script if plotting is requested.
  subroutine plot_solutions()
    if (plot_when_finished) then
      write(*, *) ""
      write(*, *) ""
      write(*, *) "Plotting results..."
      call execute_command_line("python python/process_esonas.py")
    end if
  end subroutine plot_solutions

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
