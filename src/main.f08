! =============================================================================
!> Main program for the Legolas finite element code.
!! Matrices, eigenvalues and left/right eigenvectors are defined here and passed
!! on to the different modules and submodules.
!!
!! <tt>Legolas</tt> is currently being developed by Niels Claes, Jordi De Jonghe
!! and Rony Keppens, at the Centre for mathematical Plasma-Astrophysics (CmPA),
!! KU Leuven, Belgium.
program legolas
  use mod_global_variables, only: dp, str_len, show_results, dry_run
  use mod_matrix_manager, only: build_matrices
  use mod_solvers, only: solve_evp
  use mod_output, only: datfile_name, create_datfile
  use mod_logging, only: log_message, str, print_console_info, print_whitespace
  use mod_inspections, only: handle_spurious_eigenvalues
  implicit none

  !> A matrix in eigenvalue problem wBX = AX
  complex(dp), allocatable  :: matrix_A(:, :)
  !> B matrix in eigenvalue problem wBX = AX
  real(dp), allocatable     :: matrix_B(:, :)
  !> array with eigenvalues
  complex(dp), allocatable  :: omega(:)
  !> matrix with right eigenvectors, column indices correspond to omega indices
  complex(dp), allocatable  :: eigenvecs_right(:, :)

  call initialisation()
  call print_console_info()

  call build_matrices(matrix_B, matrix_A)

  if (.not. dry_run) then
    call log_message("solving eigenvalue problem...", level="info")
    call solve_evp(matrix_A, matrix_B, omega, eigenvecs_right)
  else
    call log_message( &
      "running dry, overriding parfile and setting eigenvalues to zero", level="info"  &
    )
    omega = (0.0d0, 0.0d0)
  end if

  call handle_spurious_eigenvalues(omega)

  call create_eigenfunctions()
  call create_datfile(omega, matrix_A, matrix_B)

  call cleanup()

  if (show_results) then
    call print_whitespace(1)
    call execute_command_line("python3 pylbo_wrapper.py -i " // trim(datfile_name))
  end if

contains

  !> Subroutine responsible for all initialisations.
  !! Allocates and initialises main and global variables, then the equilibrium state
  !! and eigenfunctions are initialised and the equilibrium is set.
  subroutine initialisation()
    use mod_global_variables, only: initialise_globals, matrix_gridpts, &
      solver, number_of_eigenvalues, write_eigenfunctions, gamma, set_gamma
    use mod_input, only: read_parfile, get_parfile
    use mod_equilibrium, only: initialise_equilibrium, set_equilibrium
    use mod_logging, only: print_logo
    use mod_global_variables, only: viscosity, hall_mhd

    character(len=str_len)  :: parfile
    integer   :: nb_evs

    call initialise_globals()
    call get_parfile(parfile)
    call read_parfile(parfile)
    call set_gamma(gamma)

    call print_logo()

    allocate(matrix_A(matrix_gridpts, matrix_gridpts))
    allocate(matrix_B(matrix_gridpts, matrix_gridpts))

    if (solver == "arnoldi") then
      nb_evs = number_of_eigenvalues
    else
      nb_evs = matrix_gridpts
    end if
    call log_message("setting #eigenvalues to " // str(nb_evs), level="debug")
    allocate(omega(nb_evs))

    call initialise_equilibrium()
    call set_equilibrium()

    ! TODO: remove this warning when fully tested
    if (viscosity) then
      call log_message( &
        "using viscous MHD, note that this is not yet fully tested!", level="warning" &
      )
    end if
    if (hall_mhd) then
      call log_message( &
        "using Hall MHD, note that this does not yet work properly!", level="warning" &
      )
    end if

    ! Arnoldi solver needs this, since it always calculates an orthonormal basis
    if (write_eigenfunctions .or. solver == "arnoldi") then
      call log_message("allocating eigenvector arrays", level="debug")
      ! we need #rows = matrix dimension, #cols = #eigenvalues
      allocate(eigenvecs_right(matrix_gridpts, nb_evs))
    else
      ! @note: this is needed to prevent segfaults, since it seems that in some
      ! cases for macOS the routine zgeev references the right eigenvectors even
      ! if they are not requested.
      call log_message("allocating eigenvector arrays as dummy", level="debug")
      allocate(eigenvecs_right(2, 2))
    end if
  end subroutine initialisation


  !> Initialises and calculates the eigenfunctions if requested.
  subroutine create_eigenfunctions()
    use mod_global_variables, only: write_eigenfunctions
    use mod_eigenfunctions, only: initialise_eigenfunctions, calculate_eigenfunctions

    if (write_eigenfunctions) then
      call initialise_eigenfunctions(omega)
      call calculate_eigenfunctions(eigenvecs_right)
    end if
  end subroutine create_eigenfunctions


  !> Deallocates all main variables, then calls the cleanup
  !! routines of all relevant subroutines to do the same thing.
  subroutine cleanup()
    use mod_global_variables, only: radiative_cooling
    use mod_grid, only: grid_clean
    use mod_equilibrium, only: equilibrium_clean
    use mod_radiative_cooling, only: radiative_cooling_clean
    use mod_eigenfunctions, only: eigenfunctions_clean

    deallocate(matrix_A)
    deallocate(matrix_B)
    deallocate(omega)
    if (allocated(eigenvecs_right)) then
      deallocate(eigenvecs_right)
    end if

    call grid_clean()
    call equilibrium_clean()

    if (radiative_cooling) then
      call radiative_cooling_clean()
    end if
    call eigenfunctions_clean()
  end subroutine cleanup

end program legolas
