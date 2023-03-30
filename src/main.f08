! =============================================================================
!> Main program for the Legolas finite element code.
!! Matrices, eigenvalues and left/right eigenvectors are defined here and passed
!! on to the different modules and submodules.
!!
!! <tt>Legolas</tt> is currently being developed by Niels Claes, Jordi De Jonghe
!! and Rony Keppens, at the Centre for mathematical Plasma-Astrophysics (CmPA),
!! KU Leuven, Belgium.
program legolas
  use mod_global_variables, only: dp, str_len, initialise_globals
  use mod_matrix_structure, only: matrix_t
  use mod_equilibrium, only: set_equilibrium
  use mod_matrix_manager, only: build_matrices
  use mod_solvers, only: solve_evp
  use mod_output, only: datfile_path, create_datfile
  use mod_logging, only: logger, str
  use mod_console, only: print_console_info, print_whitespace
  use mod_timing, only: timer_t, new_timer
  use mod_settings, only: settings_t, new_settings
  use mod_background, only: background_t, new_background
  use mod_eigenfunctions, only: eigenfunctions_t, new_eigenfunctions
  use mod_physics, only: physics_t, new_physics
  implicit none

  !> A matrix in eigenvalue problem wBX = AX
  type(matrix_t) :: matrix_A
  !> B matrix in eigenvalue problem wBX = AX
  type(matrix_t) :: matrix_B
  !> timer used by the whole program
  type(timer_t) :: timer
  !> dedicated settings type
  type(settings_t) :: settings
  type(background_t) :: background
  type(eigenfunctions_t) :: eigenfunctions
  type(physics_t) :: physics
  !> array with eigenvalues
  complex(dp), allocatable  :: omega(:)
  !> matrix with right eigenvectors, column indices correspond to omega indices
  complex(dp), allocatable  :: right_eigenvectors(:, :)

  call initialise_globals()
  call logger%initialise()

  timer = new_timer()
  settings = new_settings()
  background = new_background()
  physics = new_physics(settings, background)

  call timer%start_timer()
  call initialisation()
  call set_equilibrium(settings, background, physics)
  timer%init_time = timer%end_timer()

  call print_console_info(settings)

  call timer%start_timer()
  call build_matrices(matrix_B, matrix_A, settings, background, physics)
  timer%matrix_time = timer%end_timer()

  call logger%info("solving eigenvalue problem...")
  call timer%start_timer()
  call solve_evp(matrix_A, matrix_B, settings, omega, right_eigenvectors)
  timer%evp_time = timer%end_timer()

  call timer%start_timer()
  eigenfunctions = new_eigenfunctions(settings, background)
  call get_eigenfunctions()
  timer%eigenfunction_time = timer%end_timer()

  call timer%start_timer()
  call create_datfile( &
    settings, &
    background, &
    physics, &
    omega, &
    matrix_A, &
    matrix_B, &
    right_eigenvectors, &
    eigenfunctions &
  )
  timer%datfile_time = timer%end_timer()

  call cleanup()

  call print_timelog()

  if (settings%io%show_results) then
    call print_whitespace(1)
    call execute_command_line("python3 pylbo_wrapper.py -i " // trim(datfile_path))
  end if

contains

  !> Subroutine responsible for all initialisations.
  !! Allocates and initialises main and global variables, then the equilibrium state
  !! and eigenfunctions are initialised and the equilibrium is set.
  subroutine initialisation()
    use mod_matrix_structure, only: new_matrix
    use mod_input, only: read_parfile, get_parfile
    use mod_console, only: print_logo

    character(len=5*str_len)  :: parfile
    integer   :: nb_evs

    call get_parfile(parfile)
    call read_parfile(parfile, settings)

    call print_logo()
    call logger%info("the physics type is " // settings%get_physics_type())
    call logger%info("the state vector is " // str(settings%get_state_vector()))

    select case(settings%solvers%get_solver())
    case ("arnoldi")
      nb_evs = settings%solvers%number_of_eigenvalues
    case ("inverse-iteration")
      nb_evs = 1
    case default
      nb_evs = settings%dims%get_dim_matrix()
    end select
    call logger%debug("setting #eigenvalues to " // str(nb_evs))
    allocate(omega(nb_evs))
    matrix_A = new_matrix(nb_rows=settings%dims%get_dim_matrix(), label="A")
    matrix_B = new_matrix(nb_rows=settings%dims%get_dim_matrix(), label="B")

    ! Arnoldi solver needs this, since it always calculates an orthonormal basis
    if ( &
      settings%io%should_compute_eigenvectors() &
      .or. settings%solvers%get_solver() == "arnoldi" &
    ) then
      call logger%debug("allocating eigenvector arrays")
      ! we need #rows = matrix dimension, #cols = #eigenvalues
      allocate(right_eigenvectors(settings%dims%get_dim_matrix(), nb_evs))
    else
      ! @note: this is needed to prevent segfaults, since it seems that in some
      ! cases for macOS the routine zgeev references the right eigenvectors even
      ! if they are not requested.
      call logger%debug("allocating eigenvector arrays as dummy")
      allocate(right_eigenvectors(2, 2))
    end if
  end subroutine initialisation


  !> Initialises and calculates the eigenfunctions if requested.
  subroutine get_eigenfunctions()
    if (.not. settings%io%write_eigenfunctions) return
    call eigenfunctions%initialise(omega)
    call eigenfunctions%assemble(right_eigenvectors)
  end subroutine get_eigenfunctions


  !> Deallocates all main variables, then calls the cleanup
  !! routines of all relevant subroutines to do the same thing.
  subroutine cleanup()
    use mod_grid, only: grid_clean

    call matrix_A%delete_matrix()
    call matrix_B%delete_matrix()
    deallocate(omega)
    if (allocated(right_eigenvectors)) deallocate(right_eigenvectors)

    call grid_clean()

    call settings%delete()
    call background%delete()
    call physics%delete()
    call eigenfunctions%delete()
  end subroutine cleanup


  subroutine print_timelog()
    real(dp) :: total_time

    call print_whitespace(1)
    call logger%info("---------------------------------------------")
    call logger%disable_prefix()
    total_time = timer%get_total_time()

    call logger%info("                << Time log >>")
    call logger%info("Legolas finished in " // str(total_time) // " seconds")
    call logger%info("   initialisation: " // str(timer%init_time) // " sec")
    call logger%info("   matrix construction: " // str(timer%matrix_time) // " sec")
    call logger%info("   eigenvalue problem: " // str(timer%evp_time) // " sec")
    call logger%info("   eigenfunctions: " // str(timer%eigenfunction_time) // " sec")
    call logger%info("   datfile creation: " // str(timer%datfile_time) // " sec")
    call logger%enable_prefix()
  end subroutine print_timelog

end program legolas
