!
! PROGRAM: legolas
!
!> @author: Niels Claes
!> niels.claes@kuleuven.be
!
! DESCRIPTION:
!> Main program. Finite element code to calculate eigenvalues
!! and eigenvectors of the complete non-adiabatic MHD spectrum.
!! Included physics: flow, radiative cooling, thermal conduction, resistivity
!

program legolas
  use mod_global_variables
  use mod_output
  use mod_equilibrium, only: set_equilibrium
  use mod_matrix_creation, only: create_matrices
  use mod_solvers, only: solve_QR
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

  ! allocate variables, initialise grid and physics
  call initialisation()
  ! set the equilibrium configuration
  call set_equilibrium()
  ! print configuration to console
  call startup_info_toconsole()

  ! create matrices A and B
  write(*, *) "Creating matrices..."
  call create_matrices(matrix_B, matrix_A)

  ! solve eigenvalue problem
  write(*, *) "Solving eigenvalue problem..."
  call solve_QR(matrix_A, matrix_B, omega, eigenvecs_left, eigenvecs_right)

  ! write spectrum, eigenvectors, matrices etc. to file
  call save_solutions()

  ! deallocate everything
  call cleanup()

  if (show_results) then
    write(*, *) ""
    write(*, *) "Plotting results..."
    call execute_command_line("python3 python/process_legolas.py")
  end if


contains

  !> Initialises the grid and equilibrium configuration.
  subroutine initialisation()
    use mod_input, only: read_parfile, get_parfile
    use mod_grid, only: initialise_grid
    use mod_equilibrium, only: initialise_equilibrium
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

    ! Initialise equilibrium
    call initialise_equilibrium()
    ! Initialise equilibrium derivatives
    call initialise_equilibrium_derivatives()
    ! Initialise eigenfunction arrays
    call initialise_eigenfunctions()

  end subroutine initialisation


  !> Saves the solutions
  subroutine save_solutions()
    use mod_eigenfunctions
    use mod_output

    write(*, *) "Writing eigenvalues to file..."
    call eigenvalues_tofile(omega, savename_eigenvalues)

    if (write_matrices) then
      write(*, *) "Writing matrices to file..."
      call matrices_tofile(matrix_A, matrix_B, trim(savename_matrix // 'A'), &
                           trim(savename_matrix // 'B'))
    end if

    if (write_eigenvectors) then
      write(*, *) "Writing eigenvectors to file..."
      call eigenvectors_tofile(eigenvecs_left, eigenvecs_right, &
                               trim(savename_eigenvectors // '_left'), &
                               trim(savename_eigenvectors // '_right'))
    end if

    if (write_eigenfunctions) then
      write(*, *) "Writing eigenfunctions to file..."
      call calculate_eigenfunctions(eigenvecs_right)
    end if
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
      call radiative_cooling_clean()
    end if
    call eigenfunctions_clean()

  end subroutine cleanup


end program legolas
